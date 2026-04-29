#include "mex.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

/* -------------------------------------------------------------------------
 * Pair record for dedupe
 * ------------------------------------------------------------------------- */
typedef struct {
    mwIndex i;     /* 1-based local index, canonical i < j */
    mwIndex j;     /* 1-based local index */
    double r2;
    double dr[3];  /* oriented from i -> j */
} PairRec;

/* -------------------------------------------------------------------------
 * Utilities
 * ------------------------------------------------------------------------- */

static mwSize linearize3(mwSize ix, mwSize iy, mwSize iz,
                         mwSize nx, mwSize ny)
{
    return ix + iy * nx + iz * nx * ny;
}

static void unlinearize3(mwSize lin0, mwSize nx, mwSize ny,
                         mwSize *ix, mwSize *iy, mwSize *iz)
{
    mwSize nxy = nx * ny;
    *iz = lin0 / nxy;
    {
        mwSize rem = lin0 - (*iz) * nxy;
        *iy = rem / nx;
        *ix = rem - (*iy) * nx;
    }
}

static mwSize wrap_index_int(int i, mwSize n)
{
    int nn = (int)n;
    int w = i % nn;
    if (w < 0) w += nn;
    return (mwSize)w;
}

static double round_away_from_zero(double x)
{
    if (x >= 0.0) {
        return floor(x + 0.5);
    } else {
        return ceil(x - 0.5);
    }
}

static double get_thole_a(const double *aPtr, mwSize nA, mwIndex i1)
{
    if (nA == 1) return aPtr[0];
    return aPtr[i1 - 1];
}

static void thole_f3f5_scalar(double r, double alpha_i, double alpha_j, double a,
                              double *f3, double *f5)
{
    if (r <= 0.0 || alpha_i <= 0.0 || alpha_j <= 0.0 || a <= 0.0) {
        *f3 = 1.0;
        *f5 = 1.0;
        return;
    }

    {
        double aij = pow(alpha_i * alpha_j, 1.0 / 6.0);
        double s = a * pow(r / aij, 3.0);
        double e = exp(-s);
        *f3 = 1.0 - e;
        *f5 = 1.0 - (1.0 + s) * e;
    }
}

static int pairrec_cmp_lex(const void *a, const void *b)
{
    const PairRec *pa = (const PairRec *)a;
    const PairRec *pb = (const PairRec *)b;

    if (pa->i < pb->i) return -1;
    if (pa->i > pb->i) return  1;
    if (pa->j < pb->j) return -1;
    if (pa->j > pb->j) return  1;
    return 0;
}

static void ensure_pair_capacity(PairRec **pairs, mwSize *cap, mwSize need)
{
    if (need <= *cap) return;

    {
        mwSize newCap = (*cap == 0) ? 1024 : *cap;
        while (newCap < need) {
            newCap *= 2;
        }

        PairRec *newBuf = (PairRec *)mxCalloc(newCap, sizeof(PairRec));
        if (*pairs != NULL && *cap > 0) {
            memcpy(newBuf, *pairs, (*cap) * sizeof(PairRec));
            mxFree(*pairs);
        }
        *pairs = newBuf;
        *cap = newCap;
    }
}

static void frac_mic_to_cart(const double *frac, mwSize nActive,
                             mwIndex i1, mwIndex j1,
                             const double *H,
                             double *dr)
{
    /* i1, j1 are 1-based local indices */
    double df0, df1, df2;
    const double h00 = H[0], h10 = H[1], h20 = H[2];
    const double h01 = H[3], h11 = H[4], h21 = H[5];
    const double h02 = H[6], h12 = H[7], h22 = H[8];

    df0 = frac[(j1 - 1) + 0 * nActive] - frac[(i1 - 1) + 0 * nActive];
    df1 = frac[(j1 - 1) + 1 * nActive] - frac[(i1 - 1) + 1 * nActive];
    df2 = frac[(j1 - 1) + 2 * nActive] - frac[(i1 - 1) + 2 * nActive];

    df0 -= round_away_from_zero(df0);
    df1 -= round_away_from_zero(df1);
    df2 -= round_away_from_zero(df2);

    /* H is column-major 3x3 with lattice vectors as columns */
    dr[0] = h00 * df0 + h01 * df1 + h02 * df2;
    dr[1] = h10 * df0 + h11 * df1 + h12 * df2;
    dr[2] = h20 * df0 + h21 * df1 + h22 * df2;
}

/* -------------------------------------------------------------------------
 * MEX entry
 *
 * Inputs:
 *   0  fracAct           nActive x 3 double
 *   1  gridShape         1x3 or 3x1 double
 *   2  binHead           nBins x 1 double (1-based linked-list heads)
 *   3  binNext           nActive x 1 double (1-based linked-list next)
 *   4  neighborOffsets   nOffsets x 3 double
 *   5  cellMat           3 x 3 double, lattice vectors as COLUMNS
 *   6  rcut              scalar double
 *   7  alphaEwald        scalar double
 *   8  alphaAct          [] or nActive x 1 double
 *   9  tholeA            [] or scalar or nActive x 1 double
 *
 * Outputs:
 *   0  row_ptr           (nActive+1) x 1 double
 *   1  col_idx           nDir x 1 double
 *   2  dr                nDir x 3 double
 *   3  r2_bare           nDir x 1 double
 *   4  r_bare            nDir x 1 double
 *   5  coeff_iso         nDir x 1 double
 *   6  coeff_dyad        nDir x 1 double
 *   7  nPairsUndirected  scalar double
 * ------------------------------------------------------------------------- */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    const mxArray *fracArr, *gridArr, *headArr, *nextArr, *offsArr;
    const mxArray *cellArr, *cutArr, *alphaEwaldArr, *alphaArr, *aArr;

    const double *frac;
    const double *gridD;
    const double *binHead;
    const double *binNext;
    const double *offs;
    const double *H;
    const double *alphaAct = NULL;
    const double *tholeA = NULL;

    double alphaEwald, alpha2, twoAlphaOverSqrtPi;
    mwSize nActive, nBins, nOffsets;
    mwSize nx, ny, nz;
    double rcut, cutoff2;
    bool haveThole = false;
    mwSize nA = 0;

    PairRec *pairs = NULL;
    mwSize pairsCap = 0;
    mwSize nPairs = 0;
    mwSize nPairsUnique = 0;

    double *rowPtrOut, *colIdxOut, *drOut, *r2Out, *rOut, *coeffIsoOut, *coeffDyadOut;
    double *rowCounts = NULL;
    mwSize *nextPtr = NULL;

    mwSize binLin0, oo;
    mwIndex *wrappedBinForOffset = NULL;
    mwSize zeroOffsetIndex = (mwSize)(-1);
    mwSize nNonzeroOffsets = 0;
    mwSize *nonzeroOffsetIdx = NULL;

    if (nrhs != 10) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:nrhs",
                          "Expected 10 inputs.");
    }
    if (nlhs != 8) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:nlhs",
                          "Expected 8 outputs.");
    }

    fracArr        = prhs[0];
    gridArr        = prhs[1];
    headArr        = prhs[2];
    nextArr        = prhs[3];
    offsArr        = prhs[4];
    cellArr        = prhs[5];
    cutArr         = prhs[6];
    alphaEwaldArr  = prhs[7];
    alphaArr       = prhs[8];
    aArr           = prhs[9];

    if (!mxIsDouble(fracArr) || mxIsComplex(fracArr) || mxGetN(fracArr) != 3) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:frac",
                          "fracAct must be real double nActive x 3.");
    }
    if (!mxIsDouble(gridArr) || mxGetNumberOfElements(gridArr) != 3) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:grid",
                          "gridShape must have 3 elements.");
    }
    if (!mxIsDouble(headArr) || !mxIsDouble(nextArr) ||
        !mxIsDouble(offsArr) || !mxIsDouble(cellArr) ||
        !mxIsDouble(cutArr) || !mxIsDouble(alphaEwaldArr)) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:type",
                          "Inputs 1:8 must be double.");
    }
    if (mxGetM(cellArr) != 3 || mxGetN(cellArr) != 3) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:cell",
                          "cellMat must be real double 3x3 with columns as lattice vectors.");
    }

    frac = mxGetPr(fracArr);
    nActive = mxGetM(fracArr);

    gridD = mxGetPr(gridArr);
    nx = (mwSize)gridD[0];
    ny = (mwSize)gridD[1];
    nz = (mwSize)gridD[2];

    binHead = mxGetPr(headArr);
    binNext = mxGetPr(nextArr);
    offs = mxGetPr(offsArr);
    H = mxGetPr(cellArr);

    nBins = mxGetNumberOfElements(headArr);
    nOffsets = mxGetM(offsArr);

    rcut = mxGetScalar(cutArr);
    if (!(mxIsFinite(rcut) && rcut > 0.0)) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:rcut",
                          "rcut must be positive and finite.");
    }
    cutoff2 = rcut * rcut;

    alphaEwald = mxGetScalar(alphaEwaldArr);
    if (!(mxIsFinite(alphaEwald) && alphaEwald > 0.0)) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:alphaEwald",
                          "alphaEwald must be positive and finite.");
    }
    alpha2 = alphaEwald * alphaEwald;
    twoAlphaOverSqrtPi = 2.0 * alphaEwald / sqrt(M_PI);

    if (!mxIsEmpty(alphaArr) && !mxIsEmpty(aArr)) {
        if (!mxIsDouble(alphaArr) || !mxIsDouble(aArr)) {
            mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:tholeType",
                              "alphaAct and tholeA must be double when provided.");
        }
        if (mxGetNumberOfElements(alphaArr) != nActive) {
            mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:alphaSize",
                              "alphaAct must have length nActive.");
        }
        nA = mxGetNumberOfElements(aArr);
        if (!(nA == 1 || nA == nActive)) {
            mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:aSize",
                              "tholeA must be scalar or length nActive.");
        }
        haveThole = true;
        alphaAct = mxGetPr(alphaArr);
        tholeA = mxGetPr(aArr);
    }

    /* Precompute offset partition */
    nonzeroOffsetIdx = (mwSize *)mxCalloc(nOffsets, sizeof(mwSize));
    for (oo = 0; oo < nOffsets; ++oo) {
        int offx = (int)offs[oo + 0 * nOffsets];
        int offy = (int)offs[oo + 1 * nOffsets];
        int offz = (int)offs[oo + 2 * nOffsets];
        if (offx == 0 && offy == 0 && offz == 0) {
            zeroOffsetIndex = oo;
        } else {
            nonzeroOffsetIdx[nNonzeroOffsets++] = oo;
        }
    }

    /* Precompute wrapped target bins for each (bin, offset) */
    wrappedBinForOffset = (mwIndex *)mxCalloc(nBins * nOffsets, sizeof(mwIndex));
    for (binLin0 = 0; binLin0 < nBins; ++binLin0) {
        mwSize ix, iy, iz;
        unlinearize3(binLin0, nx, ny, &ix, &iy, &iz);

        for (oo = 0; oo < nOffsets; ++oo) {
            int offx = (int)offs[oo + 0 * nOffsets];
            int offy = (int)offs[oo + 1 * nOffsets];
            int offz = (int)offs[oo + 2 * nOffsets];

            mwSize jx = wrap_index_int((int)ix + offx, nx);
            mwSize jy = wrap_index_int((int)iy + offy, ny);
            mwSize jz = wrap_index_int((int)iz + offz, nz);

            wrappedBinForOffset[binLin0 + oo * nBins] =
                (mwIndex)(linearize3(jx, jy, jz, nx, ny) + 1); /* store 1-based */
        }
    }

    /* Larger initial capacity */
    pairsCap = (mwSize)(nActive * 64);
    if (pairsCap < 1024) pairsCap = 1024;
    pairs = (PairRec *)mxCalloc(pairsCap, sizeof(PairRec));

    /* ---------------------------------------------------------------------
     * Pass 1A: true zero-offset triangular same-bin traversal
     * --------------------------------------------------------------------- */
    if (zeroOffsetIndex != (mwSize)(-1)) {
        for (binLin0 = 0; binLin0 < nBins; ++binLin0) {
            mwIndex iHead = (mwIndex)binHead[binLin0];
            if (iHead == 0) continue;

            {
                mwIndex i = iHead;
                while (i != 0) {
                    mwIndex j = (mwIndex)binNext[i - 1];
                    while (j != 0) {
                        double dr[3], r2;

                        frac_mic_to_cart(frac, nActive, i, j, H, dr);
                        r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

                        if (r2 <= cutoff2) {
                            ensure_pair_capacity(&pairs, &pairsCap, nPairs + 1);
                            pairs[nPairs].i = i;
                            pairs[nPairs].j = j;
                            pairs[nPairs].r2 = r2;
                            pairs[nPairs].dr[0] = dr[0];
                            pairs[nPairs].dr[1] = dr[1];
                            pairs[nPairs].dr[2] = dr[2];
                            nPairs += 1;
                        }

                        j = (mwIndex)binNext[j - 1];
                    }
                    i = (mwIndex)binNext[i - 1];
                }
            }
        }
    }

    /* ---------------------------------------------------------------------
     * Pass 1B: all nonzero offsets, full traversal then canonicalize
     * --------------------------------------------------------------------- */
    for (mwSize kk = 0; kk < nNonzeroOffsets; ++kk) {
        oo = nonzeroOffsetIdx[kk];

        for (binLin0 = 0; binLin0 < nBins; ++binLin0) {
            mwIndex iHead = (mwIndex)binHead[binLin0];
            mwIndex jHead;

            if (iHead == 0) continue;

            jHead = (mwIndex)binHead[wrappedBinForOffset[binLin0 + oo * nBins] - 1];
            if (jHead == 0) continue;

            {
                mwIndex i = iHead;
                while (i != 0) {
                    mwIndex j = jHead;
                    while (j != 0) {
                        double dr[3], r2;
                        mwIndex ii, jj;
                        double drCanon[3];

                        if (i == j) {
                            j = (mwIndex)binNext[j - 1];
                            continue;
                        }

                        frac_mic_to_cart(frac, nActive, i, j, H, dr);
                        r2 = dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2];

                        if (r2 <= cutoff2) {
                            if (i < j) {
                                ii = i; jj = j;
                                drCanon[0] = dr[0];
                                drCanon[1] = dr[1];
                                drCanon[2] = dr[2];
                            } else {
                                ii = j; jj = i;
                                drCanon[0] = -dr[0];
                                drCanon[1] = -dr[1];
                                drCanon[2] = -dr[2];
                            }

                            ensure_pair_capacity(&pairs, &pairsCap, nPairs + 1);
                            pairs[nPairs].i = ii;
                            pairs[nPairs].j = jj;
                            pairs[nPairs].r2 = r2;
                            pairs[nPairs].dr[0] = drCanon[0];
                            pairs[nPairs].dr[1] = drCanon[1];
                            pairs[nPairs].dr[2] = drCanon[2];
                            nPairs += 1;
                        }

                        j = (mwIndex)binNext[j - 1];
                    }
                    i = (mwIndex)binNext[i - 1];
                }
            }
        }
    }

    /* Lexicographic sort + adjacent dedupe */
    if (nPairs > 0) {
        qsort(pairs, nPairs, sizeof(PairRec), pairrec_cmp_lex);

        {
            mwSize src = 1;
            mwSize dst = 1;

            while (src < nPairs) {
                if (pairs[src].i != pairs[dst - 1].i ||
                    pairs[src].j != pairs[dst - 1].j) {
                    pairs[dst] = pairs[src];
                    dst += 1;
                }
                src += 1;
            }
            nPairsUnique = dst;
        }
    } else {
        nPairsUnique = 0;
    }

    /* ---------------------------------------------------------------------
     * Build row_ptr
     * --------------------------------------------------------------------- */
    rowCounts = (double *)mxCalloc(nActive, sizeof(double));
    for (mwSize p = 0; p < nPairsUnique; ++p) {
        rowCounts[pairs[p].i - 1] += 1.0;
        rowCounts[pairs[p].j - 1] += 1.0;
    }

    plhs[0] = mxCreateDoubleMatrix(nActive + 1, 1, mxREAL);
    rowPtrOut = mxGetPr(plhs[0]);
    rowPtrOut[0] = 1.0;

    for (mwSize i = 0; i < nActive; ++i) {
        rowPtrOut[i + 1] = rowPtrOut[i] + rowCounts[i];
    }

    /* ---------------------------------------------------------------------
     * Allocate outputs
     * --------------------------------------------------------------------- */
    {
        mwSize nDir = 2 * nPairsUnique;

        plhs[1] = mxCreateDoubleMatrix(nDir, 1, mxREAL);
        plhs[2] = mxCreateDoubleMatrix(nDir, 3, mxREAL);
        plhs[3] = mxCreateDoubleMatrix(nDir, 1, mxREAL);
        plhs[4] = mxCreateDoubleMatrix(nDir, 1, mxREAL);
        plhs[5] = mxCreateDoubleMatrix(nDir, 1, mxREAL);
        plhs[6] = mxCreateDoubleMatrix(nDir, 1, mxREAL);
        plhs[7] = mxCreateDoubleScalar((double)nPairsUnique);

        colIdxOut    = mxGetPr(plhs[1]);
        drOut        = mxGetPr(plhs[2]);
        r2Out        = mxGetPr(plhs[3]);
        rOut         = mxGetPr(plhs[4]);
        coeffIsoOut  = mxGetPr(plhs[5]);
        coeffDyadOut = mxGetPr(plhs[6]);

        nextPtr = (mwSize *)mxCalloc(nActive, sizeof(mwSize));
        for (mwSize i = 0; i < nActive; ++i) {
            nextPtr[i] = (mwSize)(rowPtrOut[i] - 1.0);
        }

        /* Fill CSR from unique pair list; reuse stored geometry */
        for (mwSize p = 0; p < nPairsUnique; ++p) {
            mwIndex i1 = pairs[p].i;
            mwIndex j1 = pairs[p].j;
            double r2 = pairs[p].r2;
            double dr[3];
            double r;
            double invR2, invR, invR3, invR5, invR4;
            double erfcar, expar2, B, C;
            double coeffIso, coeffDyad;
            double alpha_i, alpha_j;
            mwSize k;

            dr[0] = pairs[p].dr[0];
            dr[1] = pairs[p].dr[1];
            dr[2] = pairs[p].dr[2];

            r = sqrt(r2);
            invR2 = 1.0 / r2;
            invR = 1.0 / r;
            invR3 = invR * invR2;
            invR5 = invR3 * invR2;
            invR4 = invR2 * invR2;

            erfcar = erfc(alphaEwald * r);
            expar2 = exp(-alpha2 * r2);

            B = erfcar * invR3 + twoAlphaOverSqrtPi * expar2 * invR2;
            C = 3.0 * erfcar * invR5 +
                twoAlphaOverSqrtPi * (2.0 * alpha2 * invR2 + 3.0 * invR4) * expar2;

            coeffIso = -B;
            coeffDyad = +C;

            if (haveThole) {
                double f3, f5, l3, l5;
                alpha_i = alphaAct[i1 - 1];
                alpha_j = alphaAct[j1 - 1];
                thole_f3f5_scalar(r, alpha_i, alpha_j,
                                  get_thole_a(tholeA, nA, i1),
                                  &f3, &f5);
                l3 = f3 - 1.0;
                l5 = f5 - 1.0;
                coeffIso -= l3 * invR3;
                coeffDyad += 3.0 * l5 * invR5;
            }

            /* i -> j */
            k = nextPtr[i1 - 1]++;
            colIdxOut[k] = (double)j1;
            drOut[k + 0 * nDir] = dr[0];
            drOut[k + 1 * nDir] = dr[1];
            drOut[k + 2 * nDir] = dr[2];
            r2Out[k] = r2;
            rOut[k] = r;
            coeffIsoOut[k] = coeffIso;
            coeffDyadOut[k] = coeffDyad;

            /* j -> i */
            k = nextPtr[j1 - 1]++;
            colIdxOut[k] = (double)i1;
            drOut[k + 0 * nDir] = -dr[0];
            drOut[k + 1 * nDir] = -dr[1];
            drOut[k + 2 * nDir] = -dr[2];
            r2Out[k] = r2;
            rOut[k] = r;
            coeffIsoOut[k] = coeffIso;
            coeffDyadOut[k] = coeffDyad;
        }
    }

    mxFree(wrappedBinForOffset);
    mxFree(nonzeroOffsetIdx);
    mxFree(pairs);
    mxFree(rowCounts);
    mxFree(nextPtr);
}