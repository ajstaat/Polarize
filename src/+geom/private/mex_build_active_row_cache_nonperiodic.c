#include "mex.h"
#include <math.h>
#include <stdbool.h>

/* -------------------------------------------------------------------------
 * Utilities
 * ------------------------------------------------------------------------- */

static mwSize linearize3(mwSize ix, mwSize iy, mwSize iz,
                         mwSize nx, mwSize ny) {
    return ix + iy * nx + iz * nx * ny;
}

static void unlinearize3(mwSize lin, mwSize nx, mwSize ny,
                         mwSize *ix, mwSize *iy, mwSize *iz) {
    mwSize nxy = nx * ny;
    *iz = lin / nxy;
    {
        mwSize rem = lin - (*iz) * nxy;
        *iy = rem / nx;
        *ix = rem - (*iy) * nx;
    }
}

static double get_thole_a(const double *aPtr, mwSize nA, mwIndex i0) {
    if (nA == 1) {
        return aPtr[0];
    }
    return aPtr[i0];
}

static void thole_f3f5_scalar(double r, double alpha_i, double alpha_j, double a,
                              double *f3, double *f5) {
    /* Match MATLAB helper behavior:
       if any of r, alpha_i, alpha_j, a are zero/nonpositive at an entry,
       return undamped limit f3=f5=1. */
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

/* -------------------------------------------------------------------------
 * MEX entry
 *
 * Inputs:
 *   0  posAct            nActive x 3 double
 *   1  gridShape         1x3 or 3x1 double
 *   2  binHead           nBins x 1 double (1-based linked-list heads)
 *   3  binNext           nActive x 1 double (1-based linked-list next)
 *   4  neighborOffsets   nOffsets x 3 double
 *   5  rcut              scalar double
 *   6  alphaAct          [] or nActive x 1 double
 *   7  tholeA            [] or scalar or nActive x 1 double
 *
 * Outputs:
 *   0  row_ptr           (nActive+1) x 1 double
 *   1  col_idx           nDir x 1 double
 *   2  dr                nDir x 3 double
 *   3  r2_bare           nDir x 1 double
 *   4  r_bare            nDir x 1 double
 *   5  inv_r3_bare       nDir x 1 double
 *   6  inv_r5_bare       nDir x 1 double
 *   7  thole_f3          nDir x 1 double or []
 *   8  thole_f5          nDir x 1 double or []
 *   9  nPairsUndirected  scalar double
 * ------------------------------------------------------------------------- */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    const mxArray *posArr, *gridArr, *headArr, *nextArr, *offsArr, *cutArr;
    const mxArray *alphaArr, *aArr;

    const double *pos;
    const double *gridD;
    const double *binHead;
    const double *binNext;
    const double *offs;
    const double *alphaAct = NULL;
    const double *tholeA = NULL;

    mwSize nActive, nBins, nOffsets;
    mwSize nx, ny, nz;
    double rcut, cutoff2;
    bool haveThole = false;
    mwSize nA = 0;

    double *rowPtrOut, *colIdxOut, *drOut, *r2Out, *rOut, *invR3Out, *invR5Out;
    double *f3Out = NULL, *f5Out = NULL;
    double *rowCounts = NULL;
    mwSize *nextPtr = NULL;
    double nPairsUndirected = 0.0;
    double nDirDouble = 0.0;

    mwSize binLin0, oo;

    if (nrhs != 8) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_nonperiodic:nrhs",
                          "Expected 8 inputs.");
    }
    if (nlhs != 10) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_nonperiodic:nlhs",
                          "Expected 10 outputs.");
    }

    posArr   = prhs[0];
    gridArr  = prhs[1];
    headArr  = prhs[2];
    nextArr  = prhs[3];
    offsArr  = prhs[4];
    cutArr   = prhs[5];
    alphaArr = prhs[6];
    aArr     = prhs[7];

    if (!mxIsDouble(posArr) || mxIsComplex(posArr) || mxGetN(posArr) != 3) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_nonperiodic:pos",
                          "posAct must be real double nActive x 3.");
    }
    if (!mxIsDouble(gridArr) || mxGetNumberOfElements(gridArr) != 3) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_nonperiodic:grid",
                          "gridShape must have 3 elements.");
    }
    if (!mxIsDouble(headArr) || !mxIsDouble(nextArr) ||
        !mxIsDouble(offsArr) || !mxIsDouble(cutArr)) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_nonperiodic:type",
                          "Inputs 1:6 must be double.");
    }

    pos = mxGetPr(posArr);
    nActive = mxGetM(posArr);

    gridD = mxGetPr(gridArr);
    nx = (mwSize)gridD[0];
    ny = (mwSize)gridD[1];
    nz = (mwSize)gridD[2];

    binHead = mxGetPr(headArr);
    binNext = mxGetPr(nextArr);
    offs = mxGetPr(offsArr);

    nBins = mxGetNumberOfElements(headArr);
    nOffsets = mxGetM(offsArr);

    rcut = mxGetScalar(cutArr);
    if (!(mxIsFinite(rcut) && rcut > 0.0)) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_nonperiodic:rcut",
                          "rcut must be positive and finite.");
    }
    cutoff2 = rcut * rcut;

    if (!mxIsEmpty(alphaArr) && !mxIsEmpty(aArr)) {
        if (!mxIsDouble(alphaArr) || !mxIsDouble(aArr)) {
            mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_nonperiodic:tholeType",
                              "alphaAct and tholeA must be double when provided.");
        }
        if (mxGetNumberOfElements(alphaArr) != nActive) {
            mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_nonperiodic:alphaSize",
                              "alphaAct must have length nActive.");
        }
        nA = mxGetNumberOfElements(aArr);
        if (!(nA == 1 || nA == nActive)) {
            mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_nonperiodic:aSize",
                              "tholeA must be scalar or length nActive.");
        }
        haveThole = true;
        alphaAct = mxGetPr(alphaArr);
        tholeA = mxGetPr(aArr);
    }

    /* ---------------------------------------------------------------------
     * Pass 1: count rows and undirected pairs
     * --------------------------------------------------------------------- */
    rowCounts = (double *)mxCalloc(nActive, sizeof(double));

    for (binLin0 = 0; binLin0 < nBins; ++binLin0) {
        mwIndex iHead1 = (mwIndex)binHead[binLin0];
        mwSize ix, iy, iz;

        if (iHead1 == 0) continue;

        unlinearize3(binLin0, nx, ny, &ix, &iy, &iz);

        for (oo = 0; oo < nOffsets; ++oo) {
            int dx = (int)offs[oo + 0 * nOffsets];
            int dy = (int)offs[oo + 1 * nOffsets];
            int dz = (int)offs[oo + 2 * nOffsets];
            int jx = (int)ix + dx;
            int jy = (int)iy + dy;
            int jz = (int)iz + dz;
            mwSize jBin0;
            mwIndex jHead1;

            if (jx < 0 || jx >= (int)nx ||
                jy < 0 || jy >= (int)ny ||
                jz < 0 || jz >= (int)nz) {
                continue;
            }

            jBin0 = linearize3((mwSize)jx, (mwSize)jy, (mwSize)jz, nx, ny);
            jHead1 = (mwIndex)binHead[jBin0];
            if (jHead1 == 0) continue;

            if (jBin0 == binLin0) {
                mwIndex i1 = iHead1;
                while (i1 != 0) {
                    mwIndex j1 = (mwIndex)binNext[i1 - 1];
                    while (j1 != 0) {
                        mwSize i0 = (mwSize)(i1 - 1);
                        mwSize j0 = (mwSize)(j1 - 1);
                        double dxp = pos[j0 + 0 * nActive] - pos[i0 + 0 * nActive];
                        double dyp = pos[j0 + 1 * nActive] - pos[i0 + 1 * nActive];
                        double dzp = pos[j0 + 2 * nActive] - pos[i0 + 2 * nActive];
                        double r2 = dxp * dxp + dyp * dyp + dzp * dzp;

                        if (r2 <= cutoff2) {
                            rowCounts[i0] += 1.0;
                            rowCounts[j0] += 1.0;
                            nPairsUndirected += 1.0;
                        }
                        j1 = (mwIndex)binNext[j1 - 1];
                    }
                    i1 = (mwIndex)binNext[i1 - 1];
                }
            } else {
                mwIndex i1 = iHead1;
                while (i1 != 0) {
                    mwIndex j1 = jHead1;
                    while (j1 != 0) {
                        mwSize i0 = (mwSize)(i1 - 1);
                        mwSize j0 = (mwSize)(j1 - 1);
                        double dxp = pos[j0 + 0 * nActive] - pos[i0 + 0 * nActive];
                        double dyp = pos[j0 + 1 * nActive] - pos[i0 + 1 * nActive];
                        double dzp = pos[j0 + 2 * nActive] - pos[i0 + 2 * nActive];
                        double r2 = dxp * dxp + dyp * dyp + dzp * dzp;

                        if (r2 <= cutoff2) {
                            rowCounts[i0] += 1.0;
                            rowCounts[j0] += 1.0;
                            nPairsUndirected += 1.0;
                        }
                        j1 = (mwIndex)binNext[j1 - 1];
                    }
                    i1 = (mwIndex)binNext[i1 - 1];
                }
            }
        }
    }

    /* ---------------------------------------------------------------------
     * row_ptr
     * --------------------------------------------------------------------- */
    plhs[0] = mxCreateDoubleMatrix(nActive + 1, 1, mxREAL);
    rowPtrOut = mxGetPr(plhs[0]);
    rowPtrOut[0] = 1.0;

    {
        mwSize i;
        for (i = 0; i < nActive; ++i) {
            rowPtrOut[i + 1] = rowPtrOut[i] + rowCounts[i];
        }
    }

    nDirDouble = rowPtrOut[nActive] - 1.0;

    /* ---------------------------------------------------------------------
     * Allocate outputs
     * --------------------------------------------------------------------- */
    plhs[1] = mxCreateDoubleMatrix((mwSize)nDirDouble, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix((mwSize)nDirDouble, 3, mxREAL);
    plhs[3] = mxCreateDoubleMatrix((mwSize)nDirDouble, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix((mwSize)nDirDouble, 1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix((mwSize)nDirDouble, 1, mxREAL);
    plhs[6] = mxCreateDoubleMatrix((mwSize)nDirDouble, 1, mxREAL);

    colIdxOut = mxGetPr(plhs[1]);
    drOut     = mxGetPr(plhs[2]);
    r2Out     = mxGetPr(plhs[3]);
    rOut      = mxGetPr(plhs[4]);
    invR3Out  = mxGetPr(plhs[5]);
    invR5Out  = mxGetPr(plhs[6]);

    if (haveThole) {
        plhs[7] = mxCreateDoubleMatrix((mwSize)nDirDouble, 1, mxREAL);
        plhs[8] = mxCreateDoubleMatrix((mwSize)nDirDouble, 1, mxREAL);
        f3Out = mxGetPr(plhs[7]);
        f5Out = mxGetPr(plhs[8]);
    } else {
        plhs[7] = mxCreateDoubleMatrix(0, 0, mxREAL);
        plhs[8] = mxCreateDoubleMatrix(0, 0, mxREAL);
    }

    plhs[9] = mxCreateDoubleScalar(nPairsUndirected);

    /* nextPtr as 0-based write cursors into output arrays */
    nextPtr = (mwSize *)mxCalloc(nActive, sizeof(mwSize));
    {
        mwSize i;
        for (i = 0; i < nActive; ++i) {
            nextPtr[i] = (mwSize)(rowPtrOut[i] - 1.0); /* convert 1-based to 0-based */
        }
    }

    /* ---------------------------------------------------------------------
     * Pass 2: fill outputs directly
     * --------------------------------------------------------------------- */
    for (binLin0 = 0; binLin0 < nBins; ++binLin0) {
        mwIndex iHead1 = (mwIndex)binHead[binLin0];
        mwSize ix, iy, iz;

        if (iHead1 == 0) continue;

        unlinearize3(binLin0, nx, ny, &ix, &iy, &iz);

        for (oo = 0; oo < nOffsets; ++oo) {
            int dx = (int)offs[oo + 0 * nOffsets];
            int dy = (int)offs[oo + 1 * nOffsets];
            int dz = (int)offs[oo + 2 * nOffsets];
            int jx = (int)ix + dx;
            int jy = (int)iy + dy;
            int jz = (int)iz + dz;
            mwSize jBin0;
            mwIndex jHead1;

            if (jx < 0 || jx >= (int)nx ||
                jy < 0 || jy >= (int)ny ||
                jz < 0 || jz >= (int)nz) {
                continue;
            }

            jBin0 = linearize3((mwSize)jx, (mwSize)jy, (mwSize)jz, nx, ny);
            jHead1 = (mwIndex)binHead[jBin0];
            if (jHead1 == 0) continue;

            if (jBin0 == binLin0) {
                mwIndex i1 = iHead1;
                while (i1 != 0) {
                    mwIndex j1 = (mwIndex)binNext[i1 - 1];
                    while (j1 != 0) {
                        mwSize i0 = (mwSize)(i1 - 1);
                        mwSize j0 = (mwSize)(j1 - 1);
                        double dxp = pos[j0 + 0 * nActive] - pos[i0 + 0 * nActive];
                        double dyp = pos[j0 + 1 * nActive] - pos[i0 + 1 * nActive];
                        double dzp = pos[j0 + 2 * nActive] - pos[i0 + 2 * nActive];
                        double r2 = dxp * dxp + dyp * dyp + dzp * dzp;

                        if (r2 <= cutoff2) {
                            double r = sqrt(r2);
                            double invR = 1.0 / r;
                            double invR3 = invR / r2;
                            double invR5 = invR3 / r2;
                            double f3 = 0.0, f5 = 0.0;
                            mwSize k;

                            if (haveThole) {
                                thole_f3f5_scalar(r, alphaAct[i0], alphaAct[j0],
                                                  get_thole_a(tholeA, nA, i0),
                                                  &f3, &f5);
                            }

                            /* i -> j */
                            k = nextPtr[i0]++;
                            colIdxOut[k] = (double)(j0 + 1);
                            drOut[k + 0 * (mwSize)nDirDouble] = dxp;
                            drOut[k + 1 * (mwSize)nDirDouble] = dyp;
                            drOut[k + 2 * (mwSize)nDirDouble] = dzp;
                            r2Out[k] = r2;
                            rOut[k] = r;
                            invR3Out[k] = invR3;
                            invR5Out[k] = invR5;
                            if (haveThole) {
                                f3Out[k] = f3;
                                f5Out[k] = f5;
                            }

                            /* j -> i */
                            k = nextPtr[j0]++;
                            colIdxOut[k] = (double)(i0 + 1);
                            drOut[k + 0 * (mwSize)nDirDouble] = -dxp;
                            drOut[k + 1 * (mwSize)nDirDouble] = -dyp;
                            drOut[k + 2 * (mwSize)nDirDouble] = -dzp;
                            r2Out[k] = r2;
                            rOut[k] = r;
                            invR3Out[k] = invR3;
                            invR5Out[k] = invR5;
                            if (haveThole) {
                                f3Out[k] = f3;
                                f5Out[k] = f5;
                            }
                        }

                        j1 = (mwIndex)binNext[j1 - 1];
                    }
                    i1 = (mwIndex)binNext[i1 - 1];
                }
            } else {
                mwIndex i1 = iHead1;
                while (i1 != 0) {
                    mwIndex j1 = jHead1;
                    while (j1 != 0) {
                        mwSize i0 = (mwSize)(i1 - 1);
                        mwSize j0 = (mwSize)(j1 - 1);
                        double dxp = pos[j0 + 0 * nActive] - pos[i0 + 0 * nActive];
                        double dyp = pos[j0 + 1 * nActive] - pos[i0 + 1 * nActive];
                        double dzp = pos[j0 + 2 * nActive] - pos[i0 + 2 * nActive];
                        double r2 = dxp * dxp + dyp * dyp + dzp * dzp;

                        if (r2 <= cutoff2) {
                            double r = sqrt(r2);
                            double invR = 1.0 / r;
                            double invR3 = invR / r2;
                            double invR5 = invR3 / r2;
                            double f3 = 0.0, f5 = 0.0;
                            mwSize k;

                            if (haveThole) {
                                thole_f3f5_scalar(r, alphaAct[i0], alphaAct[j0],
                                                  get_thole_a(tholeA, nA, i0),
                                                  &f3, &f5);
                            }

                            /* i -> j */
                            k = nextPtr[i0]++;
                            colIdxOut[k] = (double)(j0 + 1);
                            drOut[k + 0 * (mwSize)nDirDouble] = dxp;
                            drOut[k + 1 * (mwSize)nDirDouble] = dyp;
                            drOut[k + 2 * (mwSize)nDirDouble] = dzp;
                            r2Out[k] = r2;
                            rOut[k] = r;
                            invR3Out[k] = invR3;
                            invR5Out[k] = invR5;
                            if (haveThole) {
                                f3Out[k] = f3;
                                f5Out[k] = f5;
                            }

                            /* j -> i */
                            k = nextPtr[j0]++;
                            colIdxOut[k] = (double)(i0 + 1);
                            drOut[k + 0 * (mwSize)nDirDouble] = -dxp;
                            drOut[k + 1 * (mwSize)nDirDouble] = -dyp;
                            drOut[k + 2 * (mwSize)nDirDouble] = -dzp;
                            r2Out[k] = r2;
                            rOut[k] = r;
                            invR3Out[k] = invR3;
                            invR5Out[k] = invR5;
                            if (haveThole) {
                                f3Out[k] = f3;
                                f5Out[k] = f5;
                            }
                        }

                        j1 = (mwIndex)binNext[j1 - 1];
                    }
                    i1 = (mwIndex)binNext[i1 - 1];
                }
            }
        }
    }

    mxFree(rowCounts);
    mxFree(nextPtr);
}