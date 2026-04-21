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
 * Inputs:
 *   0  pos                nSites x 3 double
 *   1  siteMask           nSites x 1 double/logical-like (nonzero = keep)
 *   2  gridShape          1x3 or 3x1 double
 *   3  binHead            nBins x 1 double
 *   4  binNext            nSites x 1 double
 *   5  neighborOffsets    nOffsets x 3 double
 *   6  rcut               scalar double
 *   7  alpha              [] or nSites x 1 double
 *   8  tholeA             [] or scalar or nSites x 1 double
 *
 * Outputs:
 *   0  pair_i
 *   1  pair_j
 *   2  dr
 *   3  r2_bare
 *   4  r_bare
 *   5  inv_r_bare
 *   6  inv_r3_bare
 *   7  inv_r5_bare
 *   8  thole_f3          [] if not requested
 *   9  thole_f5          [] if not requested
 *   10 nPairs
 * ------------------------------------------------------------------------- */

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    const mxArray *posArr, *maskArr, *gridArr, *headArr, *nextArr, *offsArr, *cutArr;
    const mxArray *alphaArr, *aArr;

    const double *pos;
    const double *siteMask;
    const double *gridD;
    const double *binHead;
    const double *binNext;
    const double *offs;
    const double *alpha = NULL;
    const double *tholeA = NULL;

    mwSize nSites, nBins, nOffsets;
    mwSize nx, ny, nz;
    double rcut, cutoff2;
    bool haveThole = false;
    mwSize nA = 0;

    double *pairIOut, *pairJOut, *drOut, *r2Out, *rOut, *invROut, *invR3Out, *invR5Out;
    double *f3Out = NULL, *f5Out = NULL;
    double nPairsDouble = 0.0;

    mwSize binLin0, oo;

    if (nrhs != 9) {
        mexErrMsgIdAndTxt("geom:mex_build_nonperiodic_pair_cache:nrhs",
                          "Expected 9 inputs.");
    }
    if (nlhs != 11) {
        mexErrMsgIdAndTxt("geom:mex_build_nonperiodic_pair_cache:nlhs",
                          "Expected 11 outputs.");
    }

    posArr   = prhs[0];
    maskArr  = prhs[1];
    gridArr  = prhs[2];
    headArr  = prhs[3];
    nextArr  = prhs[4];
    offsArr  = prhs[5];
    cutArr   = prhs[6];
    alphaArr = prhs[7];
    aArr     = prhs[8];

    if (!mxIsDouble(posArr) || mxIsComplex(posArr) || mxGetN(posArr) != 3) {
        mexErrMsgIdAndTxt("geom:mex_build_nonperiodic_pair_cache:pos",
                          "pos must be real double nSites x 3.");
    }
    if (!mxIsDouble(maskArr) || mxGetNumberOfElements(maskArr) != mxGetM(posArr)) {
        mexErrMsgIdAndTxt("geom:mex_build_nonperiodic_pair_cache:mask",
                          "siteMask must be double/logical-like with one entry per site.");
    }
    if (!mxIsDouble(gridArr) || mxGetNumberOfElements(gridArr) != 3) {
        mexErrMsgIdAndTxt("geom:mex_build_nonperiodic_pair_cache:grid",
                          "gridShape must have 3 elements.");
    }
    if (!mxIsDouble(headArr) || !mxIsDouble(nextArr) ||
        !mxIsDouble(offsArr) || !mxIsDouble(cutArr)) {
        mexErrMsgIdAndTxt("geom:mex_build_nonperiodic_pair_cache:type",
                          "Inputs 1:7 must be double.");
    }

    pos = mxGetPr(posArr);
    siteMask = mxGetPr(maskArr);
    nSites = mxGetM(posArr);

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
        mexErrMsgIdAndTxt("geom:mex_build_nonperiodic_pair_cache:rcut",
                          "rcut must be positive and finite.");
    }
    cutoff2 = rcut * rcut;

    if (!mxIsEmpty(alphaArr) && !mxIsEmpty(aArr)) {
        if (!mxIsDouble(alphaArr) || !mxIsDouble(aArr)) {
            mexErrMsgIdAndTxt("geom:mex_build_nonperiodic_pair_cache:tholeType",
                              "alpha and tholeA must be double when provided.");
        }
        if (mxGetNumberOfElements(alphaArr) != nSites) {
            mexErrMsgIdAndTxt("geom:mex_build_nonperiodic_pair_cache:alphaSize",
                              "alpha must have length nSites.");
        }
        nA = mxGetNumberOfElements(aArr);
        if (!(nA == 1 || nA == nSites)) {
            mexErrMsgIdAndTxt("geom:mex_build_nonperiodic_pair_cache:aSize",
                              "tholeA must be scalar or length nSites.");
        }
        haveThole = true;
        alpha = mxGetPr(alphaArr);
        tholeA = mxGetPr(aArr);
    }

    /* ---------------------------------------------------------------------
     * Pass 1: count accepted unordered pairs
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
                    if (siteMask[i1 - 1] != 0.0) {
                        mwIndex j1 = (mwIndex)binNext[i1 - 1];
                        while (j1 != 0) {
                            if (siteMask[j1 - 1] != 0.0) {
                                mwSize i0 = (mwSize)(i1 - 1);
                                mwSize j0 = (mwSize)(j1 - 1);
                                double dxp = pos[j0 + 0 * nSites] - pos[i0 + 0 * nSites];
                                double dyp = pos[j0 + 1 * nSites] - pos[i0 + 1 * nSites];
                                double dzp = pos[j0 + 2 * nSites] - pos[i0 + 2 * nSites];
                                double r2 = dxp * dxp + dyp * dyp + dzp * dzp;

                                if (r2 <= cutoff2) {
                                    nPairsDouble += 1.0;
                                }
                            }
                            j1 = (mwIndex)binNext[j1 - 1];
                        }
                    }
                    i1 = (mwIndex)binNext[i1 - 1];
                }
            } else {
                mwIndex i1 = iHead1;
                while (i1 != 0) {
                    if (siteMask[i1 - 1] != 0.0) {
                        mwIndex j1 = jHead1;
                        while (j1 != 0) {
                            if (siteMask[j1 - 1] != 0.0) {
                                mwSize i0 = (mwSize)(i1 - 1);
                                mwSize j0 = (mwSize)(j1 - 1);
                                double dxp = pos[j0 + 0 * nSites] - pos[i0 + 0 * nSites];
                                double dyp = pos[j0 + 1 * nSites] - pos[i0 + 1 * nSites];
                                double dzp = pos[j0 + 2 * nSites] - pos[i0 + 2 * nSites];
                                double r2 = dxp * dxp + dyp * dyp + dzp * dzp;

                                if (r2 <= cutoff2) {
                                    nPairsDouble += 1.0;
                                }
                            }
                            j1 = (mwIndex)binNext[j1 - 1];
                        }
                    }
                    i1 = (mwIndex)binNext[i1 - 1];
                }
            }
        }
    }

    /* ---------------------------------------------------------------------
     * Allocate outputs
     * --------------------------------------------------------------------- */
    plhs[0]  = mxCreateDoubleMatrix((mwSize)nPairsDouble, 1, mxREAL);
    plhs[1]  = mxCreateDoubleMatrix((mwSize)nPairsDouble, 1, mxREAL);
    plhs[2]  = mxCreateDoubleMatrix((mwSize)nPairsDouble, 3, mxREAL);
    plhs[3]  = mxCreateDoubleMatrix((mwSize)nPairsDouble, 1, mxREAL);
    plhs[4]  = mxCreateDoubleMatrix((mwSize)nPairsDouble, 1, mxREAL);
    plhs[5]  = mxCreateDoubleMatrix((mwSize)nPairsDouble, 1, mxREAL);
    plhs[6]  = mxCreateDoubleMatrix((mwSize)nPairsDouble, 1, mxREAL);
    plhs[7]  = mxCreateDoubleMatrix((mwSize)nPairsDouble, 1, mxREAL);

    pairIOut = mxGetPr(plhs[0]);
    pairJOut = mxGetPr(plhs[1]);
    drOut    = mxGetPr(plhs[2]);
    r2Out    = mxGetPr(plhs[3]);
    rOut     = mxGetPr(plhs[4]);
    invROut  = mxGetPr(plhs[5]);
    invR3Out = mxGetPr(plhs[6]);
    invR5Out = mxGetPr(plhs[7]);

    if (haveThole) {
        plhs[8] = mxCreateDoubleMatrix((mwSize)nPairsDouble, 1, mxREAL);
        plhs[9] = mxCreateDoubleMatrix((mwSize)nPairsDouble, 1, mxREAL);
        f3Out = mxGetPr(plhs[8]);
        f5Out = mxGetPr(plhs[9]);
    } else {
        plhs[8] = mxCreateDoubleMatrix(0, 0, mxREAL);
        plhs[9] = mxCreateDoubleMatrix(0, 0, mxREAL);
    }

    plhs[10] = mxCreateDoubleScalar(nPairsDouble);

    /* ---------------------------------------------------------------------
     * Pass 2: fill outputs
     * --------------------------------------------------------------------- */
    {
        mwSize k = 0;

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
                        if (siteMask[i1 - 1] != 0.0) {
                            mwIndex j1 = (mwIndex)binNext[i1 - 1];
                            while (j1 != 0) {
                                if (siteMask[j1 - 1] != 0.0) {
                                    mwSize i0 = (mwSize)(i1 - 1);
                                    mwSize j0 = (mwSize)(j1 - 1);
                                    double dxp = pos[j0 + 0 * nSites] - pos[i0 + 0 * nSites];
                                    double dyp = pos[j0 + 1 * nSites] - pos[i0 + 1 * nSites];
                                    double dzp = pos[j0 + 2 * nSites] - pos[i0 + 2 * nSites];
                                    double r2 = dxp * dxp + dyp * dyp + dzp * dzp;

                                    if (r2 <= cutoff2) {
                                        double r = sqrt(r2);
                                        double invR = 1.0 / r;
                                        double invR3 = invR / r2;
                                        double invR5 = invR3 / r2;

                                        pairIOut[k] = (double)(i0 + 1);
                                        pairJOut[k] = (double)(j0 + 1);

                                        drOut[k + 0 * (mwSize)nPairsDouble] = dxp;
                                        drOut[k + 1 * (mwSize)nPairsDouble] = dyp;
                                        drOut[k + 2 * (mwSize)nPairsDouble] = dzp;

                                        r2Out[k] = r2;
                                        rOut[k] = r;
                                        invROut[k] = invR;
                                        invR3Out[k] = invR3;
                                        invR5Out[k] = invR5;

                                        if (haveThole) {
                                            double f3, f5;
                                            thole_f3f5_scalar(r, alpha[i0], alpha[j0],
                                                              get_thole_a(tholeA, nA, i0),
                                                              &f3, &f5);
                                            f3Out[k] = f3;
                                            f5Out[k] = f5;
                                        }

                                        k++;
                                    }
                                }
                                j1 = (mwIndex)binNext[j1 - 1];
                            }
                        }
                        i1 = (mwIndex)binNext[i1 - 1];
                    }
                } else {
                    mwIndex i1 = iHead1;
                    while (i1 != 0) {
                        if (siteMask[i1 - 1] != 0.0) {
                            mwIndex j1 = jHead1;
                            while (j1 != 0) {
                                if (siteMask[j1 - 1] != 0.0) {
                                    mwSize i0 = (mwSize)(i1 - 1);
                                    mwSize j0 = (mwSize)(j1 - 1);
                                    double dxp = pos[j0 + 0 * nSites] - pos[i0 + 0 * nSites];
                                    double dyp = pos[j0 + 1 * nSites] - pos[i0 + 1 * nSites];
                                    double dzp = pos[j0 + 2 * nSites] - pos[i0 + 2 * nSites];
                                    double r2 = dxp * dxp + dyp * dyp + dzp * dzp;

                                    if (r2 <= cutoff2) {
                                        double r = sqrt(r2);
                                        double invR = 1.0 / r;
                                        double invR3 = invR / r2;
                                        double invR5 = invR3 / r2;

                                        pairIOut[k] = (double)(i0 + 1);
                                        pairJOut[k] = (double)(j0 + 1);

                                        drOut[k + 0 * (mwSize)nPairsDouble] = dxp;
                                        drOut[k + 1 * (mwSize)nPairsDouble] = dyp;
                                        drOut[k + 2 * (mwSize)nPairsDouble] = dzp;

                                        r2Out[k] = r2;
                                        rOut[k] = r;
                                        invROut[k] = invR;
                                        invR3Out[k] = invR3;
                                        invR5Out[k] = invR5;

                                        if (haveThole) {
                                            double f3, f5;
                                            thole_f3f5_scalar(r, alpha[i0], alpha[j0],
                                                              get_thole_a(tholeA, nA, i0),
                                                              &f3, &f5);
                                            f3Out[k] = f3;
                                            f5Out[k] = f5;
                                        }

                                        k++;
                                    }
                                }
                                j1 = (mwIndex)binNext[j1 - 1];
                            }
                        }
                        i1 = (mwIndex)binNext[i1 - 1];
                    }
                }
            }
        }

        if ((double)k != nPairsDouble) {
            mexErrMsgIdAndTxt("geom:mex_build_nonperiodic_pair_cache:countMismatch",
                              "Internal pair-count mismatch.");
        }
    }
}