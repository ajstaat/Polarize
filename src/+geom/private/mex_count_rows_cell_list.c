#include "mex.h"
#include <math.h>

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

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    const mxArray *posArr;
    const mxArray *gridArr;
    const mxArray *headArr;
    const mxArray *nextArr;
    const mxArray *offsArr;
    const mxArray *cutArr;

    const double *pos;
    const double *gridD;
    const double *binHead;
    const double *binNext;
    const double *offs;
    double cutoff, cutoff2;

    mwSize nActive, nx, ny, nz, nOffsets, nBins;
    double *rowCounts;
    double nPairsDouble = 0.0;

    mwSize binLin0, oo;

    if (nrhs != 6) {
        mexErrMsgIdAndTxt("geom:mex_count_rows_cell_list:nrhs",
                          "Expected 6 inputs: pos, gridShape, binHead, binNext, neighborOffsets, cutoff.");
    }
    if (nlhs > 2) {
        mexErrMsgIdAndTxt("geom:mex_count_rows_cell_list:nlhs",
                          "Returns [rowCounts, nPairs].");
    }

    posArr  = prhs[0];
    gridArr = prhs[1];
    headArr = prhs[2];
    nextArr = prhs[3];
    offsArr = prhs[4];
    cutArr  = prhs[5];

    if (!mxIsDouble(posArr) || mxIsComplex(posArr) || mxGetN(posArr) != 3) {
        mexErrMsgIdAndTxt("geom:mex_count_rows_cell_list:pos",
                          "pos must be real double N x 3.");
    }
    if (!mxIsDouble(gridArr) || mxGetNumberOfElements(gridArr) != 3) {
        mexErrMsgIdAndTxt("geom:mex_count_rows_cell_list:gridShape",
                          "gridShape must have 3 elements.");
    }
    if (!mxIsDouble(headArr) || !mxIsDouble(nextArr) ||
        !mxIsDouble(offsArr) || !mxIsDouble(cutArr)) {
        mexErrMsgIdAndTxt("geom:mex_count_rows_cell_list:type",
                          "All inputs must be double.");
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
    nOffsets = mxGetM(offsArr);

    cutoff = mxGetScalar(cutArr);
    cutoff2 = cutoff * cutoff;

    plhs[0] = mxCreateDoubleMatrix(nActive, 1, mxREAL);
    rowCounts = mxGetPr(plhs[0]);

    nBins = mxGetNumberOfElements(headArr);

    for (binLin0 = 0; binLin0 < nBins; ++binLin0) {
        mwIndex iHead1 = (mwIndex)binHead[binLin0];
        mwSize ix, iy, iz;

        if (iHead1 == 0) {
            continue;
        }

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
            if (jHead1 == 0) {
                continue;
            }

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
                            nPairsDouble += 1.0;
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
                            nPairsDouble += 1.0;
                        }

                        j1 = (mwIndex)binNext[j1 - 1];
                    }
                    i1 = (mwIndex)binNext[i1 - 1];
                }
            }
        }
    }

    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleScalar(nPairsDouble);
    }
}