#include "mex.h"
#include <math.h>
#include <stdbool.h>

static int mod_floor_int(int a, int n) {
    int r = a % n;
    return (r < 0) ? (r + n) : r;
}

static int floor_div_int(int a, int n) {
    int q = a / n;
    int r = a % n;
    if (r != 0 && ((r > 0) != (n > 0))) q -= 1;
    return q;
}

static void matvec3x3(const double *H, int nx, int ny, int nz, double *out) {
    out[0] = H[0]*nx + H[3]*ny + H[6]*nz;
    out[1] = H[1]*nx + H[4]*ny + H[7]*nz;
    out[2] = H[2]*nx + H[5]*ny + H[8]*nz;
}

static void thole_l3_l5(double r, double ai, double aj, double thole_a, double *l3, double *l5) {
    *l3 = 0.0;
    *l5 = 0.0;
    if (thole_a == 0.0 || ai <= 0.0 || aj <= 0.0 || r <= 0.0) return;

    double u = r / pow(ai * aj, 1.0/6.0);
    double au3 = thole_a * u*u*u;
    double e = exp(-au3);
    double f3 = 1.0 - e;
    double f5 = 1.0 - (1.0 + au3) * e;
    *l3 = f3 - 1.0;
    *l5 = f5 - 1.0;
}

/*
Inputs:
 0 posAct            nActive x 3 Cartesian active positions
 1 fracActWrapped    nActive x 3 wrapped fractional coordinates in [0,1)
 2 grid_shape        1x3
 3 bin_head          nBins x 1
 4 bin_next          nActive x 1
 5 neighbor_offsets  nOff x 3 integer-ish doubles
 6 H                 3x3 lattice
 7 rcut              scalar
 8 alphaEwald        scalar
 9 alphaAct          nActive x 1 (or empty)
10 tholeAAct         scalar or nActive x 1 (or empty)

Outputs:
 0 row_ptr
 1 col_idx
 2 source_idx_active
 3 dr
 4 r2_bare
 5 r_bare
 6 coeff_iso
 7 coeff_dyad
 8 nEntriesDirected
 9 nCandidatesVisited
10 nCandidatesWithinCutoff
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    const double *pos, *frac, *gridShapeD, *binHeadD, *binNextD, *offD, *H;
    const double *alphaAct = NULL, *tholeAAct = NULL;
    mwSize nActive, nBins, nOff;
    int gx, gy, gz;
    double rcut, rcut2, alphaEwald, alpha2, twoAlphaOverSqrtPi;
    bool haveAlpha = false, haveThole = false, tholeScalar = false;

    mwSize *rowCounts, *nextPtr;
    double *row_ptr, *col_idx, *source_idx, *dr, *r2_bare, *r_bare, *coeff_iso, *coeff_dyad;
    mwSize i, s;
    mwSize nEntriesDirected = 0;
    mwSize nCandidatesVisited = 0;
    mwSize nCandidatesWithinCutoff = 0;

    if (nrhs != 11) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:nrhs",
                          "Expected 11 inputs.");
    }
    if (nlhs != 11) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:nlhs",
                          "Expected 11 outputs.");
    }

    pos = mxGetPr(prhs[0]);
    nActive = mxGetM(prhs[0]);
    if (mxGetN(prhs[0]) != 3) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:pos",
                          "posAct must be nActive x 3.");
    }

    frac = mxGetPr(prhs[1]);
    if (mxGetM(prhs[1]) != nActive || mxGetN(prhs[1]) != 3) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:frac",
                          "fracActWrapped must be nActive x 3.");
    }

    gridShapeD = mxGetPr(prhs[2]);
    if (mxGetNumberOfElements(prhs[2]) != 3) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:grid",
                          "grid_shape must have 3 elements.");
    }
    gx = (int)gridShapeD[0];
    gy = (int)gridShapeD[1];
    gz = (int)gridShapeD[2];
    nBins = (mwSize)(gx * gy * gz);

    binHeadD = mxGetPr(prhs[3]);
    if (mxGetNumberOfElements(prhs[3]) != nBins) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:bin_head",
                          "bin_head length mismatch.");
    }

    binNextD = mxGetPr(prhs[4]);
    if (mxGetNumberOfElements(prhs[4]) != nActive) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:bin_next",
                          "bin_next length mismatch.");
    }

    offD = mxGetPr(prhs[5]);
    nOff = mxGetM(prhs[5]);
    if (mxGetN(prhs[5]) != 3) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:neighbor_offsets",
                          "neighbor_offsets must be nOff x 3.");
    }

    H = mxGetPr(prhs[6]);
    if (mxGetM(prhs[6]) != 3 || mxGetN(prhs[6]) != 3) {
        mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:H",
                          "H must be 3x3.");
    }

    rcut = mxGetScalar(prhs[7]);
    alphaEwald = mxGetScalar(prhs[8]);
    rcut2 = rcut * rcut;
    alpha2 = alphaEwald * alphaEwald;
    twoAlphaOverSqrtPi = 2.0 * alphaEwald / sqrt(M_PI);

    if (!mxIsEmpty(prhs[9])) {
        alphaAct = mxGetPr(prhs[9]);
        if (mxGetNumberOfElements(prhs[9]) != nActive) {
            mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:alphaAct",
                              "alphaAct must have length nActive.");
        }
        haveAlpha = true;
    }

    if (!mxIsEmpty(prhs[10])) {
        tholeAAct = mxGetPr(prhs[10]);
        if (mxGetNumberOfElements(prhs[10]) == 1) {
            tholeScalar = true;
        } else if (mxGetNumberOfElements(prhs[10]) == nActive) {
            tholeScalar = false;
        } else {
            mexErrMsgIdAndTxt("geom:mex_build_active_row_cache_periodic:tholeAAct",
                              "tholeAAct must be scalar or length nActive.");
        }
        haveThole = true;
    }

    rowCounts = (mwSize*) mxCalloc((nActive > 0 ? nActive : 1), sizeof(mwSize));

    /* Pass 1: count row entries */
    for (i = 0; i < nActive; ++i) {
        int bix = (int)floor(frac[i] * gx);
        int biy = (int)floor(frac[i + nActive] * gy);
        int biz = (int)floor(frac[i + 2*nActive] * gz);

        if (bix >= gx) bix = gx - 1;
        if (biy >= gy) biy = gy - 1;
        if (biz >= gz) biz = gz - 1;

        for (s = 0; s < nOff; ++s) {
            int dx = (int)offD[s];
            int dy = (int)offD[s + nOff];
            int dz = (int)offD[s + 2*nOff];

            int rawx = bix + dx;
            int rawy = biy + dy;
            int rawz = biz + dz;

            int imgx = floor_div_int(rawx, gx);
            int imgy = floor_div_int(rawy, gy);
            int imgz = floor_div_int(rawz, gz);

            int nbx = mod_floor_int(rawx, gx);
            int nby = mod_floor_int(rawy, gy);
            int nbz = mod_floor_int(rawz, gz);

            mwSize binLin = (mwSize)(nbx + 1) +
                            (mwSize)nby * (mwSize)gx +
                            (mwSize)nbz * (mwSize)gx * (mwSize)gy;

            mwSize j = (mwSize)binHeadD[binLin - 1];
            while (j != 0) {
                mwSize jj = j - 1;
                double shift[3], ddx, ddy, ddz, rr2;
                matvec3x3(H, imgx, imgy, imgz, shift);

                ddx = pos[jj] - pos[i] + shift[0];
                ddy = pos[jj + nActive] - pos[i + nActive] + shift[1];
                ddz = pos[jj + 2*nActive] - pos[i + 2*nActive] + shift[2];
                rr2 = ddx*ddx + ddy*ddy + ddz*ddz;

                if (!(jj == i && imgx == 0 && imgy == 0 && imgz == 0)) {
                    nCandidatesVisited += 1;
                    if (rr2 <= rcut2) {
                        rowCounts[i] += 1;
                        nEntriesDirected += 1;
                        nCandidatesWithinCutoff += 1;
                    }
                }

                j = (mwSize)binNextD[jj];
            }
        }
    }

    plhs[0] = mxCreateDoubleMatrix(nActive + 1, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nEntriesDirected, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(nEntriesDirected, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(nEntriesDirected, 3, mxREAL);
    plhs[4] = mxCreateDoubleMatrix(nEntriesDirected, 1, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(nEntriesDirected, 1, mxREAL);
    plhs[6] = mxCreateDoubleMatrix(nEntriesDirected, 1, mxREAL);
    plhs[7] = mxCreateDoubleMatrix(nEntriesDirected, 1, mxREAL);
    plhs[8] = mxCreateDoubleScalar((double)nEntriesDirected);
    plhs[9]  = mxCreateDoubleScalar((double)nCandidatesVisited);
    plhs[10] = mxCreateDoubleScalar((double)nCandidatesWithinCutoff);

    row_ptr = mxGetPr(plhs[0]);
    col_idx = mxGetPr(plhs[1]);
    source_idx = mxGetPr(plhs[2]);
    dr = mxGetPr(plhs[3]);
    r2_bare = mxGetPr(plhs[4]);
    r_bare = mxGetPr(plhs[5]);
    coeff_iso = mxGetPr(plhs[6]);
    coeff_dyad = mxGetPr(plhs[7]);

    row_ptr[0] = 1.0;
    for (i = 0; i < nActive; ++i) {
        row_ptr[i + 1] = row_ptr[i] + (double)rowCounts[i];
    }

    nextPtr = (mwSize*) mxCalloc((nActive > 0 ? nActive : 1), sizeof(mwSize));
    for (i = 0; i < nActive; ++i) {
        nextPtr[i] = (mwSize)(row_ptr[i] - 1.0);
    }

    /* Pass 2: fill */
    for (i = 0; i < nActive; ++i) {
        int bix = (int)floor(frac[i] * gx);
        int biy = (int)floor(frac[i + nActive] * gy);
        int biz = (int)floor(frac[i + 2*nActive] * gz);

        double ai = haveAlpha ? alphaAct[i] : 0.0;
        double thi = haveThole ? (tholeScalar ? tholeAAct[0] : tholeAAct[i]) : 0.0;

        if (bix >= gx) bix = gx - 1;
        if (biy >= gy) biy = gy - 1;
        if (biz >= gz) biz = gz - 1;

        for (s = 0; s < nOff; ++s) {
            int dx = (int)offD[s];
            int dy = (int)offD[s + nOff];
            int dz = (int)offD[s + 2*nOff];

            int rawx = bix + dx;
            int rawy = biy + dy;
            int rawz = biz + dz;

            int imgx = floor_div_int(rawx, gx);
            int imgy = floor_div_int(rawy, gy);
            int imgz = floor_div_int(rawz, gz);

            int nbx = mod_floor_int(rawx, gx);
            int nby = mod_floor_int(rawy, gy);
            int nbz = mod_floor_int(rawz, gz);

            mwSize binLin = (mwSize)(nbx + 1) +
                            (mwSize)nby * (mwSize)gx +
                            (mwSize)nbz * (mwSize)gx * (mwSize)gy;

            mwSize j = (mwSize)binHeadD[binLin - 1];
            while (j != 0) {
                mwSize jj1 = j;
                mwSize jj = j - 1;

                double shift[3], ddx, ddy, ddz, rr2;

                matvec3x3(H, imgx, imgy, imgz, shift);

                ddx = pos[jj] - pos[i] + shift[0];
                ddy = pos[jj + nActive] - pos[i + nActive] + shift[1];
                ddz = pos[jj + 2*nActive] - pos[i + 2*nActive] + shift[2];
                rr2 = ddx*ddx + ddy*ddy + ddz*ddz;

                if (!(jj == i && imgx == 0 && imgy == 0 && imgz == 0) && rr2 <= rcut2) {
                    double rr = sqrt(rr2);
                    double invr = 1.0 / rr;
                    double invr2 = 1.0 / rr2;
                    double invr3 = invr * invr2;
                    double invr5 = invr3 * invr2;
                    double invr4 = invr2 * invr2;
                    double erfcar = erfc(alphaEwald * rr);
                    double expar2 = exp(-alpha2 * rr2);
                    double B = erfcar * invr3 + twoAlphaOverSqrtPi * expar2 * invr2;
                    double C = 3.0 * erfcar * invr5 +
                               twoAlphaOverSqrtPi * (2.0 * alpha2 * invr2 + 3.0 * invr4) * expar2;

                    double ci = -B;
                    double cd = +C;

                    if (haveAlpha && haveThole) {
                        double aj = alphaAct[jj];
                        double thj = tholeScalar ? tholeAAct[0] : tholeAAct[jj];
                        double l3, l5;
                        double th = thi;
                        if (thj != thi) th = thi;

                        thole_l3_l5(rr, ai, aj, th, &l3, &l5);
                        ci -= l3 * invr3;
                        cd += 3.0 * l5 * invr5;
                    }

                    mwSize k = nextPtr[i]++;
                    col_idx[k] = (double)jj1;
                    source_idx[k] = (double)jj1;

                    dr[k] = ddx;
                    dr[k + nEntriesDirected] = ddy;
                    dr[k + 2*nEntriesDirected] = ddz;

                    r2_bare[k] = rr2;
                    r_bare[k] = rr;
                    coeff_iso[k] = ci;
                    coeff_dyad[k] = cd;
                }

                j = (mwSize)binNextD[jj];
            }
        }
    }

    mxFree(rowCounts);
    mxFree(nextPtr);
}