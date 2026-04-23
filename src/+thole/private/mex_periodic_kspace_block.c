#include "mex.h"
#include <string.h>

/*
 * mex_periodic_kspace_block.c
 *
 * Two modes:
 *
 * 1) apply_row
 *    Erow = mex_periodic_kspace_block('apply_row', ...
 *             Ablk, Bblk, cos_i, sin_i, two_pref_blk, kvecs_blk)
 *
 *    Inputs:
 *      Ablk         1 x nk or nk x 1 double
 *      Bblk         1 x nk or nk x 1 double
 *      cos_i        1 x nk or nk x 1 double
 *      sin_i        1 x nk or nk x 1 double
 *      two_pref_blk 1 x nk or nk x 1 double
 *      kvecs_blk    nk x 3 double
 *
 *    Output:
 *      Erow         1 x 3 double
 *
 * 2) update_ab
 *    [Ablk_out, Bblk_out] = mex_periodic_kspace_block('update_ab', ...
 *             Ablk, Bblk, delta_i, cos_i, sin_i, kvecsT_blk)
 *
 *    Inputs:
 *      Ablk         1 x nk or nk x 1 double
 *      Bblk         1 x nk or nk x 1 double
 *      delta_i      1 x 3 or 3 x 1 double
 *      cos_i        1 x nk or nk x 1 double
 *      sin_i        1 x nk or nk x 1 double
 *      kvecsT_blk   3 x nk double
 *
 *    Outputs:
 *      Ablk_out     same shape as Ablk input
 *      Bblk_out     same shape as Bblk input
 */

static void require_double_real(const mxArray *arr, const char *name)
{
    if (!mxIsDouble(arr) || mxIsComplex(arr)) {
        mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:TypeError",
                          "%s must be a real double array.", name);
    }
}

static mwSize vector_length(const mxArray *arr, const char *name)
{
    mwSize m = mxGetM(arr);
    mwSize n = mxGetN(arr);

    if (!(m == 1 || n == 1)) {
        mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:ShapeError",
                          "%s must be a vector.", name);
    }
    return (m > n) ? m : n;
}

static void copy_shape_like(const mxArray *src, mxArray **dst, double **dstp)
{
    mwSize m = mxGetM(src);
    mwSize n = mxGetN(src);
    *dst = mxCreateDoubleMatrix(m, n, mxREAL);
    *dstp = mxGetPr(*dst);
}

static void do_apply_row(int nlhs, mxArray *plhs[],
                         int nrhs, const mxArray *prhs[])
{
    /* prhs:
       0 mode string
       1 Ablk
       2 Bblk
       3 cos_i
       4 sin_i
       5 two_pref_blk
       6 kvecs_blk
    */
    const mxArray *Aarr = prhs[1];
    const mxArray *Barr = prhs[2];
    const mxArray *Carr = prhs[3];
    const mxArray *Sarr = prhs[4];
    const mxArray *Parr = prhs[5];
    const mxArray *Karr = prhs[6];

    require_double_real(Aarr, "Ablk");
    require_double_real(Barr, "Bblk");
    require_double_real(Carr, "cos_i");
    require_double_real(Sarr, "sin_i");
    require_double_real(Parr, "two_pref_blk");
    require_double_real(Karr, "kvecs_blk");

    mwSize nkA = vector_length(Aarr, "Ablk");
    mwSize nkB = vector_length(Barr, "Bblk");
    mwSize nkC = vector_length(Carr, "cos_i");
    mwSize nkS = vector_length(Sarr, "sin_i");
    mwSize nkP = vector_length(Parr, "two_pref_blk");

    mwSize kRows = mxGetM(Karr);
    mwSize kCols = mxGetN(Karr);

    if (kCols != 3) {
        mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:ShapeError",
                          "kvecs_blk must be nk x 3.");
    }

    if (!(nkA == nkB && nkA == nkC && nkA == nkS && nkA == nkP && nkA == kRows)) {
        mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:SizeMismatch",
                          "All block vectors must have the same length nk, and kvecs_blk must be nk x 3.");
    }

    const double *A = mxGetPr(Aarr);
    const double *B = mxGetPr(Barr);
    const double *C = mxGetPr(Carr);
    const double *S = mxGetPr(Sarr);
    const double *P = mxGetPr(Parr);
    const double *K = mxGetPr(Karr);

    plhs[0] = mxCreateDoubleMatrix(1, 3, mxREAL);
    double *E = mxGetPr(plhs[0]);

    /* MATLAB column-major for nk x 3:
       col 0: K[0 ... nk-1]         = kx
       col 1: K[nk ... 2nk-1]       = ky
       col 2: K[2nk ... 3nk-1]      = kz
    */
    double Ex = 0.0, Ey = 0.0, Ez = 0.0;

    for (mwSize j = 0; j < nkA; ++j) {
        double phase_factor = C[j] * A[j] + S[j] * B[j];
        double w = phase_factor * P[j];

        Ex += w * K[j];
        Ey += w * K[j + nkA];
        Ez += w * K[j + 2 * nkA];
    }

    E[0] = Ex;
    E[1] = Ey;
    E[2] = Ez;

    if (nlhs > 1) {
        mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:TooManyOutputs",
                          "apply_row returns exactly one output.");
    }
}

static void do_update_ab(int nlhs, mxArray *plhs[],
                         int nrhs, const mxArray *prhs[])
{
    /* prhs:
       0 mode string
       1 Ablk
       2 Bblk
       3 delta_i
       4 cos_i
       5 sin_i
       6 kvecsT_blk
    */
    const mxArray *Aarr = prhs[1];
    const mxArray *Barr = prhs[2];
    const mxArray *Darr = prhs[3];
    const mxArray *Carr = prhs[4];
    const mxArray *Sarr = prhs[5];
    const mxArray *KTarr = prhs[6];

    require_double_real(Aarr, "Ablk");
    require_double_real(Barr, "Bblk");
    require_double_real(Darr, "delta_i");
    require_double_real(Carr, "cos_i");
    require_double_real(Sarr, "sin_i");
    require_double_real(KTarr, "kvecsT_blk");

    mwSize nkA = vector_length(Aarr, "Ablk");
    mwSize nkB = vector_length(Barr, "Bblk");
    mwSize nkC = vector_length(Carr, "cos_i");
    mwSize nkS = vector_length(Sarr, "sin_i");
    mwSize nD  = vector_length(Darr, "delta_i");

    mwSize ktRows = mxGetM(KTarr);
    mwSize ktCols = mxGetN(KTarr);

    if (!(nD == 3)) {
        mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:ShapeError",
                          "delta_i must have length 3.");
    }
    if (!(ktRows == 3)) {
        mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:ShapeError",
                          "kvecsT_blk must be 3 x nk.");
    }
    if (!(nkA == nkB && nkA == nkC && nkA == nkS && nkA == ktCols)) {
        mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:SizeMismatch",
                          "Ablk, Bblk, cos_i, sin_i must have length nk, and kvecsT_blk must be 3 x nk.");
    }

    const double *A = mxGetPr(Aarr);
    const double *B = mxGetPr(Barr);
    const double *D = mxGetPr(Darr);
    const double *C = mxGetPr(Carr);
    const double *S = mxGetPr(Sarr);
    const double *KT = mxGetPr(KTarr);

    double dx = D[0];
    double dy = D[1];
    double dz = D[2];

    double *Aout;
    double *Bout;
    copy_shape_like(Aarr, &plhs[0], &Aout);
    copy_shape_like(Barr, &plhs[1], &Bout);

    /* MATLAB column-major for 3 x nk:
       col j = KT[0 + 3*j], KT[1 + 3*j], KT[2 + 3*j]
    */
    for (mwSize j = 0; j < nkA; ++j) {
        double delta_v = dx * KT[0 + 3 * j] +
                         dy * KT[1 + 3 * j] +
                         dz * KT[2 + 3 * j];

        Aout[j] = A[j] + C[j] * delta_v;
        Bout[j] = B[j] + S[j] * delta_v;
    }

    if (nlhs != 2) {
        mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:BadNumOutputs",
                          "update_ab returns exactly two outputs.");
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs < 1) {
        mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:BadNumInputs",
                          "First input must be a mode string.");
    }
    if (!mxIsChar(prhs[0])) {
        mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:BadMode",
                          "First input must be a mode string.");
    }

    char mode[64];
    if (mxGetString(prhs[0], mode, sizeof(mode)) != 0) {
        mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:BadMode",
                          "Could not read mode string.");
    }

    if (strcmp(mode, "apply_row") == 0) {
        if (nrhs != 7) {
            mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:BadNumInputs",
                              "apply_row expects 7 inputs total.");
        }
        if (nlhs != 1) {
            mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:BadNumOutputs",
                              "apply_row returns exactly one output.");
        }
        do_apply_row(nlhs, plhs, nrhs, prhs);
        return;
    }

    if (strcmp(mode, "update_ab") == 0) {
        if (nrhs != 7) {
            mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:BadNumInputs",
                              "update_ab expects 7 inputs total.");
        }
        if (nlhs != 2) {
            mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:BadNumOutputs",
                              "update_ab returns exactly two outputs.");
        }
        do_update_ab(nlhs, plhs, nrhs, prhs);
        return;
    }

    mexErrMsgIdAndTxt("thole:mex_periodic_kspace_block:UnknownMode",
                      "Unknown mode. Supported modes are 'apply_row' and 'update_ab'.");
}