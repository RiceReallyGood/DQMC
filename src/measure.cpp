#include "measure.h"

#include <cstring>

#include "mkl.h"
#include "util.h"

TDM::TDM(int n_, int M_, int npack_)
    : n(n_), M(M_), npack(npack_), nblk(M_ / npack_), count(0), sign(0) {
    inc = npack_ / 2;
    dec = npack_ - inc - 1;
    gtup = new double[M * n * n];
    gtdn = new double[M * n * n];

    SzSz = new double[M * n * n];
    denden = new double[M * n * n];

    nunu = new double[M * n * n];
    nund = new double[M * n * n];
    ndnd = new double[M * n * n];
}

TDM::~TDM() {
    delete[] denden;
    delete[] SzSz;
    delete[] ndnd;
    delete[] nund;
    delete[] nunu;
    delete[] gtup;
    delete[] gtdn;
}

void TDM::measure(Gtau& gt_up, Gtau& gt_dn) {
    // std::cout << "inc = " << inc << ", dec = " << dec << std::endl;
    for (int blk0 = 0; blk0 < nblk; blk0++) {
        for (int blkt = 0; blkt < nblk; blkt++) {
            gt_up.Load_G_to_g(blk0, blkt);
            gt_dn.Load_G_to_g(blk0, blkt);
            record(gt_up, gt_dn);

            for (int i = 0; i < inc; i++) {
                gt_up.Increase_g_time();
                gt_dn.Increase_g_time();
                record(gt_up, gt_dn);
            }

            if (dec > 0) {
                gt_up.Load_G_to_g(blk0, blkt);
                gt_dn.Load_G_to_g(blk0, blkt);
            }

            for (int i = 0; i < dec; i++) {
                gt_up.Decrease_g_time();
                gt_dn.Decrease_g_time();
                record(gt_up, gt_dn);
            }
        }
    }
    count += nblk;
}

void TDM::record(Gtau& gt_up, Gtau& gt_dn) {
    double curr_sign = gt_up.sign * gt_dn.sign;
    sign += curr_sign;
    /* only rows */
    // int l;
    // // gtau
    // if (gt_up.lt >= gt_up.l0) {
    //     l = gt_up.lt - gt_up.l0;
    //     cblas_daxpy(n * n, 1., gt_up.gt0, 1, &gtup[l * n * n], 1);
    //     cblas_daxpy(n * n, 1., gt_dn.gt0, 1, &gtdn[l * n * n], 1);
    // }
    // else {
    //     l = M + gt_up.lt - gt_up.l0;
    //     cblas_daxpy(n * n, -1., gt_up.gt0, 1, &gtup[l * n * n], 1);
    //     cblas_daxpy(n * n, -1., gt_dn.gt0, 1, &gtdn[l * n * n], 1);
    // }

    /* rows and cols */
    int l;
    if (gt_up.lt == gt_up.l0) {
        l = 0;
        cblas_daxpy(n * n, 1. * curr_sign, gt_up.gt0, 1, gtup, 1);
        cblas_daxpy(n * n, 1. * curr_sign, gt_dn.gt0, 1, gtdn, 1);
    }
    else {
        double factor;
        if (gt_up.lt > gt_up.l0) {
            l = gt_up.lt - gt_up.l0;
            factor = 0.5 * curr_sign;
        }
        else {
            l = M + gt_up.lt - gt_up.l0;
            factor = -0.5 * curr_sign;
        }
        int invl = (M - l) % M;

        cblas_daxpy(n * n, factor, gt_up.gt0, 1, &gtup[l * n * n], 1);
        cblas_daxpy(n * n, -factor, gt_up.g0t, 1, &gtup[invl * n * n], 1);

        cblas_daxpy(n * n, factor, gt_dn.gt0, 1, &gtdn[l * n * n], 1);
        cblas_daxpy(n * n, -factor, gt_dn.g0t, 1, &gtdn[invl * n * n], 1);
    }
    // std::cout << "gt_up.lt = " << gt_up.lt << ", gt_up.l0 = " << gt_up.l0
    //           << ", l = " << l << std::endl;
    // printmat("meas gt0", n, n, n, gt_up.gt0);

    // nunu
    double* destuu = &nunu[l * n * n];
    double* destud = &nund[l * n * n];
    double* destdd = &ndnd[l * n * n];
    for (int x = 0; x < n; x++) {
        for (int y = 0; y < n; y++) {
            destuu[x * n + y] +=
                curr_sign * (gt_up.gtt[x * n + x] * gt_up.g00[y * n + y] -
                             gt_up.gt0[x * n + y] * gt_up.g0t[y * n + x]);
            destdd[x * n + y] +=
                curr_sign * (gt_dn.gtt[x * n + x] * gt_dn.g00[y * n + y] -
                             gt_dn.gt0[x * n + y] * gt_dn.g0t[y * n + x]);
            destud[x * n + y] +=
                curr_sign * (gt_up.gtt[x * n + x] * gt_dn.g00[y * n + y]);
        }
    }
}

void TDM::finalize() {
    sign /= double(count * M);
    double* WS = new double[M * n * n];

    // Symmetrization of G up
    for (int l = 0; l < M; l++) {
        Trans(n, &WS[l * n * n], &gtup[l * n * n]);
    }
    vdAdd(M * n * n, gtup, WS, gtup);
    cblas_dscal(M * n * n, 0.5 / (double(count) * sign), gtup, 1);

    // Symmetrization of G down
    for (int l = 0; l < M; l++) {
        Trans(n, &WS[l * n * n], &gtdn[l * n * n]);
    }
    vdAdd(M * n * n, gtdn, WS, gtdn);
    cblas_dscal(M * n * n, 0.5 / (double(count) * sign), gtdn, 1);

    cblas_dscal(M * n * n, 1. / (double(count) * sign), nunu, 1);
    cblas_dscal(M * n * n, 1. / (double(count) * sign), ndnd, 1);
    cblas_dscal(M * n * n, 1. / (double(count) * sign), nund, 1);

    for (int l = 0; l < M; l++) {
        Trans(n, &WS[l * n * n], &nund[(M - l) % M * n * n]);
    }

    memcpy(SzSz, nunu, M * n * n * sizeof(double));
    vdAdd(M * n * n, SzSz, ndnd, SzSz);
    vdSub(M * n * n, SzSz, nund, SzSz);
    vdSub(M * n * n, SzSz, WS, SzSz);

    memcpy(denden, nunu, M * n * n * sizeof(double));
    vdAdd(M * n * n, denden, ndnd, denden);
    vdAdd(M * n * n, denden, nund, denden);
    vdAdd(M * n * n, denden, WS, denden);

    // Symmetrization of SzSz
    for (int l = 0; l < M; l++) {
        Trans(n, &WS[l * n * n], &SzSz[l * n * n]);
    }
    vdAdd(M * n * n, SzSz, WS, SzSz);
    cblas_dscal(M * n * n, 0.5, SzSz, 1);

    // Symmetrization of denden
    for (int l = 0; l < M; l++) {
        Trans(n, &WS[l * n * n], &denden[l * n * n]);
    }
    vdAdd(M * n * n, denden, WS, denden);
    cblas_dscal(M * n * n, 0.5, denden, 1);

    delete[] WS;
}

void TDM::clear() {
    memset(gtup, 0, M * n * n * sizeof(double));
    memset(gtdn, 0, M * n * n * sizeof(double));
    memset(nunu, 0, M * n * n * sizeof(double));
    memset(nund, 0, M * n * n * sizeof(double));
    memset(ndnd, 0, M * n * n * sizeof(double));
    count = 0;
    sign = 0.;
}