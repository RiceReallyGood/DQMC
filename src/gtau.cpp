#include "gtau.h"

#include <cstring>

#include "mkl.h"
#include "util.h"

Gtau::Gtau(const Model* model, double*& emV_, int npack_)
    : n(model->nsites),
      M(model->M),
      emV(emV_),
      npack(npack_),
      nblk(model->M / npack_),
      dim(model->nsites * model->M / npack_),
      mb(model),
      l0(0),
      lt(0),
      sign(1) {
    G = new double[dim * dim];
    g00 = new double[n * n];
    gt0 = new double[n * n];
    g0t = new double[n * n];
    gtt = new double[n * n];

    WS = new double[n * n];
    pvt = new int[dim];
}

Gtau::~Gtau() {
    delete[] pvt;
    delete[] WS;

    delete[] gtt;
    delete[] g0t;
    delete[] gt0;
    delete[] g00;
    delete[] G;
}

void Gtau::Get_FullG() {
    Eye(dim, G);
    int l = 0;
    for (int jblk = 0; jblk < nblk; jblk++) {
        mb.GetB(&emV[l * n], WS);
        l = (l + 1) % M;
        for (int i = 1; i < npack; i++) {
            mb.MultB_Left(&emV[l * n], WS);
            l = (l + 1) % M;
        }
        if (jblk != nblk - 1) cblas_dscal(n * n, -1, WS, 1);

        int iblk = (jblk + 1) % nblk;
        double* dest = &G[iblk * n * dim + jblk * n];
        for (int i = 0; i < n; i++) {
            // memcpy(&dest[i * dim], &WS[i * n], n * sizeof(double));

            // also right for npack == M
            cblas_daxpy(n, 1., &WS[i * n], 1, &dest[i * dim], 1);
        }
    }

    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, G, dim, pvt);
    sign = Sgn(dim, G, pvt);
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, dim, G, dim, pvt);
}

void Gtau::Load_G_to_g(int blk0, int blkt) {
    l0 = blk0 * npack, lt = blkt * npack;

    double* src;
    src = &G[blk0 * n * dim + blk0 * n];
    for (int i = 0; i < n; i++) {
        memcpy(&g00[i * n], &src[i * dim], n * sizeof(double));
    }
    for (int i = 0; i < n; i++) g00[i * n + i] -= 1.;

    src = &G[blkt * n * dim + blk0 * n];
    for (int i = 0; i < n; i++) {
        memcpy(&gt0[i * n], &src[i * dim], n * sizeof(double));
    }

    src = &G[blk0 * n * dim + blkt * n];
    for (int i = 0; i < n; i++) {
        memcpy(&g0t[i * n], &src[i * dim], n * sizeof(double));
    }

    src = &G[blkt * n * dim + blkt * n];
    for (int i = 0; i < n; i++) {
        memcpy(&gtt[i * n], &src[i * dim], n * sizeof(double));
    }
    for (int i = 0; i < n; i++) gtt[i * n + i] -= 1.;

    if (blk0 == blkt) {
        for (int i = 0; i < n; i++) g0t[i * n + i] -= 1.;
    }
}

void Gtau::Increase_g_time() {
    mb.MultB_Left(&emV[lt * n], gt0);
    mb.MultBi_Right(&emV[lt * n], g0t);
    mb.MultB_Left(&emV[lt % M * n], gtt);
    mb.MultBi_Right(&emV[lt % M * n], gtt);

    lt = (lt + 1) % M;
    if (lt == 0) {
        cblas_dscal(n * n, -1., gt0, 1);
        cblas_dscal(n * n, -1., g0t, 1);
    }

    if (lt == l0) std::swap(g0t, gt0);
}

void Gtau::Decrease_g_time() {
    if (lt == l0) std::swap(gt0, g0t);

    if (lt == 0) {
        cblas_dscal(n * n, -1., gt0, 1);
        cblas_dscal(n * n, -1., g0t, 1);
    }
    lt = lt == 0 ? M - 1 : lt - 1;

    mb.MultBi_Left(&emV[lt * n], gt0);
    mb.MultB_Right(&emV[lt * n], g0t);
    mb.MultBi_Left(&emV[lt * n], gtt);
    mb.MultB_Right(&emV[lt * n], gtt);
}