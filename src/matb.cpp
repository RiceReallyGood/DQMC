#include "matb.h"

#include <cmath>
#include <cstring>

#include "mkl.h"
#include "util.h"

MatB::MatB(const Model* model) : n(model->nsites) {
    B0 = new double[n * n];
    B0i = new double[n * n];
    WS = new double[n * n];

    // alias
    const int* const start = model->Adj.start;
    const int* const col = model->Adj.col;
    const int* const kind = model->Adj.value;
    const double* const t = model->t;
    const double& dtau = model->dtau;
    const double& mu = model->mu;

    memset(WS, 0, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = start[i]; j < start[i + 1]; j++) {
            WS[i * n + col[j]] += dtau * t[kind[j]];
        }
        WS[i * n + i] += dtau * mu;
    }

    double* w = new double[n];
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, WS, n, w);

    double* TempMat = new double[n * n];
    memcpy(TempMat, WS, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        cblas_dscal(n, std::exp(w[i]), &TempMat[i], n);
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1., TempMat, n, WS, n,
                0, B0, n);

    memcpy(TempMat, WS, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        cblas_dscal(n, std::exp(-w[i]), &TempMat[i], n);
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1., TempMat, n, WS, n,
                0, B0i, n);

    // printmat("matb", n, n, n, B0);

    delete[] TempMat;
    delete[] w;
}

MatB::~MatB() {
    delete[] WS;
    delete[] B0i;
    delete[] B0;
}

void MatB::MultB_Left(const double* emV, double* M) {
    memcpy(WS, M, n * n * sizeof(double));
    ScaleRows(n, WS, emV);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., B0, n, WS, n, 0,
                M, n);
}

void MatB::MultB_Right(const double* emV, double* M) {
    memcpy(WS, M, n * n * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., WS, n, B0, n, 0,
                M, n);
    ScaleCols(n, M, emV);
}

void MatB::MultBi_Left(const double* emV, double* M) {
    memcpy(WS, M, n * n * sizeof(double));
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., B0i, n, WS, n,
                0, M, n);
    ScaleRowsInv(n, M, emV);
}

void MatB::MultBi_Right(const double* emV, double* M) {
    memcpy(WS, M, n * n * sizeof(double));
    ScaleColsInv(n, WS, emV);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., WS, n, B0i, n,
                0, M, n);
}

void MatB::GetB(const double* emV, double* M) {
    memcpy(M, B0, n * n * sizeof(double));
    ScaleCols(n, M, emV);
}

void MatB::GetBi(const double* emV, double* M) {
    memcpy(M, B0i, n * n * sizeof(double));
    ScaleRowsInv(n, M, emV);
}