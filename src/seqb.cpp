#include "seqb.h"

#include <cstring>
#include <iostream>

#include "mkl.h"
#include "util.h"

SeqB::SeqB(int n_, int M_, int npack_) : n(n_), M(M_), npack(npack_) {
    U = new double[n * n];
    D = new double[n];
    T = new double[n * n];
    jpvt = new int[n];
    tau = new double[n];
    WS = new double[n * n];
    reset();
}

SeqB::~SeqB() {
    delete[] WS;
    delete[] tau;
    delete[] jpvt;
    delete[] T;
    delete[] D;
    delete[] U;
}

void SeqB::reset() {
    Eye(n, U);
    Eye(n, T);
    for (int i = 0; i < n; i++) D[i] = 1.;
}

void SeqB::Update_UDT_Form() {
    ScaleCols(n, U, D);
    memset(jpvt, 0, n * sizeof(int));
    LAPACKE_dgeqp3(LAPACK_ROW_MAJOR, n, n, U, n, jpvt, tau);

    for (int i = 0; i < n; i++) {
        D[i] = U[i * n + i];
        cblas_dscal(n - i, 1. / D[i], &U[i * n + i], 1);
    }

    for (int i = 0; i < n; i++) {
        memcpy(&WS[i * n], &T[(jpvt[i] - 1) * n], n * sizeof(double));
    }

    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasUnit, n, n, 1.,
                U, n, WS, n);
    memcpy(T, WS, n * n * sizeof(double));

    LAPACKE_dorgqr(LAPACK_ROW_MAJOR, n, n, n, U, n, tau);
}

void SeqB::Get_SeqMultB(int ll, int lr, MatB &mb, const double *emV) {
    // This subroutine computes A = B_{ll}B_{ll-1}...B_{lr}
    // and returns A's UDT decomposition.
    int interval = ll >= lr ? ll - lr + 1 : ll + M - lr + 1;
    int l = lr;

    reset();
    for (int i = 0; i < interval; i++) {
        mb.MultB_Left(&emV[l * n], U);

        if ((i + 1) % npack == 0) Update_UDT_Form();

        l = (l + 1) % M;
    }

    if (interval % npack != 0) Update_UDT_Form();
}

void SeqB::Get_SeqMultBi(int ll, int lr, MatB &mb, const double *emV) {
    // This functions computes A = inv(B_{ll}B_{ll-1}...B_{lr})
    // and returns A's UDT decomposition.
    int interval = ll >= lr ? ll - lr + 1 : ll + M - lr + 1;
    int l = lr;

    reset();
    for (int i = 0; i < interval; i++) {
        mb.MultBi_Right(&emV[l * n], U);

        if ((i + 1) % npack == 0) Update_UDT_Form();

        l = (l + 1) % M;
    }

    if (interval % npack != 0) Update_UDT_Form();
}