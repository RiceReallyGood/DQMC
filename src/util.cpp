#include "util.h"

#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>

#include "mkl.h"

void ScaleRows(int n, double *A, const double *D) {
    for (int i = 0; i < n; i++) {
        cblas_dscal(n, D[i], &A[i * n], 1);
    }
}

void ScaleCols(int n, double *A, const double *D) {
    for (int i = 0; i < n; i++) {
        cblas_dscal(n, D[i], &A[i], n);
    }
}

void ScaleRowsInv(int n, double *A, const double *D) {
    for (int i = 0; i < n; i++) {
        cblas_dscal(n, 1. / D[i], &A[i * n], 1);
    }
}

void ScaleColsInv(int n, double *A, const double *D) {
    for (int i = 0; i < n; i++) {
        cblas_dscal(n, 1. / D[i], &A[i], n);
    }
}

void Get_Avrg_Err(int nBins, int ndata, const double *data, double *avrg, double *err) {
    memset(avrg, 0, ndata * sizeof(double));
    for (int iBin = 0; iBin < nBins; iBin++) {
        vdAdd(ndata, avrg, &data[iBin * ndata], avrg);
    }
    cblas_dscal(ndata, 1. / double(nBins), avrg, 1);

    double *temp = new double[ndata];
    memset(err, 0, ndata * sizeof(double));
    for (int iBin = 0; iBin < nBins; iBin++) {
        vdSub(ndata, &data[iBin * ndata], avrg, temp);
        vdMul(ndata, temp, temp, temp);
        vdAdd(ndata, err, temp, err);
    }

    cblas_dscal(ndata, 1. / (double(nBins) * double(nBins - 1)), err, 1);

    vdSqrt(ndata, err, err);

    delete[] temp;
}

void Eye(int n, double *A) {
    memset(A, 0, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        A[i * n + i] = 1.;
    }
}

void randmat(int n, double *A, std::default_random_engine &dre,
             std::uniform_real_distribution<double> &urd) {
    for (int i = 0; i < n * n; i++) A[i] = urd(dre);
}

void printmat(const std::string &name, int m, int n, int lda, const double *mat,
              int precision_, int floatfield_) {
    std::ios::fmtflags OldFlags = std::cout.flags();
    std::cout << name << " :" << std::endl;

    std::cout.precision(precision_);
    int width = 0;
    if (floatfield_ == FIXED) {
        std::cout.setf(std::ios::fixed);
        width = precision_ + 3;
    }
    else if (floatfield_ == SCIENTIFIC) {
        std::cout.setf(std::ios::scientific);
        width = precision_ + 8;
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << std::setw(width) << mat[i * lda + j];
            if (j != n - 1) std::cout << "   ";
        }
        std::cout << std::endl;
    }

    std::cout.flags(OldFlags);
    std::cout.precision(6);
}

double MatDiff(int n, const double *A, const double *B) {
    double diff = 0.;
    for (int i = 0; i < n * n; i++) {
        diff = std::max(diff, std::abs(A[i] - B[i]));
    }
    return diff;
}

void Trans(int n, double *dest, double *src) {
    for (int i = 0; i < n; i++) {
        cblas_dcopy(n, &src[i * n], 1, &dest[i], n);
    }
}

int Sgn(int n, double *LUf, int *pvt) {
    int s = 1;
    for (int i = 0; i < n; i++) {
        if (pvt[i] != i + 1) s = -s;
        if (LUf[i * n + i] < 0.) s = -s;
    }
    return s;
}

double LogDet(int n, double *LUf) {
    double d = 0.;
    for (int i = 0; i < n; i++) {
        d += std::log(std::abs(LUf[i * n + i]));
    }
    return d;
}