#include <cstring>

#include "../util.h"
#include "mkl.h"
using namespace std;

int main() {
    int n = 4;
    double *A = new double[n * n];
    double *ExpA = new double[n * n];

    memset(A, 0, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        A[i * n + (i - 1 + n) % n] = 1.;
        A[i * n + (i + 1) % n] = 1.;
    }

    printmat("Matrix A", n, n, n, A);

    double *w = new double[n];
    LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', n, A, n, w);

    printmat("\nEigenvalues", 1, n, n, w);
    printmat("\nEigenvectors", n, n, n, A);

    double *TempMat = new double[n * n];
    memcpy(TempMat, A, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        cblas_dscal(n, std::exp(w[i]), &TempMat[i], n);
    }
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1., TempMat, n, A, n,
                0., ExpA, n);

    printmat("\nExp A", n, n, n, ExpA);

    delete[] TempMat;
    delete[] w;
    delete[] ExpA;
    delete[] A;
}