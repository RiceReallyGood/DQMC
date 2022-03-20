#include <random>

#include "../RingHubbard.h"
#include "../matb.h"
#include "../util.h"
#include "mkl.h"

int main() {
    int n = 4;
    double t = 1.;
    double U = 4.;
    double mu = 2.;
    double dtau = 0.1;
    int M = 10;
    RingHubbard r(n, t, U, mu, dtau, M);

    MatB mb(&r);
    double* matB = new double[n * n];
    double* matBi = new double[n * n];
    double* WS = new double[n * n];
    double* emV = new double[n];
    double temp = std::exp(dtau * U / 2);
    double lambda = std::log(temp + std::sqrt(temp * temp - 1));
    for (int i = 0; i < n; i++) {
        emV[i] = std::exp(-(i % 2 * 2 - 1) * lambda);
    }

    mb.GetB(emV, matB);
    printmat("Matrix B:", n, n, n, matB);
    mb.GetBi(emV, matBi);
    printmat("\nInverse of Matrix B", n, n, n, matBi);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., matB, n, matBi,
                n, 0, WS, n);
    printmat("\nmatB * matBi", n, n, n, WS);

    std::default_random_engine dre(time(NULL));
    std::uniform_real_distribution<double> urd(0, 1);

    double* rm = new double[n * n];
    for (int i = 0; i < n * n; i++) rm[i] = urd(dre);
    randmat(n, rm, dre, urd);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., matB, n, rm, n,
                0, WS, n);
    printmat("\nB * randmat", n, n, n, WS);

    mb.MultB_Left(emV, rm);
    printmat("\n call mb.MultB_Left(emV, randmat)", n, n, n, rm);

    randmat(n, rm, dre, urd);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., rm, n, matBi, n,
                0, WS, n);
    printmat("\nrandmat * Bi", n, n, n, WS);

    mb.MultBi_Right(emV, rm);
    printmat("\nCall mb.MultBi_Right(emV, randmat)", n, n, n, rm);

    randmat(n, rm, dre, urd);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., rm, n, matB, n,
                0, WS, n);
    printmat("\nrandmat * B", n, n, n, WS);

    mb.MultB_Right(emV, rm);
    printmat("\nCall mb.MultB_Right(emV, randmat)", n, n, n, rm);

    randmat(n, rm, dre, urd);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., matBi, n, rm, n,
                0, WS, n);
    printmat("\nCall mb.MultB_Right(emV, randmat)", n, n, n, WS);

    mb.MultBi_Left(emV, rm);
    printmat("\nCall mb.MultBi_Left(emV, randmat)", n, n, n, rm);

    delete[] rm;
    delete[] WS;
    delete[] emV;
    delete[] matBi;
    delete[] matB;
}
