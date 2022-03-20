#include <cstring>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "../RingHubbard.h"
#include "../gfun.h"
#include "../matb.h"
#include "../seqb.h"
#include "../util.h"
#include "mkl.h"

int main(int argc, char *argv[]) {
    int n = 8, M = 160, npack = 1;
    double t = 1., U = 4, mu = 0., dtau = 0.05;
    RingHubbard r(n, t, U, mu, dtau, M);
    std::default_random_engine dre(time(NULL));
    std::uniform_real_distribution<double> urd(0, 1);
    SeqB A(n, M, npack);

    randmat(n, A.U, dre, urd);

    printmat("randmat", n, n, n, A.U);

    A.Update_UDT_Form();

    ScaleCols(n, A.U, A.D);
    double *prod = new double[n * n];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., A.U, n, A.T, n,
                0, prod, n);

    printmat("U * D * T", n, n, n, prod);

    MatB B(&r);
    double *emV = new double[M * n];
    double temp = std::exp(dtau * U / 2);
    double lambda = std::log(temp + std::sqrt(temp * temp - 1));
    std::cout << "lambda = " << lambda << std::endl;
    // std::uniform_int_distribution<int> uid(0, 1);
    for (int i = 0; i < M * n; i++) {
        emV[i] = std::exp(-(2 * ((i % M) % 2) - 1) * lambda);
    }

    // get sequential multiplication of matrix b by UDT
    A.Get_SeqMultB(M - 1, 0, B, emV);
    ScaleCols(n, A.U, A.D);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., A.U, n, A.T, n,
                0, prod, n);
    printmat("U * D * T of seqB", n, n, n, prod, 15, SCIENTIFIC);

    double *WS1 = new double[n * n];
    double *WS2 = new double[n * n];

    Eye(n, prod);
    for (int l = 0; l < M; l++) {
        B.GetB(&emV[l * n], WS1);
        memcpy(WS2, prod, n * n * sizeof(double));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1., WS1, n, WS2,
                    n, 0, prod, n);
    }
    printmat("sequential Multiplication of seqB", n, n, n, prod, 15, SCIENTIFIC);

    int nWrap = 15;
    Gfun gf(&r, emV, nWrap, npack, 1e-6);

    gf.ComputeG(0);
    printmat("G from ComputeG()", n, n, n, gf.G);

    std::cout << "Log(Abs(Det(G))): " << gf.logdet << std::endl;
    std::cout << "Abs(Det(G)): " << std::exp(gf.logdet) << std::endl;
    std::cout << "Sign(Det(G)): " << gf.sign << std::endl;

    for (int i = 0; i < n; i++) prod[i * n + i] += 1.;

    int *pvt = new int[n];
    double ldet = 1.;
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, prod, n, pvt);
    for (int i = 0; i < n; i++) ldet -= std::log(std::abs(prod[i * n + i]));
    LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, prod, n, pvt);

    printmat("inverse of 1 + seqB", n, n, n, prod);
    std::cout << "LogDet of inverse of 1 + seqB: " << ldet << std::endl;

    // update time slice 0
    for (int j = 0; j < n; j += 2) {
        double gamma = std::exp(-2. * lambda) - 1.;
        double r = 1. + (1. - gf.GetGjj(j)) * gamma;
        gf.UpdateG(j, gamma, r);
        emV[j] = 1. / emV[j];
    }

    for (int j = 0; j < n; j++) {
        prod[j] = gf.GetGjj(j);
    }

    printmat("\nGjj", 1, n, n, prod);
    gf.ApplyUpdate();
    printmat("\nUpdated G", n, n, n, gf.G);
    std::cout << "Updated logdet: " << gf.logdet << std::endl;
    std::cout << "Updated sign: " << gf.sign << std::endl;

    gf.ComputeG(0);
    printmat("\nComputed G", n, n, n, gf.G);
    std::cout << "Computed logdet: " << gf.logdet << std::endl;
    std::cout << "Computed sign: " << gf.sign << std::endl;

    std::ostringstream s;
    for (int l = 1; l < nWrap; l++) {
        gf.GetG(l);
        s.str("");
        s << "\nWraped G(" << l << ")";
        printmat(s.str(), n, n, n, gf.G);
    }
    gf.ComputeG(nWrap - 1);
    s.str("");
    s << "\nComputed G(" << nWrap - 1 << ")";
    printmat(s.str(), n, n, n, gf.G);

    delete[] WS1;
    delete[] WS2;
    delete[] emV;
    delete[] prod;
}