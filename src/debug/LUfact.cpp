#include <iostream>
#include <random>

#include "../util.h"
#include "mkl.h"

int main() {
    int n = 8;
    std::default_random_engine dre(time(nullptr));
    std::uniform_real_distribution<double> urd(0, 1);
    double* mat = new double[n * n];
    randmat(n, mat, dre, urd);

    int* pvt = new int[n];
    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, mat, n, pvt);

    for(int i = 0; i < n; i++)
        std::cout << pvt[i] << std::endl;
}