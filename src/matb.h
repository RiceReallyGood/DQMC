#ifndef DQMC_MAT_B_H
#define DQMC_MAT_B_H

#include "model.h"

class MatB {
public:
    MatB(const Model *model);
    MatB(const MatB &) = delete;
    MatB &operator=(const MatB &) = delete;
    ~MatB();
    void MultB_Left(const double *emV, double *M);    // M = B * M
    void MultB_Right(const double *emV, double *M);   // M = M * B
    void MultBi_Left(const double *emV, double *M);   // M = B^{-1} * M
    void MultBi_Right(const double *emV, double *M);  // M = M * B^{-1}
    void GetB(const double *emV, double *M);          // M = B = B0 * emV
    void GetBi(const double *emV, double *M);         // M = B^{-1} = emV^{-1} * B0i

private:
    int n;        // Dimension of matrix B
    double *B0;   // Matrix B0 = exp(-dtau * K)
    double *B0i;  // Inverse of matrix B0, B0i = exp(dtau * K)
    double *WS;   // Work space
};

#endif  // DQMC_MAT_B_H