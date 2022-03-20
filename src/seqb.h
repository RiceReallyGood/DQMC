#ifndef DQMC_SEQB_H
#define DQMC_SEQB_H

#include "matb.h"

class SeqB {
public:
    SeqB(int n_, int M_, int npack_);
    SeqB(const SeqB &) = delete;
    SeqB &operator=(const SeqB &) = delete;
    ~SeqB();
    double *U;
    double *D;
    double *T;

    void Update_UDT_Form();
    void Get_SeqMultB(int ll, int lr, MatB &mb, const double *emV);
    void Get_SeqMultBi(int ll, int lr, MatB &mb, const double *emV);

private:
    int n;        // Dimension of the matrix
    int M;        // Number of time slices
    int npack;    //
    int *jpvt;    // working array in QR factor
    double *tau;  // working array in QR factor
    double *WS;   // working space for matrix multiplication

    void reset();
};

#endif  // DQMC_SEQB_H