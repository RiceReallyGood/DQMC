#ifndef DQMC_GFUN_H
#define DQMC_GFUN_H

#include "matb.h"
#include "model.h"
#include "seqb.h"

class Gfun {
public:
    Gfun(const Model *model, double *&emV_, int nWrap_, int npack_, double toldiff_);
    Gfun(const Gfun &) = delete;
    Gfun &operator=(const Gfun &) = delete;
    ~Gfun();

    int sign;       // Sign of det(G)
    double logdet;  // Log of Abs(det(G))
    double *G;      // Matrix of Green's function
    int l;          // Index of right most B

    void GetG(int lr);
    void ComputeG(int lr);
    void UpdateG(int j, double gamma, double r);
    void ApplyUpdate();
    double GetGjj(int j);

private:
    int n;           // Number of sites
    int M;           // Number of time slices
    double *(&emV);  // Matrix info of emV(0)...emV(M - 1)

    // For numerical stabilization
    int nWrap;       // Number of max Warps
    int wps;         // Number of warps can be done before next check
    double toldiff;  // Tolerable matrix diff
    MatB mb;
    SeqB sb;

    // For delayed update
    int ndelay;
    double *V;  // V is used as workspace when GetG() is called
    double *W;  // W is also used as workspace when Compute G from scratch

    // Work space
    double *D2;
    int *pvt1, *pvt2;
};

#endif  // DQMC_GFUN_H