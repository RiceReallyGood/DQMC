#ifndef DQMC_MEASURE_H
#define DQMC_MEASURE_H

#include "gtau.h"

// Time dependent measurements
class TDM {
public:
    TDM(int n, int M, int npack);
    TDM(const TDM&) = delete;
    TDM& operator=(const TDM&) = delete;
    ~TDM();

    double sign;

    double* gtup;
    double* gtdn;

    double* SzSz;
    double* denden;

    void measure(Gtau& gt_up, Gtau& gt_dn);
    void finalize();
    void clear();

private:
    int n;
    int M;
    int npack;
    int nblk;
    int inc;
    int dec;
    int count;

    double* nunu;
    double* nund;
    double* ndnd;

    void record(Gtau& gt_up, Gtau& gt_dn);
};

#endif  // DQMC_MEASURE_H