#ifndef DQMC_H
#define DQMC_H

#include <iostream>
#include <random>

#include "gfun.h"
#include "gtau.h"
#include "measure.h"
#include "model.h"

class DQMC {
public:
    DQMC(const Model* model, int nWarm, int nBins, int BinSize, int MeasFreq, int nWrap,
         int npack, double toldiff, std::ostream& os);
    DQMC(const DQMC&) = delete;
    DQMC& operator=(const DQMC&) = delete;
    ~DQMC();

    void PrintProfile();
    void Simulate();
    void PrintResult();

private:
    const Model* model;
    int M;
    int n;
    int nWarm;
    int nBins;
    int BinSize;
    int MeasFreq;
    std::ostream& os;

    int* HSF;
    double* explookup;
    double* emV_up;
    double* emV_dn;
    TDM tdm;
    Gtau gt_up;
    Gtau gt_dn;
    Gfun gf_up;
    Gfun gf_dn;

    double sign;
    double *giup, *giup_err;
    double *gidn, *gidn_err;
    std::complex<double>*gkup, *gkup_err;
    std::complex<double>*gkdn, *gkdn_err;
    double *giavrg, *giavrg_err;
    std::complex<double>*gkavrg, *gkavrg_err;
    double *SzSzi, *SzSzi_err;
    std::complex<double>*SzSzk, *SzSzk_err;
    double *dendeni, *dendeni_err;
    std::complex<double>*dendenk, *dendenk_err;

    std::mt19937 e;  // mersenne_twister_engine
    std::uniform_int_distribution<int>
        uid;  // uniform int distribution in {0, 1, 2, ..., M - 1}
    std::uniform_real_distribution<double>
        urd;  // uniform real distribution in [0., 1.)

    void Sweep();
    void Measure();
    void GetTrans(int iBin);
    void GetFT(int iBin);
    void Get_Statistics();
};

#endif