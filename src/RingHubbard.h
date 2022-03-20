#ifndef DQMC_RING_HUBBARD_H
#define DQMC_RING_HUBBARD_H

#include "model.h"

class RingHubbard : public Model {
public:
    RingHubbard(int n_, double t_, double U_, double mu_, double dtau_, int M_);
    RingHubbard(const RingHubbard &) = delete;
    RingHubbard &operator=(const RingHubbard &) = delete;
    ~RingHubbard();

    virtual void Print_Parms(std::ostream &os = std::cout) const override;
    virtual void Trans_Avrg(double *datamat, double *datavec) const override;
    virtual void FT(double *data, std::complex<double> *ftdata) const override;
    virtual void Print_ETM(const std::string &name, double *data, double *data_err,
                           std::ostream &os = std::cout) const override;
    virtual void Print_TDM(const std::string &name, double *data, double *data_err,
                           std::ostream &os = std::cout) const override;
    virtual void Print_ETM_FT(const std::string &name, std::complex<double> *data,
                              std::complex<double> *data_err,
                              std::ostream &os = std::cout) const override;
    virtual void Print_TDM_FT(const std::string &name, std::complex<double> *data,
                              std::complex<double> *data_err,
                              std::ostream &os = std::cout) const override;
    int n;
    std::complex<double> *FT_Factor;  // Fourier tranform factor
};

#endif  // DQMC_RING_HUBBARD_H