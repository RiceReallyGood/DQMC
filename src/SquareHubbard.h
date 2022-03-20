#ifndef DQMC_SQUARE_HUBBARD_H
#define DQMC_SQUARE_HUBBARD_H

#include "model.h"

class SquareHubbard : public Model {
public:
    SquareHubbard(int n_, double t_, double tp_, double tpp_, double U_, double mu_,
                  double dtau_, int M_);
    SquareHubbard(const SquareHubbard &) = delete;
    SquareHubbard &operator=(const SquareHubbard &) = delete;
    ~SquareHubbard();

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
    int n;                            // Lattice in each direction
    std::complex<double> *FT_Factor;  // Fourier tranform factor
};

#endif  // DQMC_SQUARE_HUBBARD_H