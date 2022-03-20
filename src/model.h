#ifndef DQMC_MODEL_H
#define DQMC_MODEL_H

#include <complex>
#include <iostream>
#include <string>

class CSR_Mat {
public:
    int n;       // dimension of matrix
    int nnz;     // number of nonzeros
    int *value;  // nonzeros elements, dim = nnz
    int *col;    // col indices of nonzeros elements, dim = nnz
    int *start;  // start[i] is the index of the element in value array
                 // that is first non-zero element in a row j of sparse matrix
                 // dim = n + 1

    CSR_Mat(int dim_, int nnz_) : n(dim_), nnz(nnz_) {
        value = new int[nnz];
        col = new int[nnz];
        start = new int[n + 1];
    }

    ~CSR_Mat() {
        delete[] start;
        delete[] col;
        delete[] value;
    }
};

class Model {
public:
    Model(int nsites_, int n_t_, int nnz_, double U_, double mu_, double dtau_, int M_)
        : nsites(nsites_),
          n_t(n_t_),
          Adj(nsites_, nnz_),
          U(U_),
          mu(mu_),
          dtau(dtau_),
          M(M_) {
        t = new double[n_t];
    }
    Model(const Model &) = delete;
    Model &operator=(const Model &) = delete;
    virtual ~Model() { delete[] t; }

    int nsites;   // Number of sites
    int n_t;      // Number of hopping types
    double *t;    // Hopping energy of different types
    CSR_Mat Adj;  // Adjacent Matrix T in CSR Storage
                  // T[i,j] = hopping type
    double U;     // Onsite interaction
    double mu;    // Chemical potential
    double dtau;  // Length of timeslices
    int M;        // Number of timeslices
    virtual void Print_Parms(std::ostream &os = std::cout) const = 0;
    virtual void Trans_Avrg(double *datamat, double *datavec) const = 0;
    virtual void FT(double *data, std::complex<double> *ftdata) const = 0;
    virtual void Print_ETM(const std::string &name, double *data, double *data_err,
                           std::ostream &os = std::cout) const = 0;
    virtual void Print_TDM(const std::string &name, double *data, double *data_err,
                           std::ostream &os = std::cout) const = 0;
    virtual void Print_ETM_FT(const std::string &name, std::complex<double> *data,
                              std::complex<double> *data_err,
                              std::ostream &os = std::cout) const = 0;
    virtual void Print_TDM_FT(const std::string &name, std::complex<double> *data,
                              std::complex<double> *data_err,
                              std::ostream &os = std::cout) const = 0;
};

#endif  // DQMC_MODEL_H