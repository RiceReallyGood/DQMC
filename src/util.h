#ifndef DQMC_UTIL_H
#define DQMC_UTIL_H

#include <random>
#include <string>
#define FIXED 0
#define SCIENTIFIC 1

void ScaleRows(int n, double *A, const double *D);     // A = D * A
void ScaleCols(int n, double *A, const double *D);     // A = A * D
void ScaleRowsInv(int n, double *A, const double *D);  // A = Inv(D) * A
void ScaleColsInv(int n, double *A, const double *D);  // A = A * Inv(D)
void Get_Avrg_Err(int nBin, int ndata, const double *data, double *avrg, double *err);
void Eye(int n, double *A);  // Identity Matrix
inline double sq(double x) { return x * x; }
void printmat(const std::string &name, int m, int n, int lda, const double *mat,
              int precision_ = 6, int floatfield_ = FIXED);
void randmat(int n, double *A, std::default_random_engine &dre,
             std::uniform_real_distribution<double> &urd);
double MatDiff(int n, const double *A, const double *B);
void Trans(int n, double *dest, double *src);
int Sgn(int n, double *LUf, int *pvt);
double LogDet(int n, double *LUf);
#endif  // DQMC_UTIL_H
