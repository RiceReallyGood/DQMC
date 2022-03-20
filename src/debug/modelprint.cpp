#include "../RingHubbard.h"
#include "../SquareHubbard.h"
using namespace std;

int main() {
    int n = 4;
    double t = 1.;
    double U = 4.;
    double mu = 2.;
    double dtau = 0.1;
    int M = 10;
    Model* mod;
    RingHubbard r(n, t, U, mu, dtau, M);

    double* rrdata = new double[n];
    double* rrdata_err = new double[n];
    std::complex<double>* rcdata = new std::complex<double>[n];
    std::complex<double>* rcdata_err = new std::complex<double>[n];

    double* rrtdata = new double[n * M];
    double* rrtdata_err = new double[n * M];
    std::complex<double>* rctdata = new std::complex<double>[n * M];
    std::complex<double>* rctdata_err = new std::complex<double>[n * M];

    for (int i = 0; i < n; i++) {
        rrdata[i] = rrdata_err[i] = i;
        rcdata[i] = rcdata_err[i] = std::complex<double>(i, i);
        for (int l = 0; l < M; l++) {
            rrtdata[l * n + i] = rrtdata_err[l * n + i] = l * n + i;
            rctdata[l * n + i] = rctdata_err[l * n + i] =
                std::complex<double>(l * n + i, l * n + i);
        }
    }

    mod = &r;
    mod->Print_ETM("rdata", rrdata, rrdata_err);
    std::cout << "=======================================================" << std::endl;
    mod->Print_ETM_FT("cdata", rcdata, rcdata_err);
    std::cout << "=======================================================" << std::endl;
    mod->Print_TDM("rtdata", rrtdata, rrtdata_err);
    std::cout << "=======================================================" << std::endl;
    mod->Print_TDM_FT("ctdata", rctdata, rctdata_err);

    SquareHubbard s(n, t, 0, 0, U, mu, dtau, M);
    double* srdata = new double[n * n];
    double* srdata_err = new double[n * n];
    std::complex<double>* scdata = new std::complex<double>[n * n];
    std::complex<double>* scdata_err = new std::complex<double>[n * n];

    double* srtdata = new double[n * n * M];
    double* srtdata_err = new double[n * n * M];
    std::complex<double>* sctdata = new std::complex<double>[n * n * M];
    std::complex<double>* sctdata_err = new std::complex<double>[n * n * M];

    for (int i = 0; i < n * n; i++) {
        srdata[i] = srdata_err[i] = i;
        scdata[i] = scdata_err[i] = std::complex<double>(i, i);
        for (int l = 0; l < M; l++) {
            srtdata[l * n * n + i] = srtdata_err[l * n * n + i] = l * n * n + i;
            sctdata[l * n * n + i] = sctdata_err[l * n * n + i] =
                std::complex<double>(l * n * n + i, l * n * n + i);
        }
    }

    mod = &s;
    mod->Print_ETM("srdata", srdata, srdata_err);
    std::cout
        << "=========================================================================="
        << std::endl;
    mod->Print_ETM_FT("scdata", scdata, scdata_err);
    std::cout
        << "=========================================================================="
        << std::endl;
    mod->Print_TDM("srtdata", srtdata, srtdata_err);
    std::cout
        << "=========================================================================="
        << std::endl;
    mod->Print_TDM_FT("sctdata", sctdata, sctdata_err);

    delete[] sctdata_err;
    delete[] sctdata;
    delete[] scdata_err;
    delete[] scdata;
    delete[] srtdata_err;
    delete[] srtdata;
    delete[] srdata_err;
    delete[] srdata;

    delete[] rctdata_err;
    delete[] rctdata;
    delete[] rcdata_err;
    delete[] rcdata;
    delete[] rrtdata_err;
    delete[] rrtdata;
    delete[] rrdata_err;
    delete[] rrdata;
}