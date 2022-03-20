#include "RingHubbard.h"

#include <iomanip>

RingHubbard::RingHubbard(int n_, double t_, double U_, double mu_, double dtau_, int M_)
    : Model(n_, 1, 2 * n_, U_, mu_, dtau_, M_), n(n_) {
    t[0] = t_;

    // Fill in adjcent Matrix
    for (int i = 0; i < n; i++) {
        Adj.start[i] = 2 * i;
        Adj.col[2 * i] = (i - 1 + n) % n;
        Adj.col[2 * i + 1] = (i + 1) % n;
        Adj.value[2 * i] = 0;
        Adj.value[2 * i + 1] = 0;
    }
    Adj.start[n] = 2 * n;

    const double pi = 3.1415926535897932385;
    const std::complex<double> I(0, 1);
    FT_Factor = new std::complex<double>[n];
    for (int i = 0; i < n; i++) {
        FT_Factor[i] = std::exp(-I * 2. * pi * double(i) / double(n));
    }
}

RingHubbard::~RingHubbard() { delete[] FT_Factor; }

void RingHubbard::Print_Parms(std::ostream &os) const {
    // save current format flags
    std::ios::fmtflags oldFlags = os.flags();
    os.setf(std::ios::fixed);

    os << "                         Model:  1D Ring Hubbard" << std::endl;
    os << "                             n:  " << n << std::endl;
    os << "               Number of sites:  " << n << std::endl;
    os << "                             t:  " << t[0] << std::endl;
    os << "                             U:  " << U << std::endl;
    os << "                            mu:  " << mu << std::endl;
    os << "                          dtau:  " << dtau << std::endl;
    os << "         Number of time slices:  " << M << std::endl;
    os << "                          beta:  " << M * dtau << std::endl;

    // restore saved format flags
    os.flags(oldFlags);
}

void RingHubbard::Trans_Avrg(double *datamat, double *datavec) const {
    for (int d = 0; d < n; d++) {
        datavec[d] = 0.;
        for (int i = 0; i < n; i++) {
            datavec[d] += datamat[i * n + (i + d) % n];
        }
        datavec[d] /= n;
    }
}

void RingHubbard::FT(double *data, std::complex<double> *ftdata) const {
    for (int k = 0; k < n; k++) {
        ftdata[k] = 0.;
        for (int i = 0; i < n; i++) {
            ftdata[k] += data[i] * FT_Factor[(i * k) % n];
        }
        // ftdata[k] /= double(n);
    }
}

void RingHubbard::Print_ETM(const std::string &name, double *data, double *data_err,
                            std::ostream &os) const {
    // save current format flags
    std::ios::fmtflags oldFlags = os.flags();

    os.setf(std::ios::scientific);

    os << name << ":" << std::endl;
    for (int i = 0; i < n; i++) {
        os << std::setw(3) << i << "    " << std::setw(13) << data[i] << " +- "
           << std::setw(12) << data_err[i] << std::endl;
    }

    // restore saved format flags
    os.flags(oldFlags);
}

void RingHubbard::Print_TDM(const std::string &name, double *data, double *data_err,
                            std::ostream &os) const {
    // save current format flags
    std::ios::fmtflags oldFlags = os.flags();

    for (int i = 0; i < n; i++) {
        os << name << " i = " << i << ":" << std::endl;
        for (int l = 0; l < M; l++) {
            os << std::setiosflags(std::ios::fixed) << std::setw(10) << l * dtau
               << "    " << std::resetiosflags(std::ios::floatfield)
               << std::setiosflags(std::ios::scientific) << std::setw(13)
               << data[l * n + i] << " +- " << std::setw(12) << data_err[l * n + i]
               << std::resetiosflags(std::ios::floatfield) << std::endl;
        }
        os << std::endl;
    }

    // restore saved format flags
    os.flags(oldFlags);
}

void RingHubbard::Print_ETM_FT(const std::string &name, std::complex<double> *data,
                               std::complex<double> *data_err, std::ostream &os) const {
    // save current format flags
    std::ios::fmtflags oldFlags = os.flags();

    os.setf(std::ios::scientific);

    os << name << ":" << std::endl;
    for (int k = 0; k < n; k++) {
        os << std::setw(3) << k << "    (" << std::setw(13) << data[k].real() << " +- "
           << std::setw(12) << data_err[k].real() << ") + I(" << std::setw(13)
           << data[k].imag() << " +- " << std::setw(12) << data_err[k].imag() << ")"
           << std::endl;
    }

    // restore saved format flags
    os.flags(oldFlags);
}

void RingHubbard::Print_TDM_FT(const std::string &name, std::complex<double> *data,
                               std::complex<double> *data_err, std::ostream &os) const {
    // save current format flags
    std::ios::fmtflags oldFlags = os.flags();

    for (int k = 0; k < n; k++) {
        os << name << " k = " << k << ":" << std::endl;
        for (int l = 0; l < M; l++) {
            os << std::setiosflags(std::ios::fixed) << std::setw(10) << l * dtau
               << "    (" << std::resetiosflags(std::ios::floatfield)
               << std::setiosflags(std::ios::scientific) << std::setw(13)
               << data[l * n + k].real() << " +- " << std::setw(12)
               << data_err[l * n + k].real() << ") + I(" << std::setw(13)
               << data[l * n + k].imag() << " +- " << std::setw(12)
               << data_err[l * n + k].imag() << ")"
               << std::resetiosflags(std::ios::floatfield) << std::endl;
        }
        os << std::endl;
    }

    // restore saved format flags
    os.flags(oldFlags);
}