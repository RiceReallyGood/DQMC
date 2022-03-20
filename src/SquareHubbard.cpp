#include "SquareHubbard.h"

#include <iomanip>

SquareHubbard::SquareHubbard(int n_, double t_, double tp_, double tpp_, double U_,
                             double mu_, double dtau_, int M_)
    : Model(n_ * n_, 3, 12 * n_ * n_, U_, mu_, dtau_, M_), n(n_) {
    t[0] = t_;
    t[1] = tp_;
    t[2] = tpp_;

    // Fill in adjcent Matrix
    for (int i = 0; i < nsites; i++) {
        Adj.start[i] = 12 * i;
        int x = i / n, y = i % n;
        int l = (x - 1 + n) % n, r = (x + 1) % n;
        int u = (y + 1) % n, d = (y - 1 + n) % n;
        int ll = (x - 2 + n) % n, rr = (x + 2) % n;
        int uu = (y + 2) % n, dd = (y - 2 + n) % n;
        // left
        Adj.col[12 * i + 0] = l * n + y;
        Adj.value[12 * i + 0] = 0;
        // right
        Adj.col[12 * i + 1] = r * n + y;
        Adj.value[12 * i + 1] = 0;
        // up
        Adj.col[12 * i + 2] = x * n + u;
        Adj.value[12 * i + 2] = 0;
        // down
        Adj.col[12 * i + 3] = x * n + d;
        Adj.value[12 * i + 3] = 0;

        // left-up
        Adj.col[12 * i + 4] = l * n + u;
        Adj.value[12 * i + 4] = 1;
        // right-up
        Adj.col[12 * i + 5] = r * n + u;
        Adj.value[12 * i + 5] = 1;
        // right-down
        Adj.col[12 * i + 6] = r * n + d;
        Adj.value[12 * i + 6] = 1;
        // left-down
        Adj.col[12 * i + 7] = l * n + d;
        Adj.value[12 * i + 7] = 1;

        // left-left
        Adj.col[12 * i + 8] = ll * n + y;
        Adj.value[12 * i + 8] = 2;
        // right-right
        Adj.col[12 * i + 9] = rr * n + y;
        Adj.value[12 * i + 9] = 2;
        // up-up
        Adj.col[12 * i + 10] = x * n + uu;
        Adj.value[12 * i + 10] = 2;
        // down-down
        Adj.col[12 * i + 11] = x * n + dd;
        Adj.value[12 * i + 11] = 2;
    }
    Adj.start[nsites] = 12 * nsites;

    const double pi = 3.1415926535897932385;
    const std::complex<double> I(0, 1);
    FT_Factor = new std::complex<double>[n];
    for (int i = 0; i < n; i++) {
        FT_Factor[i] = std::exp(-I * 2. * pi * double(i) / double(n));
    }
}

SquareHubbard::~SquareHubbard() { delete[] FT_Factor; }

void SquareHubbard::Print_Parms(std::ostream &os) const {
    // save current format flags
    std::ios::fmtflags oldFlags = os.flags();
    os.setf(std::ios::fixed);

    os << "                         Model:  2D Square Hubbard" << std::endl;
    os << "                            nx:  " << n << std::endl;
    os << "                            ny:  " << n << std::endl;
    os << "               Number of sites:  " << nsites << std::endl;
    os << "                             t:  " << t[0] << std::endl;
    os << "                            tp:  " << t[1] << std::endl;
    os << "                           tpp:  " << t[2] << std::endl;
    os << "                             U:  " << U << std::endl;
    os << "                            mu:  " << mu << std::endl;
    os << "                          dtau:  " << dtau << std::endl;
    os << "         Number of time slices:  " << M << std::endl;
    os << "                          beta:  " << M * dtau << std::endl;

    // restore saved format flags
    os.flags(oldFlags);
}

void SquareHubbard::Trans_Avrg(double *datamat, double *datavec) const {
    for (int d = 0; d < nsites; d++) {
        int dx = d / n, dy = d % n;
        datavec[d] = 0.;
        for (int i = 0; i < nsites; i++) {
            int x = i / n, y = i % n;
            int target = (x + dx) % n * n + (y + dy) % n;
            datavec[d] += datamat[i * nsites + target];
        }
        datavec[d] /= nsites;
    }
}

void SquareHubbard::FT(double *data, std::complex<double> *ftdata) const {
    for (int k = 0; k < nsites; k++) {
        int kx = k / n, ky = k % n;
        ftdata[k] = 0.;
        for (int i = 0; i < nsites; i++) {
            int x = i / n, y = i % n;
            ftdata[k] += data[i] * FT_Factor[(kx * x + ky * y) % n];
        }
        // ftdata[k] /= double(nsites);
    }
}

void SquareHubbard::Print_ETM(const std::string &name, double *data, double *data_err,
                              std::ostream &os) const {
    // save current format flags
    std::ios::fmtflags oldFlags = os.flags();

    os.setf(std::ios::scientific);

    os << name << ":" << std::endl;
    for (int i = 0; i < nsites; i++) {
        int x = i / n, y = i % n;
        os << std::setw(3) << x << "  " << std::setw(3) << y << "    " << std::setw(13)
           << data[i] << " +- " << std::setw(12) << data_err[i] << std::endl;
    }

    // restore saved format flags
    os.flags(oldFlags);
}

void SquareHubbard::Print_TDM(const std::string &name, double *data, double *data_err,
                              std::ostream &os) const {
    // save current format flags
    std::ios::fmtflags oldFlags = os.flags();

    for (int i = 0; i < nsites; i++) {
        int x = i / n, y = i % n;
        os << name << " x = " << x << "  "
           << "y = " << y << ":" << std::endl;
        for (int l = 0; l < M; l++) {
            os << std::setiosflags(std::ios::fixed) << std::setw(10) << l * dtau
               << "    " << std::resetiosflags(std::ios::floatfield)
               << std::setiosflags(std::ios::scientific) << std::setw(13)
               << data[l * nsites + i] << " +- " << std::setw(12)
               << data_err[l * nsites + i] << std::resetiosflags(std::ios::floatfield)
               << std::endl;
        }
        os << std::endl;
    }

    // restore saved format flags
    os.flags(oldFlags);
}

void SquareHubbard::Print_ETM_FT(const std::string &name, std::complex<double> *data,
                                 std::complex<double> *data_err,
                                 std::ostream &os) const {
    // save current format flags
    std::ios::fmtflags oldFlags = os.flags();

    os.setf(std::ios::scientific);

    os << name << ":" << std::endl;
    for (int k = 0; k < nsites; k++) {
        int kx = k / n, ky = k % n;
        os << std::setw(3) << kx << "  " << std::setw(3) << ky << "    ("
           << std::setw(13) << data[k].real() << " +- " << std::setw(12)
           << data_err[k].real() << ") + I(" << std::setw(13) << data[k].imag()
           << " +- " << std::setw(12) << data_err[k].imag() << ")" << std::endl;
    }

    // restore saved format flags
    os.flags(oldFlags);
}

void SquareHubbard::Print_TDM_FT(const std::string &name, std::complex<double> *data,
                                 std::complex<double> *data_err,
                                 std::ostream &os) const {
    // save current format flags
    std::ios::fmtflags oldFlags = os.flags();

    for (int k = 0; k < nsites; k++) {
        int kx = k / n, ky = k % n;
        os << name << " kx = " << kx << "  "
           << "ky = " << ky << ":" << std::endl;
        for (int l = 0; l < M; l++) {
            os << std::setiosflags(std::ios::fixed) << std::setw(10) << l * dtau
               << "    (" << std::resetiosflags(std::ios::floatfield)
               << std::setiosflags(std::ios::scientific) << std::setw(13)
               << data[l * nsites + k].real() << " +- " << std::setw(12)
               << data_err[l * nsites + k].real() << ") + I(" << std::setw(13)
               << data[l * nsites + k].imag() << " +- " << std::setw(12)
               << data_err[l * nsites + k].imag() << ")"
               << std::resetiosflags(std::ios::floatfield) << std::endl;
        }
        os << std::endl;
    }

    // restore saved format flags
    os.flags(oldFlags);
}