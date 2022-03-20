#include "gfun.h"

#include <cstring>

#include "mkl.h"
#include "util.h"

Gfun::Gfun(const Model* model, double*& emV_, int nWrap_, int npack_, double toldiff_)
    : n(model->nsites),
      M(model->M),
      emV(emV_),
      l(-2),
      nWrap(nWrap_),
      wps(nWrap_),
      toldiff(toldiff_),
      mb(model),
      sb(model->nsites, model->M, npack_),
      ndelay(0) {
    G = new double[n * n];
    V = new double[n * n];
    W = new double[n * n];
    D2 = new double[n];
    pvt1 = new int[n];
    pvt2 = new int[n];
}

Gfun::~Gfun() {
    delete[] pvt2;
    delete[] pvt1;
    delete[] D2;
    delete[] W;
    delete[] V;
    delete[] G;
}

void Gfun::GetG(int lr) {
    bool match = (lr == (l + 1) % M);
    if (nWrap > 0 && match) {
        mb.MultB_Left(&emV[l * n], G);
        mb.MultBi_Right(&emV[l * n], G);
        wps--;

        if (wps == 0) {
            memcpy(V, G, n * n * sizeof(double));

            ComputeG(lr);

            double diff = MatDiff(n, V, G);

            if (diff > toldiff) nWrap--;
            wps = nWrap;
        }
    }
    else {
        ComputeG(lr);
        wps = nWrap;
    }

    l = lr;
    //printmat("gf", n, n, n, G);
}

void Gfun::ComputeG(int lr) {
    logdet = 0.;
    sign = 1;
    sb.Get_SeqMultB((lr - 1 + M) % M, lr, mb, emV);

    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, sb.T, n, pvt1);
    logdet += LogDet(n, sb.T);
    sign *= Sgn(n, sb.T, pvt1);

    memcpy(W, sb.U, n * n * sizeof(double));

    LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'T', n, n, sb.T, n, pvt1, W, n);

    // sb.D --> inv(sb.D) * D2
    for (int i = 0; i < n; i++) {
        if (std::abs(sb.D[i]) > 1.) {
            D2[i] = sb.D[i] > 0 ? 1 : -1;
            sb.D[i] = 1. / abs(sb.D[i]);
        }
        else {
            D2[i] = sb.D[i];
            sb.D[i] = 1.;
        }
    }

    ScaleCols(n, W, sb.D);
    for (int i = 0; i < n; i++) W[i * n + i] += D2[i];

    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, W, n, pvt2);
    logdet += LogDet(n, W);
    sign *= Sgn(n, W, pvt2);

    Trans(n, G, sb.U);
    ScaleRows(n, G, sb.D);

    LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'T', n, n, W, n, pvt2, G, n);

    LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', n, n, sb.T, n, pvt1, G, n);

    LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, sb.U, n, pvt2);
    sign *= Sgn(n, sb.U, pvt2);

    logdet = -logdet;
    for (int i = 0; i < n; i++) {
        logdet += std::log(std::abs(sb.D[i]));
    }
}

void Gfun::UpdateG(int j, double gamma, double r) {
    double* x = &W[ndelay * n];
    double* y = &V[ndelay * n];

    cblas_dcopy(n, &G[j], n, x, 1);
    x[j] -= 1.;
    cblas_dcopy(n, &G[j * n], 1, y, 1);

    if (ndelay > 0) {
        cblas_dgemv(CblasRowMajor, CblasTrans, ndelay, n, 1., W, n, &V[j], n, 1., x, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, ndelay, n, 1., V, n, &W[j], n, 1., y, 1);
    }
    cblas_dscal(n, gamma / r, x, 1);
    logdet -= std::log(std::abs(r));
    sign = r > 0 ? sign : -sign;
    ndelay++;
}

void Gfun::ApplyUpdate() {
    if (ndelay > 0)
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, ndelay, 1., W, n, V,
                    n, 1., G, n);
    ndelay = 0;
}

double Gfun::GetGjj(int j) {
    double gjj = G[j * n + j];
    if (ndelay > 0) gjj += cblas_ddot(ndelay, &W[j], n, &V[j], n);
    return gjj;
}