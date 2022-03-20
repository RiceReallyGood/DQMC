#ifndef DQMC_GTAU_H
#define DQMC_GTAU_H

#include "matb.h"
#include "model.h"

class Gtau {
public:
    Gtau(const Model* model, double*& emV_, int npack_);
    Gtau(const Gtau&) = delete;
    Gtau& operator=(const Gtau&) = delete;
    ~Gtau();

    int npack;
    int nblk;  // Number of blocks
    int dim;   // Dimension of matrix G
    double* G;
    double* g00;  //[g00]_{ij} = -<c^{\dagger}_{j}(l0) c_{i}(l0)>
    double* gt0;  //[gt0]_{ij} = <c_{i}(lt) c^{\dagger}_{j}(l0) > for lt >= l0
                  //           = -<c^{\dagger}_{j}(l0) c_{i}(lt)> for lt <  l0
    double* g0t;  //[g0t]_{ij} = <c_{i}(l0) c^{\dagger}_{j}(lt) > for lt <  l0
                  //           = -<c^{\dagger}_{j}(lt) c_{i}(l0)> for lt >= l0
    double* gtt;  //[gtt]_{ij} = -<c^{\dagger}_{j}(lt) c_{i}(lt)>
    int l0;
    int lt;

    int sign;  // Sign of Det(G)

    void Get_FullG();
    void Load_G_to_g(int blk0, int blkt);
    void Increase_g_time();
    void Decrease_g_time();

private:
    int n;          // Number of sites
    int M;          // Number of time slices
    double*(&emV);  // Matrix info of emV(0)...emV(M - 1)

    MatB mb;

    // Work Space
    double* WS;
    int* pvt;
};

#endif  // DQMC_GTAU_H