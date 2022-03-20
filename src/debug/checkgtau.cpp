#include <random>
#include <sstream>

#include "../RingHubbard.h"
#include "../gfun.h"
#include "../gtau.h"
#include "../util.h"

int main() {
    int n = 8, M = 10;
    double t = 1., U = 4., mu = 0., dtau = 0.1;
    RingHubbard r(n, t, U, mu, dtau, M);

    int nWrap = 5, npack = 2;
    double *emV;
    Gfun gf(&r, emV, nWrap, 1, 1e-6);
    Gtau gt1(&r, emV, 1);
    Gtau gtpack(&r, emV, npack);

    emV = new double[M * n];
    double temp = std::exp(dtau * U / 2);
    double lambda = std::log(temp + std::sqrt(temp * temp - 1));

    std::default_random_engine dre(time(NULL));
    std::uniform_int_distribution<int> uid(0, 1);
    for (int i = 0; i < M * n; i++) {
        emV[i] = std::exp(-(2 * uid(dre) - 1) * lambda);
    }

    gt1.Get_FullG();

    for (int blk = 0; blk < gt1.nblk; blk++) {
        int l = blk * gt1.npack;
        std::ostringstream s;
        s << "G(" << l << ")";
        printmat("\n" + s.str() + " from Gtau 1", n, n, gt1.dim,
                 &(gt1.G[blk * n * gt1.dim + blk * n]));

        gf.ComputeG(l);
        printmat(s.str() + " from Gfun", n, n, n, gf.G);
    }

    gtpack.Get_FullG();
    for (int blk = 0; blk < gtpack.nblk; blk++) {
        int l = blk * gtpack.npack;
        std::ostringstream s;
        s << "G(" << l << ")";
        printmat("\n" + s.str() + " from Gtau pack", n, n, gtpack.dim,
                 &(gtpack.G[blk * n * gtpack.dim + blk * n]));

        gf.ComputeG(l);
        printmat(s.str() + " from Gfun", n, n, n, gf.G);
    }

    int inc = npack / 2, dec = npack - 1 - inc;
    for (int blk0 = 0; blk0 < gtpack.nblk; blk0++) {
        for (int blkt = 0; blkt < gtpack.nblk; blkt++) {
            gtpack.Load_G_to_g(blk0, blkt);
            const int &lt = gtpack.lt;
            const int &l0 = gtpack.l0;
            std::ostringstream sp5, sp1;
            sp5 << "gtpack(" << lt << "," << l0 << ")";
            sp1 << "gt1(" << lt << "," << l0 << ")";
            printmat(sp5.str(), n, n, n, gtpack.gt0);
            printmat(sp1.str(), n, n, gt1.dim,
                     &gt1.G[lt * n * gt1.dim + l0 * n]);

            for (int i = 0; i < inc; i++) {
                gtpack.Increase_g_time();
                sp5.str("");
                sp1.str("");
                sp5 << "gtpack(" << lt << "," << l0 << ")";
                sp1 << "gt1(" << lt << "," << l0 << ")";
                printmat(sp5.str(), n, n, n, gtpack.gt0);
                printmat(sp1.str(), n, n, gt1.dim,
                         &gt1.G[lt * n * gt1.dim + l0 * n]);
            }

            if (dec > 0) {
                gtpack.Load_G_to_g(blk0, blkt);
            }

            for (int i = 0; i < dec; i++) {
                gtpack.Decrease_g_time();
                sp5.str("");
                sp1.str("");
                sp5 << "gtpack5(" << lt << "," << l0 << ")";
                sp1 << "gtpack1(" << lt << "," << l0 << ")";
                printmat(sp5.str(), n, n, n, gtpack.gt0);
                printmat(sp1.str(), n, n, gt1.dim,
                         &gt1.G[lt * n * gt1.dim + l0 * n]);
            }
        }
    }

    // gt.Load_G_to_g(0, 0);

    // for (int i = 1; i < 2 * M; i++) {
    //     int lt = i % M;
    //     std::ostringstream st0, s0t, stt;
    //     st0 << "Gt0(" << lt << ", 0)";
    //     s0t << "G0t(0, " << lt << ")";
    //     stt << "Gtt(" << lt << ", " << lt << ")";
    //     gt.Increase_g_time();
    //     printmat("\n" + st0.str() + " from Increase", n, n, n, gt.gt0);
    //     printmat("\n" + st0.str() + " from Full Gtau", n, n, gt.dim,
    //              &(gt.G[lt * n * gt.dim]));

    //     printmat("\n" + s0t.str() + " from Increase", n, n, n, gt.g0t);
    //     printmat("\n" + s0t.str() + " from Full Gtau", n, n, gt.dim, &(gt.G[lt *
    //     n]));

    //     printmat("\n" + stt.str() + " from Increase", n, n, n, gt.gtt);
    //     printmat("\n" + stt.str() + " from Full Gtau", n, n, gt.dim,
    //              &(gt.G[lt * n * gt.dim + lt * n]));
    // }

    // gt.Load_G_to_g(0, 0);

    // for (int i = 1; i < 2 * M; i++) {
    //     int lt = (2 * M - i) % M;
    //     std::ostringstream st0, s0t, stt;
    //     st0 << "Gt0(" << lt << ", 0)";
    //     s0t << "G0t(0, " << lt << ")";
    //     stt << "Gtt(" << lt << ", " << lt << ")";
    //     gt.Decrease_g_time();
    //     printmat("\n" + st0.str() + " from Decrease", n, n, n, gt.gt0);
    //     printmat("\n" + st0.str() + " from Full Gtau", n, n, gt.dim,
    //              &(gt.G[lt * n * gt.dim]));

    //     printmat("\n" + s0t.str() + " from Decrease", n, n, n, gt.g0t);
    //     printmat("\n" + s0t.str() + " from Full Gtau", n, n, gt.dim, &(gt.G[lt *
    //     n]));

    //     printmat("\n" + stt.str() + " from Decrease", n, n, n, gt.gtt);
    //     printmat("\n" + stt.str() + " from Full Gtau", n, n, gt.dim,
    //              &(gt.G[lt * n * gt.dim + lt * n]));
    // }

    delete[] emV;
}