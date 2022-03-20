#include <fstream>

#include "RingHubbard.h"
#include "SquareHubbard.h"
#include "cfg.h"
#include "dqmc.h"
#include "model.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Please input configure file name!" << std::endl;
        std::cout << "Usage: ./dqmc filename" << std::endl;
        exit(1);
    }

    Config cfg(argv[1]);

    std::string mname = cfg.GetCfg("model");
    if (mname != "1D Ring Hubbard" && mname != "2D Square Hubbard") {
        std::cout << "Please specify which model to simulate!\n"
                  << "Add "
                  << "\"model = 1D Ring Hubbard\""
                  << " or "
                  << "\"model = 2D Square Hubbard\""
                  << " in input file." << std::endl;
        exit(1);
    }

    std::string val;
    val = cfg.GetCfg("n");
    int n;
    if (val == "") {
        std::cout << "n not specified, use default n = 4" << std::endl;
        n = 4;
    }
    else {
        n = std::stoi(val);
    }

    val = cfg.GetCfg("t");
    double t;
    if (val == "") {
        std::cout << "t not specified, use default t = 1" << std::endl;
        t = 1.;
    }
    else {
        t = std::stod(val);
    }

    val = cfg.GetCfg("U");
    double U;
    if (val == "") {
        std::cout << "U not specified, use default U = 4" << std::endl;
        U = 4.;
    }
    else {
        U = std::stod(val);
    }

    val = cfg.GetCfg("mu");
    double mu;
    if (val == "") {
        std::cout << "mu not specified, use default mu = 0" << std::endl;
        mu = 0.;
    }
    else {
        mu = std::stod(val);
    }

    val = cfg.GetCfg("dtau");
    double dtau;
    if (val == "") {
        std::cout << "dtau not specified, use default mdau = 0.1" << std::endl;
        dtau = 0.1;
    }
    else {
        dtau = std::stod(val);
    }

    val = cfg.GetCfg("M");
    int M;
    if (val == "") {
        std::cout << "M not specified, use default mdau = 10" << std::endl;
        M = 10;
    }
    else {
        M = std::stod(val);
    }

    Model* model;
    if (mname == "1D Ring Hubbard") {
        model = new RingHubbard(n, t, U, mu, dtau, M);
    }
    else {
        double tp = std::stod(cfg.GetCfg("tp"));
        double tpp = std::stod(cfg.GetCfg("tpp"));
        model = new SquareHubbard(n, t, tp, tpp, U, mu, dtau, M);
    }

    val = cfg.GetCfg("nWarm");
    int nWarm;
    if (val == "") {
        std::cout << "nWarm not specified, use default nWarm = 10000" << std::endl;
        nWarm = 10000;
    }
    else {
        nWarm = std::stoi(val);
    }

    val = cfg.GetCfg("nBins");
    int nBins;
    if (val == "") {
        std::cout << "nBins not specified, use default nBins = 32" << std::endl;
        nBins = 32;
    }
    else {
        nBins = std::stoi(val);
    }

    val = cfg.GetCfg("BinSize");
    int BinSize;
    if (val == "") {
        std::cout << "BinSize not specified, use default BinSize = 10000" << std::endl;
        BinSize = 10000;
    }
    else {
        BinSize = std::stoi(val);
    }

    val = cfg.GetCfg("MeasFreq");
    int MeasFreq;
    if (val == "") {
        std::cout << "MeasFreq not specified, use default MeasFreq = 5" << std::endl;
        MeasFreq = 5;
    }
    else {
        MeasFreq = std::stoi(val);
    }

    val = cfg.GetCfg("nWrap");
    int nWrap;
    if (val == "") {
        std::cout << "nWrap not specified, use default nWrap = 10" << std::endl;
        nWrap = 10;
    }
    else {
        nWrap = std::stoi(val);
    }

    val = cfg.GetCfg("npack");
    int npack;
    if (val == "") {
        std::cout << "npack not specified, use default npack = 5" << std::endl;
        npack = 5;
    }
    else {
        npack = std::stoi(val);
    }

    val = cfg.GetCfg("toldiff");
    double toldiff;
    if (val == "") {
        std::cout << "toldiff not specified, use default toldiff = 1e-6" << std::endl;
        toldiff = 1e-6;
    }
    else {
        toldiff = std::stod(val);
    }

    val = cfg.GetCfg("ofile");
    std::ostream* osp;
    if (val == "") {
        std::cout << "ofile not specified, use default out put stream cout"
                  << std::endl;
        osp = &std::cout;
    }
    else {
        osp = new std::ofstream(val);
        if (!(*osp)) {
            std::cout << "Can not open file " << val << std::endl;
            exit(1);
        }
    }

    DQMC mysystem(model, nWarm, nBins, BinSize, MeasFreq, nWrap, npack, toldiff, *osp);
    mysystem.PrintProfile();
    mysystem.Simulate();
    mysystem.PrintResult();

    delete model;
}
