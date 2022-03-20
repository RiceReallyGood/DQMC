#include <iostream>

#include "../cfg.h"
using namespace std;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Please input configure file name" << endl;
        cout << "Usage: ./checkcfg filename" << endl;
        exit(1);
    }
    Config cfg(argv[1]);

    cout << "ofile = " << cfg.GetCfg("ofile") << endl;
    cout << "gfile = " << cfg.GetCfg("gfile") << endl;
    cout << "L     = " << cfg.GetCfg("L") << endl;
    cout << "nbin  = " << cfg.GetCfg("nbin") << endl;
    cout << "hhhh  = " << cfg.GetCfg("hhhh") << endl;
    cout << "tdm   = " << cfg.GetCfg("tdm") << endl;
    cout << "nx    = " << cfg.GetCfg("nx") << endl;
}