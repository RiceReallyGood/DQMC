#include "dqmc.h"

#include "mkl.h"
#include "util.h"

DQMC::DQMC(const Model* model_, int nWarm_, int nBins_, int BinSize_, int MeasFreq_,
           int nWrap_, int npack_, double toldiff_, std::ostream& os_)
    : model(model_),
      M(model_->M),
      n(model_->nsites),
      nWarm(nWarm_),
      nBins(nBins_),
      BinSize(BinSize_),
      MeasFreq(MeasFreq_),
      tdm(model_->nsites, model_->M, npack_),
      gt_up(model_, emV_up, npack_),
      gt_dn(model_, emV_dn, npack_),
      gf_up(model_, emV_up, nWrap_, npack_, toldiff_),
      gf_dn(model_, emV_dn, npack_, npack_, toldiff_),
      e(time(nullptr)),
      uid(0, model_->M - 1),
      urd(0, 1),
      os(os_) {
    std::uniform_int_distribution<int> zo(0, 1);
    HSF = new int[M * n];
    for (int i = 0; i < M * n; i++) HSF[i] = 2 * zo(e) - 1;
    // for (int i = 0; i < M * n; i++)
    //     std::cout << "HSF[" << i / n << ", " << i % n << "] = " << HSF[i] <<
    //     std::endl;

    double temp = std::exp(model_->dtau * model_->U / 2.);
    double lambda = std::log(temp + std::sqrt(temp * temp - 1));
    // std::cout << "lambda = " << lambda << std::endl;

    explookup = new double[5];
    explookup = explookup + 2;
    for (int i = -2; i <= 2; i++) {
        explookup[i] = std::exp(i * lambda);
        // std::cout << "explookup[" << i << "] = " << explookup[i] << std::endl;
    }

    emV_up = new double[M * n];
    emV_dn = new double[M * n];

    for (int i = 0; i < M * n; i++) {
        emV_up[i] = explookup[-HSF[i]];
        emV_dn[i] = explookup[HSF[i]];
    }

    // printmat("emV_up", M, n, n, emV_up);
    // printmat("emV_dn", M, n, n, emV_dn);

    giup = new double[(nBins + 1) * M * n];
    giup_err = new double[M * n];
    gkup = new std::complex<double>[(nBins + 1) * M * n];
    gkup_err = new std::complex<double>[M * n];
    gidn = new double[(nBins + 1) * M * n];
    gidn_err = new double[M * n];
    gkdn = new std::complex<double>[(nBins + 1) * M * n];
    gkdn_err = new std::complex<double>[M * n];
    giavrg = new double[(nBins + 1) * M * n];
    giavrg_err = new double[M * n];
    gkavrg = new std::complex<double>[(nBins + 1) * M * n];
    gkavrg_err = new std::complex<double>[M * n];
    SzSzi = new double[(nBins + 1) * M * n];
    SzSzi_err = new double[M * n];
    SzSzk = new std::complex<double>[(nBins + 1) * M * n];
    SzSzk_err = new std::complex<double>[M * n];
    dendeni = new double[(nBins + 1) * M * n];
    dendeni_err = new double[M * n];
    dendenk = new std::complex<double>[(nBins + 1) * M * n];
    dendenk_err = new std::complex<double>[M * n];
}

void DQMC::Sweep() {
    for (int l = 0; l < M; l++) {
        gf_up.GetG(l);
        gf_dn.GetG(l);
        for (int i = 0; i < n; i++) {
            int& hsf = HSF[l * n + i];
            double gammaup = explookup[2 * hsf] - 1.;
            double gammadn = explookup[-2 * hsf] - 1.;
            double rup = 1. + (1. - gf_up.GetGjj(i)) * gammaup;
            double rdn = 1. + (1. - gf_dn.GetGjj(i)) * gammadn;

            double r = rup * rdn;

            if (urd(e) < std::abs(r)) {
                hsf = -hsf;
                std::swap(emV_up[l * n + i], emV_dn[l * n + i]);
                gf_up.UpdateG(i, gammaup, rup);
                gf_dn.UpdateG(i, gammadn, rdn);
            }
        }
        gf_up.ApplyUpdate();
        gf_dn.ApplyUpdate();
        // printmat("G updated", n, n, n, gf_up.G);
        // gf_up.ComputeG(l);
        // printmat("G computed", n, n, n, gf_up.G);
    }
}

void DQMC::Measure() {
    gt_up.Get_FullG();
    gt_dn.Get_FullG();
    tdm.measure(gt_up, gt_dn);
}

void DQMC::GetTrans(int iBin) {
    double* datamat;
    double* datavec;
    for (int l = 0; l < M; l++) {
        // gup
        datamat = &(tdm.gtup[l * n * n]);
        datavec = &giup[iBin * M * n + l * n];
        model->Trans_Avrg(datamat, datavec);

        // gdn
        datamat = &(tdm.gtdn[l * n * n]);
        datavec = &gidn[iBin * M * n + l * n];
        model->Trans_Avrg(datamat, datavec);

        // SzSz
        datamat = &(tdm.SzSz[l * n * n]);
        datavec = &SzSzi[iBin * M * n + l * n];
        model->Trans_Avrg(datamat, datavec);

        // denden
        datamat = &(tdm.denden[l * n * n]);
        datavec = &dendeni[iBin * M * n + l * n];
        model->Trans_Avrg(datamat, datavec);
    }
}

void DQMC::GetFT(int iBin) {
    double* datar;
    std::complex<double>* datac;
    for (int l = 0; l < M; l++) {
        // gup
        datar = &giup[iBin * M * n + l * n];
        datac = &gkup[iBin * M * n + l * n];
        model->FT(datar, datac);

        // gdn
        datar = &gidn[iBin * M * n + l * n];
        datac = &gkdn[iBin * M * n + l * n];
        model->FT(datar, datac);

        // Sdenden
        datar = &dendeni[iBin * M * n + l * n];
        datac = &dendenk[iBin * M * n + l * n];
        model->FT(datar, datac);

        // SzSz
        datar = &SzSzi[iBin * M * n + l * n];
        datac = &SzSzk[iBin * M * n + l * n];
        model->FT(datar, datac);
    }
}

void DQMC::Get_Statistics() {
    sign /= double(nBins);
    vdAdd(nBins * M * n, &giup[M * n], &gidn[M * n], &giavrg[M * n]);
    cblas_dscal(nBins * M * n, 0.5, &giavrg[M * n], 1);
    vdAdd(2 * nBins * M * n, reinterpret_cast<double*>(&gkup[M * n]),
          reinterpret_cast<double*>(&gkdn[M * n]),
          reinterpret_cast<double*>(&gkavrg[M * n]));
    cblas_dscal(2 * nBins * M * n, 0.5, reinterpret_cast<double*>(&gkavrg[M * n]), 1);

    Get_Avrg_Err(nBins, M * n, &giavrg[M * n], giavrg, giavrg_err);
    Get_Avrg_Err(nBins, M * n, &giup[M * n], giup, giup_err);
    Get_Avrg_Err(nBins, M * n, &gidn[M * n], gidn, gidn_err);

    Get_Avrg_Err(nBins, 2 * M * n, reinterpret_cast<double*>(&gkavrg[M * n]),
                 reinterpret_cast<double*>(gkavrg),
                 reinterpret_cast<double*>(gkavrg_err));
    Get_Avrg_Err(nBins, 2 * M * n, reinterpret_cast<double*>(&gkup[M * n]),
                 reinterpret_cast<double*>(gkup), reinterpret_cast<double*>(gkup_err));
    Get_Avrg_Err(nBins, 2 * M * n, reinterpret_cast<double*>(&gkdn[M * n]),
                 reinterpret_cast<double*>(gkdn), reinterpret_cast<double*>(gkdn_err));

    Get_Avrg_Err(nBins, M * n, &SzSzi[M * n], SzSzi, SzSzi_err);
    Get_Avrg_Err(nBins, 2 * M * n, reinterpret_cast<double*>(&SzSzk[M * n]),
                 reinterpret_cast<double*>(SzSzk),
                 reinterpret_cast<double*>(SzSzk_err));

    Get_Avrg_Err(nBins, M * n, &dendeni[M * n], dendeni, dendeni_err);
    Get_Avrg_Err(nBins, 2 * M * n, reinterpret_cast<double*>(&dendenk[M * n]),
                 reinterpret_cast<double*>(dendenk),
                 reinterpret_cast<double*>(dendenk_err));
}

void DQMC::PrintProfile() {
    model->Print_Parms(os);

    // save current format flags
    std::ios::fmtflags oldFlags = os.flags();
    os.setf(std::ios::fixed);

    os << "       Number of warmup sweeps:  " << nWarm << std::endl;
    os << "                Number of Bins:  " << nBins << std::endl;
    os << "              Size of each Bin:  " << BinSize << std::endl;
    os << "      Frequency of measurement:  " << MeasFreq << std::endl;
    os << std::string(80, '=') << std::endl;

    // restore saved format flags
    os.flags(oldFlags);
}

void DQMC::Simulate() {
    for (int i = 0; i < nWarm; i++) {
        Sweep();
        std::cout << "Warm sweep: " << i << std::endl;
    }

    for (int iBin = 1; iBin <= nBins; iBin++) {
        tdm.clear();

        for (int nMeas = 0; nMeas < BinSize; nMeas++) {
            for (int pass = 0; pass < MeasFreq; pass++) {
                Sweep();
            }

            std::cout << "iBin, nMeas : " << iBin << " " << nMeas << std::endl;
            Measure();
        }

        tdm.finalize();

        sign += tdm.sign;
        GetTrans(iBin);
        GetFT(iBin);
    }
    Get_Statistics();
}

void DQMC::PrintResult() {
    // save current format flags
    std::ios::fmtflags oldFlags = os.flags();
    os.setf(std::ios::fixed);
    os << "                  Average Sign:  " << sign << std::endl;
    // restore saved format flags
    os.flags(oldFlags);

    model->Print_TDM("G tau average  ", giavrg, giavrg_err, os);
    model->Print_TDM("G tau up       ", giup, giup_err, os);
    model->Print_TDM("G tau down     ", gidn, gidn_err, os);

    model->Print_TDM("SzSz           ", SzSzi, SzSzi_err, os);
    model->Print_TDM("denden         ", dendeni, dendeni_err, os);

    model->Print_TDM_FT("G tau average k  ", gkavrg, gkavrg_err, os);
    model->Print_TDM_FT("G tau up k  ", gkup, gkup_err, os);
    model->Print_TDM_FT("G tau down k", gkdn, gkdn_err, os);

    model->Print_TDM_FT("SzSz k      ", SzSzk, SzSzk_err, os);
    model->Print_TDM_FT("denden k    ", dendenk, dendenk_err, os);
}

DQMC::~DQMC() {
    delete[] dendenk;
    delete[] dendenk_err;
    delete[] dendeni;
    delete[] dendeni_err;

    delete[] SzSzk;
    delete[] SzSzk_err;
    delete[] SzSzi;
    delete[] SzSzi_err;

    delete[] gkavrg;
    delete[] gkavrg_err;

    delete[] giavrg;
    delete[] giavrg_err;

    delete[] gkdn;
    delete[] gkdn_err;
    delete[] gkup;
    delete[] gkup_err;

    delete[] gidn;
    delete[] gidn_err;
    delete[] giup;
    delete[] giup_err;

    delete[] emV_dn;
    delete[] emV_up;

    explookup = explookup - 2;
    delete[] explookup;

    delete[] HSF;
}