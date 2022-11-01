#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"
#include "TROOT.h"
#include "TChain.h"
#include "TNtuple.h"

using namespace std;

Bool_t reject;
Double_t fline(Double_t *x, Double_t *par)
{
    if (reject && x[0] > 1.880 && x[0] < 1.898) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0];
}

Double_t fquadratic(Double_t *x, Double_t *par)
{
    if (reject && x[0] > 1.880 && x[0] < 1.898) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t fcubic(Double_t *x, Double_t *par)
{
    if (reject && x[0] > 1.880 && x[0] < 1.898) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}



void Danalyza_pp()
{

    //   TFile* data = new TFile("2021-06-13_13-41_D0_4536.picoD0AnaMaker.root");
    TChain *ntp = new TChain("ntp_signal");
    ntp->Add("1merged_output_3_000.root");
    ntp->Add("1merged_output_3_001.root");
    ntp->Add("1merged_output_3_002.root");
    ntp->Add("1merged_output_3_003.root");
    ntp->Add("1merged_output_3_004.root");
    ntp->Add("1merged_output_3_005.root");
    ntp->Add("1merged_output_3_006.root");

    ntp->Add("2merged_output_3_000.root");
    ntp->Add("2merged_output_3_001.root");
    ntp->Add("2merged_output_3_002.root");
    ntp->Add("2merged_output_3_003.root");
    ntp->Add("2merged_output_3_004.root");
    ntp->Add("2merged_output_3_005.root");
    ntp->Add("2merged_output_3_006.root");

    ntp->Add("3merged_output_3_000.root");
    ntp->Add("3merged_output_3_001.root");
    ntp->Add("3merged_output_3_002.root");
    ntp->Add("3merged_output_3_003.root");
    ntp->Add("3merged_output_3_004.root");

//    ntp->Add("2022-10-11_05-19_D0_31.picoD0AnaMaker.root");


    TChain *bntp = new TChain("ntp_background");
    bntp->Add("1merged_output_3_000.root");
    bntp->Add("1merged_output_3_001.root");
    bntp->Add("1merged_output_3_002.root");
    bntp->Add("1merged_output_3_003.root");
    bntp->Add("1merged_output_3_004.root");
    bntp->Add("1merged_output_3_005.root");
    bntp->Add("1merged_output_3_006.root");

    bntp->Add("2merged_output_3_000.root");
    bntp->Add("2merged_output_3_001.root");
    bntp->Add("2merged_output_3_002.root");
    bntp->Add("2merged_output_3_003.root");
    bntp->Add("2merged_output_3_004.root");
    bntp->Add("2merged_output_3_005.root");
    bntp->Add("2merged_output_3_006.root");

    bntp->Add("3merged_output_3_000.root");
    bntp->Add("3merged_output_3_001.root");
    bntp->Add("3merged_output_3_002.root");
    bntp->Add("3merged_output_3_003.root");
    bntp->Add("3merged_output_3_004.root");

//    bntp->Add("2022-10-11_05-19_D0_31.picoD0AnaMaker.root");


    TChain *Dntp = new TChain("Dntp_signal");
    Dntp->Add("1merged_output_3_000.root");
    Dntp->Add("1merged_output_3_001.root");
    Dntp->Add("1merged_output_3_002.root");
    Dntp->Add("1merged_output_3_003.root");
    Dntp->Add("1merged_output_3_004.root");
    Dntp->Add("1merged_output_3_005.root");
    Dntp->Add("1merged_output_3_006.root");

    Dntp->Add("2merged_output_3_000.root");
    Dntp->Add("2merged_output_3_001.root");
    Dntp->Add("2merged_output_3_002.root");
    Dntp->Add("2merged_output_3_003.root");
    Dntp->Add("2merged_output_3_004.root");
    Dntp->Add("2merged_output_3_005.root");
    Dntp->Add("2merged_output_3_006.root");

    Dntp->Add("3merged_output_3_000.root");
    Dntp->Add("3merged_output_3_001.root");
    Dntp->Add("3merged_output_3_002.root");
    Dntp->Add("3merged_output_3_003.root");
    Dntp->Add("3merged_output_3_004.root");

//    Dntp->Add("2022-10-11_05-19_D0_31.picoD0AnaMaker.root");


    TChain *Dbntp = new TChain("Dntp_background");
    Dbntp->Add("1merged_output_3_000.root");
    Dbntp->Add("1merged_output_3_001.root");
    Dbntp->Add("1merged_output_3_002.root");
    Dbntp->Add("1merged_output_3_003.root");
    Dbntp->Add("1merged_output_3_004.root");
    Dbntp->Add("1merged_output_3_005.root");
    Dbntp->Add("1merged_output_3_006.root");

    Dbntp->Add("2merged_output_3_000.root");
    Dbntp->Add("2merged_output_3_001.root");
    Dbntp->Add("2merged_output_3_002.root");
    Dbntp->Add("2merged_output_3_003.root");
    Dbntp->Add("2merged_output_3_004.root");
    Dbntp->Add("2merged_output_3_005.root");
    Dbntp->Add("2merged_output_3_006.root");

    Dbntp->Add("3merged_output_3_000.root");
    Dbntp->Add("3merged_output_3_001.root");
    Dbntp->Add("3merged_output_3_002.root");
    Dbntp->Add("3merged_output_3_003.root");
    Dbntp->Add("3merged_output_3_004.root");

//    Dbntp->Add("2022-10-11_05-19_D0_31.picoD0AnaMaker.root");



    //  TNtuple* ntp = (TNtuple*)myChain -> Get("ntp_signal");
    //  TList* list = (TList*) myChain -> Get("picoD0AnaMaker;1");
    Float_t D0_theta, D0_mass, D0_pt, D0_decayL, k_pt, pi1_pt, pi1_dca, k_dca, k_nSigma, pi1_nSigma, pi1_TOFinvbeta, k_TOFinvbeta, dcaMax, pi1_eventId, k_eventId, dcaDaughters, D_cosThetaStar, dcaD0ToPv, primVz, primVzVpd, k_nHitFit, pi1_nHitFit, pi1_p, k_p;
    ntp -> SetBranchAddress("D_mass", &D0_mass);
    ntp -> SetBranchAddress("D_decayL", &D0_decayL);
    ntp -> SetBranchAddress("D_theta", &D0_theta);
    ntp -> SetBranchAddress("D_pt", &D0_pt);
    ntp -> SetBranchAddress("pi1_pt", &pi1_pt);
    ntp -> SetBranchAddress("k_pt", &k_pt);
    ntp -> SetBranchAddress("pi1_dca", &pi1_dca);
    ntp -> SetBranchAddress("k_dca", &k_dca);
    ntp -> SetBranchAddress("k_nSigma", &k_nSigma);
    ntp -> SetBranchAddress("pi1_nSigma", &pi1_nSigma);
    ntp -> SetBranchAddress("pi1_TOFinvbeta", &pi1_TOFinvbeta);
    ntp -> SetBranchAddress("k_TOFinvbeta", &k_TOFinvbeta);
    ntp -> SetBranchAddress("D_cosThetaStar", &D_cosThetaStar);
    ntp -> SetBranchAddress("dcaD0ToPv", &dcaD0ToPv);
    ntp -> SetBranchAddress("dcaDaughters", &dcaDaughters);
    ntp -> SetBranchAddress("primVz", &primVz);
    ntp -> SetBranchAddress("primVzVpd", &primVzVpd);
    ntp -> SetBranchAddress("k_nHitFit", &k_nHitFit);
    ntp -> SetBranchAddress("pi1_nHitFit", &pi1_nHitFit);
    ntp -> SetBranchAddress("pi1_p", &pi1_p);
    ntp -> SetBranchAddress("k_p", &k_p);


    Float_t bD0_theta, bD0_mass, bD0_pt, bD0_decayL, bk_pt, bpi1_pt, bpi1_dca, bk_dca, bk_nSigma, bpi1_nSigma, bpi1_TOFinvbeta, bk_TOFinvbeta, bdcaMax, bpi1_eventId, bk_eventId, bdcaDaughters, bD_cosThetaStar, bdcaD0ToPv, bprimVz, bprimVzVpd, bk_nHitFit, bpi1_nHitFit, bk_p;
    //   TNtuple* bntp = (TNtuple*)myChain -> Get("ntp_background");
    bntp -> SetBranchAddress("D_mass", &bD0_mass);
    bntp -> SetBranchAddress("D_decayL", &bD0_decayL);
    bntp -> SetBranchAddress("D_theta", &bD0_theta);
    bntp -> SetBranchAddress("D_pt", &bD0_pt);
    bntp -> SetBranchAddress("pi1_pt", &bpi1_pt);
    bntp -> SetBranchAddress("k_pt", &bk_pt);
    bntp -> SetBranchAddress("pi1_dca", &bpi1_dca);
    bntp -> SetBranchAddress("k_dca", &bk_dca);
    bntp -> SetBranchAddress("k_nSigma", &bk_nSigma);
    bntp -> SetBranchAddress("pi1_nSigma", &bpi1_nSigma);
    bntp -> SetBranchAddress("pi1_TOFinvbeta", &bpi1_TOFinvbeta);
    bntp -> SetBranchAddress("k_TOFinvbeta", &bk_TOFinvbeta);
    bntp -> SetBranchAddress("dcaDaughters", &bdcaDaughters);
    bntp -> SetBranchAddress("dcaD0ToPv", &bdcaD0ToPv);
    bntp -> SetBranchAddress("D_cosThetaStar", &bD_cosThetaStar);
    bntp -> SetBranchAddress("primVz", &bprimVz);
    bntp -> SetBranchAddress("primVzVpd", &bprimVzVpd);
    bntp -> SetBranchAddress("k_nHitFit", &bk_nHitFit);
    bntp -> SetBranchAddress("pi1_nHitFit", &bpi1_nHitFit);
    bntp -> SetBranchAddress("k_p", &bk_p);


    Float_t Dstar_mass, Dstar_pt, Sk_pt, Spi1_pt, Spi1_p, Sk_p, Spi2_pt, Spi2_p, Dstar_D0;
    Dntp -> SetBranchAddress("Dstar_mass", &Dstar_mass);
    Dntp -> SetBranchAddress("Dstar_pt", &Dstar_pt);
    Dntp -> SetBranchAddress("pi1_pt", &Spi1_pt);
    Dntp -> SetBranchAddress("k_pt", &Sk_pt);
    Dntp -> SetBranchAddress("pi1_p", &Spi1_p);
    Dntp -> SetBranchAddress("k_p", &Sk_p);
    Dntp -> SetBranchAddress("pi2_pt", &Spi2_pt);
    Dntp -> SetBranchAddress("pi2_p", &Spi2_p);
    Dntp -> SetBranchAddress("triplet_pair", &Dstar_D0);


    Float_t bDstar_mass, bDstar_pt, bSk_pt, bSpi1_pt, bSpi1_p, bSk_p, bSpi2_pt, bSpi2_p, bDstar_D0;
    Dbntp -> SetBranchAddress("Dstar_mass", &bDstar_mass);
    Dbntp -> SetBranchAddress("Dstar_pt", &bDstar_pt);
    Dbntp -> SetBranchAddress("pi1_pt", &bSpi1_pt);
    Dbntp -> SetBranchAddress("k_pt", &bSk_pt);
    Dbntp -> SetBranchAddress("pi1_p", &bSpi1_p);
    Dbntp -> SetBranchAddress("k_p", &bSk_p);
    Dbntp -> SetBranchAddress("pi2_pt", &bSpi2_pt);
    Dbntp -> SetBranchAddress("pi2_p", &bSpi2_p);
    Dbntp -> SetBranchAddress("triplet_pair", &bDstar_D0);



    //  TH1D *hEventStat1 = (TH1D*)list->FindObject("hEventStat1");


    //vzor TH1* h1 = new TH1I("h1", "h1 title", 100-počet binů, 0.0, 4.0 -rozsah);
    TH2F *costhetastar  = new TH2F("costhetastar","costhetastar",100,0,7,100,-1,1);
    TH2F *bcosthetastar  = new TH2F("bcosthetastar","bcosthetastar",100,0,7,100,-1,1);
    TH2F *kaonnsigma  = new TH2F("kaonnsigma","kaonnsigma",500,0,3.5,100,-10,10);
    TH2F *pionnsigma  = new TH2F("pionnsigma","pionnsigma",500,0,3.5,100,-10,10);

    TH2F *kaonnsigma50  = new TH2F("kaonnsigma50","kaonnsigma50",50,0,2,100,-10,10);
  //  TH2F *pionnsigma50  = new TH2F("pionnsigma50","pionnsigma50",50,0,2,100,-10,10);


    TH1D* hInvMassSignD0 = new TH1D("D0 signal", "D0 signal", 32, 1.7, 2.05);
    TH1D* hInvMassBackD0 = new TH1D("D0 background", "D0 background", 32, 1.7, 2.05);

    TH1D* hInvMassSignD0w = new TH1D("D0 signalw", "D0 signalw", 200, 0.5, 2.5);
    TH1D* hInvMassBackD0w = new TH1D("D0 backgroundw", "D0 backgroundw", 200, 0.5, 2.5);
   // TH1D* hStat = (TH1D*) list -> FindObject("hEventStat1");

    TH1D* hInvMassSignD0test = new TH1D("D0 signaltest", "D0 signaltest", 100, 1.4, 2.2);


    TH1D* hInvMassSignD012 = new TH1D("D0 S12", "D0 S12", 32, 1.7, 2.05);
    TH1D* hInvMassSignD023 = new TH1D("D0 S23", "D0 S23", 32, 1.7, 2.05);
    TH1D* hInvMassSignD034 = new TH1D("D0 S34", "D0 S34", 32, 1.7, 2.05);
    TH1D* hInvMassSignD045 = new TH1D("D0 S45", "D0 S45", 32, 1.7, 2.05);

    TH1D* hInvMassSignD012test = new TH1D("D0 S12t", "D0 S12t", 100, 1.4, 2.2);
    TH1D* hInvMassSignD023test = new TH1D("D0 S23t", "D0 S23t", 100, 1.4, 2.2);
    TH1D* hInvMassSignD034test = new TH1D("D0 S34t", "D0 S34t", 100, 1.4, 2.2);
    TH1D* hInvMassSignD045test = new TH1D("D0 S45t", "D0 S45t", 100, 1.4, 2.2);

    TH1D* hInvMassBackD012 = new TH1D("D0 B12", "D0 B12", 32, 1.7, 2.05);
    TH1D* hInvMassBackD023 = new TH1D("D0 B23", "D0 B23", 32, 1.7, 2.05);
    TH1D* hInvMassBackD034 = new TH1D("D0 B34", "D0 B34", 32, 1.7, 2.05);
    TH1D* hInvMassBackD045 = new TH1D("D0 B45", "D0 B45", 32, 1.7, 2.05);

    TH1D* hInvMassSignD012w = new TH1D("D0 S12w", "D0 S12w", 200, 0.5, 2.5);
    TH1D* hInvMassSignD023w = new TH1D("D0 S23w", "D0 S23w", 200, 0.5, 2.5);
    TH1D* hInvMassSignD034w = new TH1D("D0 S34w", "D0 S34w", 200, 0.5, 2.5);
    TH1D* hInvMassSignD045w = new TH1D("D0 S45w", "D0 S45w", 200, 0.5, 2.5);

    TH1D* hInvMassBackD012w = new TH1D("D0 B12w", "D0 B12w", 200, 0.5, 2.5);
    TH1D* hInvMassBackD023w = new TH1D("D0 B23w", "D0 B23w", 200, 0.5, 2.5);
    TH1D* hInvMassBackD034w = new TH1D("D0 B34w", "D0 B34w", 200, 0.5, 2.5);
    TH1D* hInvMassBackD045w = new TH1D("D0 B45w", "D0 B45w", 200, 0.5, 2.5);




    TH1D* hk_pt = new TH1D("hk_pt", "hk_pt", 100, 0, 4);
    TH1D* hpi1_pt = new TH1D("hpi1_pt", "hpi1_pt", 100, 0, 4);
    TH1D* hD0_pt = new TH1D("hD0_pt", "hD0_pt", 100, 0, 7);
    TH1D* hdecayLength = new TH1D("hdecayLength", "hdecayLength", 100, 0, 3);
    TH1D* hpi1_dca = new TH1D("hpi1_dca", "hpi1_dca", 100, 0, 1.5);
    TH1D* hk_dca = new TH1D("hk_dca", "hk_dca", 100, 0, 1.5);
    TH1D* hdcaDaughters = new TH1D("hdcaDaughters", "hdcaDaughters", 100, 0, 2);
    TH1D* hcosTheta = new TH1D("hcosTheta", "hcosTheta",100 , 0, 1);
    TH1D* hdcaD0ToPv = new TH1D("hdcaD0ToPv", "hdcaD0ToPv",100 , 0, 3);

    TH1D* hk_nSigma = new TH1D("hk_nSigma", "hk_nSigma",100 , -6, 6);
    TH1D* hpi1_nSigma = new TH1D("hpi1_nSigma", "hpi1_nSigma",100 , -6, 6);
    TH1D* hpi1_TOFinvbeta = new TH1D("hpi1_TOFinvbeta", "hpi1_TOFinvbeta",100 , -0.09, 0.09);
    TH1D* hk_TOFinvbeta = new TH1D("hk_TOFinvbeta", "hk_TOFinvbeta",100 , -0.09, 0.09);
    TH1D* hD_cosThetaStar = new TH1D("hD_cosThetaStar", "hD_cosThetaStar",100 , -1, 1);
    TH1D* hprimVzVpd = new TH1D("hprimVzVpd", "hprimVzVpd",100 , -40, 40);
    TH1D* hprimVz = new TH1D("hprimVz", "hprimVz",100 , -6, 6);
    TH1D* hk_nHitFit = new TH1D("hk_nHitFit", "hk_nHitFit",100 , 0, 45);
    TH1D* hpi1_nHitFit = new TH1D("hpi1_nHitFit", "hpi1_nHitFit",100 , 0, 45);


    TH1D* hvzvz = new TH1D("hvzvz", "hvzvz",100 , 0, 6);





    TH1D* hbk_pt = new TH1D("hbk_pt", "hbk_pt", 100, 0, 4);
    TH1D* hbpi1_pt = new TH1D("hbpi1_pt", "hbpi1_pt", 100, 0, 4);
    TH1D* hbD0_pt = new TH1D("hbD0_pt", "hbD0_pt", 100, 0, 7);
    TH1D* hbdecayLength = new TH1D("hbdecayLength", "hbdecayLength", 100, 0, 3);
    TH1D* hbpi1_dca = new TH1D("hbpi1_dca", "hbpi1_dca", 100, 0, 1.5);
    TH1D* hbk_dca = new TH1D("hbk_dca", "hbk_dca", 100, 0, 1.5);
    TH1D* hbdcaDaughters = new TH1D("hbdcaDaughters", "hbdcaDaughters", 100, 0, 2);
    TH1D* hbcosTheta = new TH1D("hbcosTheta", "hbcosTheta",100 , 0, 1);
    TH1D* hbdcaD0ToPv = new TH1D("hbdcaD0ToPv", "hbddcaD0ToPv",100 , 0, 3);

    TH1D* hbk_nSigma = new TH1D("hbk_nSigma", "hbk_nSigma",100 , -6, 6);
    TH1D* hbpi1_nSigma = new TH1D("hbpi1_nSigma", "hbpi1_nSigma",100 , -6, 6);
    TH1D* hbpi1_TOFinvbeta = new TH1D("hbpi1_TOFinvbeta", "hbpi1_TOFinvbeta",100 , -0.09, 0.09);
    TH1D* hbk_TOFinvbeta = new TH1D("hbk_TOFinvbeta", "hbk_TOFinvbeta",100 , -0.09, 0.09);
    TH1D* hbD_cosThetaStar = new TH1D("hbD_cosThetaStar", "hbD_cosThetaStar",100 , -1, 1);
    TH1D* hbprimVzVpd = new TH1D("hbprimVzVpd", "hbprimVzVpd",100 , -40, 40);
    TH1D* hbprimVz = new TH1D("hbprimVz", "hbprimVz",100 , -6, 6);
    TH1D* hbk_nHitFit = new TH1D("hbk_nHitFit", "hbk_nHitFit",100 , 0, 45);
    TH1D* hbpi1_nHitFit = new TH1D("hbpi1_nHitFit", "hbpi1_nHitFit",100 , 0, 45);

    TH1D* hbvzvz = new TH1D("hbvzvz", "hbvzvz",100 , 0, 6);

    TH1D* hVzVPD_VzTPC = new TH1D("hVzVPD_VzTPC", "hVzVPD_VzTPC",100 , -40, 40);
  //  hVzVPD_VzTPC->Add(hprimVzVpd,hprimVz,1,-1);

  // Dstar histogramy
    TH1D* hDstar_pt = new TH1D("hDstar_pt", "hDstar_pt", 100, 0, 7);


    TH1D* hInvMassSignDstar = new TH1D("Dstar signal", "Dstar signal", 32, 1.85, 2.2);
    TH1D* hInvMassBackDstar = new TH1D("Dstar background", "Dstar background", 32, 1.85, 2.2);

    TH1D* hInvMassSignDstarw = new TH1D("Dstar signalw", "Dstar signalw", 200, 0.5, 2.5);
    TH1D* hInvMassBackDstarw = new TH1D("Dstar backgroundw", "Dstar backgroundw", 200, 0.5, 2.5);

    TH1D* hInvMassSignDstar12 = new TH1D("Dstar S12", "Dstar S12", 32, 1.85, 2.2);
    TH1D* hInvMassSignDstar23 = new TH1D("Dstar S23", "Dstar S23", 32, 1.85, 2.2);
    TH1D* hInvMassSignDstar34 = new TH1D("Dstar S34", "Dstar S34", 32, 1.85, 2.2);
    TH1D* hInvMassSignDstar45 = new TH1D("Dstar S45", "Dstar S45", 32, 1.85, 2.2);

    TH1D* hInvMassBackDstar12 = new TH1D("bDstar B12", "bDstar B12", 32, 1.85, 2.2);
    TH1D* hInvMassBackDstar23 = new TH1D("bDstar B23", "bDstar B23", 32, 1.85, 2.2);
    TH1D* hInvMassBackDstar34 = new TH1D("bDstar B34", "bDstar B34", 32, 1.85, 2.2);
    TH1D* hInvMassBackDstar45 = new TH1D("bDstar B45", "bDstar B45", 32, 1.85, 2.2);

    TH1D* hInvMassSignDstar12w = new TH1D("Dstar S12w", "Dstar S12w", 200, 0.5, 2.5);
    TH1D* hInvMassSignDstar23w = new TH1D("Dstar S23w", "Dstar S23w", 200, 0.5, 2.5);
    TH1D* hInvMassSignDstar34w = new TH1D("Dstar S34w", "Dstar S34w", 200, 0.5, 2.5);
    TH1D* hInvMassSignDstar45w = new TH1D("Dstar S45w", "Dstar S45w", 200, 0.5, 2.5);

    TH1D* hInvMassBackDstar12w = new TH1D("bDstar B12w", "bDstar B12w", 200, 0.5, 2.5);
    TH1D* hInvMassBackDstar23w = new TH1D("bDstar B23w", "bDstar B23w", 200, 0.5, 2.5);
    TH1D* hInvMassBackDstar34w = new TH1D("bDstar B34w", "bDstar B34w", 200, 0.5, 2.5);
    TH1D* hInvMassBackDstar45w = new TH1D("bDstar B45w", "bDstar B45w", 200, 0.5, 2.5);



    TH1D* hDstar_D0 = new TH1D("Dstar-D0", "Dstar-D0", 32, 0.14, 0.16);
    TH1D* hbDstar_D0 = new TH1D("bDstar-D0", "bDstar-D0", 32, 0.14, 0.16);

    TH1D* hDstar_D012 = new TH1D("Dstar-D0 S12", "Dstar-D0 S12", 32, 0.14, 0.16);
    TH1D* hDstar_D023 = new TH1D("Dstar-D0 S23", "Dstar-D0 S23", 32, 0.14, 0.16);
    TH1D* hDstar_D034 = new TH1D("Dstar-D0 S34", "Dstar-D0 S34", 32, 0.14, 0.16);
    TH1D* hDstar_D045 = new TH1D("Dstar-D0 S45", "Dstar-D0 S45", 32, 0.14, 0.16);

    TH1D* hbDstar_D012 = new TH1D("bDstar-D0 B12", "bDstar-D0 B12", 32, 0.14, 0.16);
    TH1D* hbDstar_D023 = new TH1D("bDstar-D0 B23", "bDstar-D0 B23", 32, 0.14, 0.16);
    TH1D* hbDstar_D034 = new TH1D("bDstar-D0 B34", "bDstar-D0 B34", 32, 0.14, 0.16);
    TH1D* hbDstar_D045 = new TH1D("bDstar-D0 B45", "bDstar-D0 B45", 32, 0.14, 0.16);


    float Pion_nsigma = 3.0;
    float Kaon_nsigma = 2.0;
    float Pion_invbeta = 0.03;
    float Kaon_invbeta = 0.03;
    float Pion_nHits = 20;
    float Kaon_nHits = 20;
    float VzTPCVzVPD = 6;
    float cosThetaStarCut = 0.8;  //minimum

    hInvMassSignD0 -> Sumw2();
    hInvMassSignD0w -> Sumw2();
    hInvMassSignD012 -> Sumw2();
    hInvMassSignD012w -> Sumw2();
    hInvMassSignD023 -> Sumw2();
    hInvMassSignD023w -> Sumw2();
    hInvMassSignD034 -> Sumw2();
    hInvMassSignD034w -> Sumw2();
    hInvMassSignD045 -> Sumw2();
    hInvMassSignD045w -> Sumw2();
    hInvMassBackD0 -> Sumw2();
    hInvMassBackD0w -> Sumw2();
    hInvMassBackD012 -> Sumw2();
    hInvMassBackD012w -> Sumw2();
    hInvMassBackD023 -> Sumw2();
    hInvMassBackD023w -> Sumw2();
    hInvMassBackD034 -> Sumw2();
    hInvMassBackD034w -> Sumw2();
    hInvMassBackD045 -> Sumw2();
    hInvMassBackD045w -> Sumw2();
    hInvMassSignD0test -> Sumw2();
    hInvMassSignD012test -> Sumw2();
    hInvMassSignD023test -> Sumw2();
    hInvMassSignD034test -> Sumw2();
    hInvMassSignD045test -> Sumw2();
    hInvMassSignDstar -> Sumw2();
    hInvMassBackDstar -> Sumw2();
    hInvMassSignDstarw -> Sumw2();
    hInvMassBackDstarw -> Sumw2();
    hInvMassSignDstar12 -> Sumw2();
    hInvMassSignDstar23 -> Sumw2();
    hInvMassSignDstar34 -> Sumw2();
    hInvMassSignDstar45 -> Sumw2();
    hInvMassBackDstar12 -> Sumw2();
    hInvMassBackDstar23 -> Sumw2();
    hInvMassBackDstar34 -> Sumw2();
    hInvMassBackDstar45 -> Sumw2();
    hInvMassSignDstar12w -> Sumw2();
    hInvMassSignDstar23w -> Sumw2();
    hInvMassSignDstar34w -> Sumw2();
    hInvMassSignDstar45w -> Sumw2();
    hInvMassBackDstar12w -> Sumw2();
    hInvMassBackDstar23w -> Sumw2();
    hInvMassBackDstar34w -> Sumw2();
    hInvMassBackDstar45w -> Sumw2();
    hDstar_D0 -> Sumw2();
    hbDstar_D0 -> Sumw2();
    hDstar_D012 -> Sumw2();
    hDstar_D023 -> Sumw2();
    hDstar_D034 -> Sumw2();
    hDstar_D045 -> Sumw2();
    hbDstar_D012 -> Sumw2();
    hbDstar_D023 -> Sumw2();
    hbDstar_D034 -> Sumw2();
    hbDstar_D045 -> Sumw2();




    Long64_t equal = 0;
    Float_t dca_d0 = 9999;
    Long64_t numberEntr = ntp -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<<"Signal D0 "<<i<<endl;}
        ntp -> GetEntry(i);
       // if (k_eventId == pi1_eventId) equal++;
        if (cos(D0_theta)>0.0) { //cuty
            if ((D0_mass > 0.0) && (D0_mass < 4)) {
                float npi1_TOFinvbeta = pi1_TOFinvbeta / 0.012;
                float nk_TOFinvbeta = k_TOFinvbeta / 0.012;
                float f_res = pow(0.929095 + 0.0779541 / (k_p - 0.113628), 1.62916);  //sigma
                float f_pos = pow(-0.0538389 + 0.0439373 / (k_p - 0.0651247), 2.27704);  //mean

                float kaon_higher = 3 * f_res + f_pos;
                float kaon_lower = -2 * f_res + f_pos;

                float VzVPD_VzTPC = primVzVpd - primVz;
                hk_nSigma->Fill(k_nSigma);
                hpi1_nSigma->Fill(pi1_nSigma);
                hpi1_TOFinvbeta->Fill(pi1_TOFinvbeta);
                hk_TOFinvbeta->Fill(k_TOFinvbeta);
                hD_cosThetaStar->Fill(D_cosThetaStar);
                hprimVzVpd->Fill(primVzVpd);
                hprimVz->Fill(primVz);
                hk_nHitFit->Fill(k_nHitFit);
                hpi1_nHitFit->Fill(pi1_nHitFit);
                hvzvz->Fill(abs(primVz - primVzVpd));
                costhetastar->Fill(D0_pt, D_cosThetaStar);
                kaonnsigma->Fill(k_p, nk_TOFinvbeta);
                pionnsigma->Fill(pi1_p, npi1_TOFinvbeta);
                //  kaonnsigma50->Fill(k_p,nk_TOFinvbeta);
                //  pionnsigma50->Fill(pi1_p,npi1_TOFinvbeta);
                hVzVPD_VzTPC->Fill(VzVPD_VzTPC);


                if (abs(nk_TOFinvbeta) < 3) {
                    kaonnsigma50->Fill(k_p, nk_TOFinvbeta);
                }
                    if((abs(pi1_nSigma) < Pion_nsigma) && (abs(k_nSigma) < Kaon_nsigma)/*&& k_TOFinvbeta<kaon_higher && k_TOFinvbeta>kaon_lower (abs(k_TOFinvbeta)<Kaon_invbeta) && (abs(pi1_TOFinvbeta) < Pion_invbeta)*/&& (k_nHitFit > Kaon_nHits) && (pi1_nHitFit > Pion_nHits)&& (abs(primVz - primVzVpd) < VzTPCVzVPD) /*&& (D_cosThetaStar < cosThetaStarCut)*/){

                   // dca_d0 = D0_decayL * sqrt(1 - cos(D0_theta) * cos(D0_theta));

                    hInvMassSignD0->Fill(D0_mass);
                    hInvMassSignD0test->Fill(D0_mass);
                    hInvMassSignD0w->Fill(D0_mass);
                    hdecayLength->Fill(D0_decayL);
                    hpi1_dca->Fill(pi1_dca);
                    hk_dca->Fill(k_dca);
                    hdcaDaughters->Fill(dcaDaughters);
                    hdcaD0ToPv->Fill(dcaD0ToPv);
                    hk_pt->Fill(k_pt);
                    hpi1_pt->Fill(pi1_pt);
                    hD0_pt->Fill(D0_pt);



                    if ((D0_pt > 1) && (D0_pt < 2)) {
                    //  if (dca_d0 < 0.00493)
                    hInvMassSignD012->Fill(D0_mass);
                    hInvMassSignD012test->Fill(D0_mass);
                    hInvMassSignD012w->Fill(D0_mass);
                    }

                    if ((D0_pt > 2) && (D0_pt < 3) ) {
                    //  if (dca_d0 < 0.0046)
                    hInvMassSignD023->Fill(D0_mass);
                    hInvMassSignD023test->Fill(D0_mass);
                    hInvMassSignD023w->Fill(D0_mass);
                    }

                    if ((D0_pt > 3) && (D0_pt < 4) ) {
                    // if (dca_d0 < 0.00435)
                    hInvMassSignD034->Fill(D0_mass);
                    hInvMassSignD034test->Fill(D0_mass);
                    hInvMassSignD034w->Fill(D0_mass);
                    }

                    if (D0_pt > 4 ) {
                    // if (dca_d0 < 0.00435)
                    hInvMassSignD045->Fill(D0_mass);
                    hInvMassSignD045test->Fill(D0_mass);
                    hInvMassSignD045w->Fill(D0_mass);
                    }
                }
            }
        }
    }


  //  TH1D* kaonprojection = kaonnsigma->ProjectionY();Pion_invbeta
   // TH1D* pionprojection = pionnsigma->ProjectionY();



    kaonnsigma50->FitSlicesY();
    TH1D *sigmakaon = (TH1D*)gDirectory->Get("kaonnsigma50_2");
    TH1D *meankaon = (TH1D*)gDirectory->Get("kaonnsigma50_1");



    //  pionnsigma50->FitSlicesY();
   // TH1D *pionfit = (TH1D*)gDirectory->Get("pionnsigma50_2");





    numberEntr = bntp -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<< "Background D0 "<<i<<endl;}
        bntp -> GetEntry(i);
      //  if (k_eventId == bpi1_eventId) equal++;
        if (cos(bD0_theta)>0.0) { //cuty
            if ((bD0_mass > 0.0) && (bD0_mass < 4)){
                hbk_nSigma-> Fill(bk_nSigma);
                hbpi1_nSigma-> Fill(bpi1_nSigma);
                hbpi1_TOFinvbeta-> Fill(bpi1_TOFinvbeta);
                hbk_TOFinvbeta-> Fill(bk_TOFinvbeta);
                hbD_cosThetaStar-> Fill(bD_cosThetaStar);
                hbprimVzVpd-> Fill(bprimVzVpd);
                hbprimVz-> Fill(bprimVz);
                hbk_nHitFit-> Fill(bk_nHitFit);
                hbpi1_nHitFit-> Fill(bpi1_nHitFit);
                bcosthetastar->Fill(bD0_pt,bD_cosThetaStar);

                float bf_res = (0.929095+0.0779541/(bk_p -0.113628),1.62916);  //sigma
                float bf_pos = (-0.0538389+0.0439373/(bk_p -0.0651247),2.27704);  //mean

                float bkaon_higher = 3*bf_res + bf_pos;
                float bkaon_lower = -2*bf_res + bf_pos;

                if ((abs(bpi1_nSigma)<Pion_nsigma)&&(abs(bk_nSigma)<Kaon_nsigma) /*&& (abs(bpi1_TOFinvbeta) < Pion_invbeta)*/ && (bk_nHitFit>Kaon_nHits) && (bpi1_nHitFit>Pion_nHits) && (abs(bprimVz-bprimVzVpd)<VzTPCVzVPD)/*&& (bD_cosThetaStar < cosThetaStarCut)*/) {
                  //  dca_d0 = bD0_decayL * sqrt(1 - cos(bD0_theta) * cos(bD0_theta));

                    hInvMassBackD0->Fill(bD0_mass);
                    hInvMassBackD0w->Fill(bD0_mass);
                    hbdecayLength->Fill(bD0_decayL);
                    hbpi1_dca->Fill(bpi1_dca);
                    hbk_dca->Fill(bk_dca);
                    hbdcaDaughters->Fill(bdcaDaughters);
                    hbdcaD0ToPv->Fill(bdcaD0ToPv);
                    hbk_pt->Fill(bk_pt);
                    hbpi1_pt->Fill(bpi1_pt);
                    hbD0_pt->Fill(bD0_pt);


                    if ((bD0_pt > 1) && (bD0_pt < 2) ) {
                        //if (dca_d0 < 0.00493)
                        hInvMassBackD012->Fill(bD0_mass);
                        hInvMassBackD012w->Fill(bD0_mass);
                    }

                    if ((bD0_pt > 2) && (bD0_pt < 3) ) {
                        //if (dca_d0 < 0.0046)
                        hInvMassBackD023->Fill(bD0_mass);
                        hInvMassBackD023w->Fill(bD0_mass);
                    }

                    if ((bD0_pt > 3) && (bD0_pt < 4) ) {
                        //if (dca_d0 < 0.00435)
                        hInvMassBackD034->Fill(bD0_mass);
                        hInvMassBackD034w->Fill(bD0_mass);
                    }

                    if (bD0_pt > 4 ) {
                        //if (dca_d0 < 0.00435)
                        hInvMassBackD045->Fill(bD0_mass);
                        hInvMassBackD045w->Fill(bD0_mass);
                    }
                }
            }
        }
    }
   // cout<<"Number of equal runIds: "<<equal<<endl;





    Long64_t DnumberEntr = Dntp -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<DnumberEntr<<endl;
    for (Long64_t i = 0; i < DnumberEntr; i++) {
        if (i%10000000==0) {cout<<"Signal Dstar "<<i<<endl;}
        Dntp -> GetEntry(i);
        hDstar_pt->Fill(Dstar_pt);

        if ((Dstar_mass > 0.0) && (Dstar_mass < 4)&&(Dstar_pt > 2)) {
           // if(/*(abs(pi1_nSigma) < Pion_nsigma) && (abs(k_nSigma) < Kaon_nsigma)&& k_TOFinvbeta<kaon_higher && k_TOFinvbeta>kaon_lower (abs(k_TOFinvbeta)<Kaon_invbeta) && (abs(pi1_TOFinvbeta) < Pion_invbeta)&& (k_nHitFit > Kaon_nHits) && (pi1_nHitFit > Pion_nHits)&& (abs(primVz - primVzVpd) < VzTPCVzVPD) && (D_cosThetaStar < cosThetaStarCut)*/){

                    // dca_d0 = D0_decayL * sqrt(1 - cos(D0_theta) * cos(D0_theta));

                    hInvMassSignDstar->Fill(Dstar_mass);
                    hInvMassSignDstarw->Fill(Dstar_mass);
                    hDstar_D0->Fill(Dstar_D0);


                    if ((Dstar_pt > 1) && (Dstar_pt < 2)) {
                        hInvMassSignDstar12->Fill(Dstar_mass);
                        hInvMassSignDstar12w->Fill(Dstar_mass);
                        hDstar_D012->Fill(Dstar_D0);
                    }

                    if ((Dstar_pt > 2) && (Dstar_pt < 3) ) {
                        hInvMassSignDstar23->Fill(Dstar_mass);
                        hInvMassSignDstar23w->Fill(Dstar_mass);
                        hDstar_D023->Fill(Dstar_D0);
                    }

                    if ((Dstar_pt > 3) && (Dstar_pt < 4) ) {
                        hInvMassSignDstar34->Fill(Dstar_mass);
                        hInvMassSignDstar34w->Fill(Dstar_mass);
                        hDstar_D034->Fill(Dstar_D0);
                    }

                    if (Dstar_pt > 4 ) {
                        hInvMassSignDstar45->Fill(Dstar_mass);
                        hInvMassSignDstar45w->Fill(Dstar_mass);
                        hDstar_D045->Fill(Dstar_D0);
                    }
           // }
        }

    }


    DnumberEntr = Dbntp -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<DnumberEntr<<endl;
    for (Long64_t i = 0; i < DnumberEntr; i++) {
        if (i%10000000==0) {cout<< "Background Dstar "<<i<<endl;}
        Dbntp -> GetEntry(i);
        if ((bDstar_mass > 0.0) && (bDstar_mass < 4) && (bDstar_pt > 2)){

           // if (/*(abs(bpi1_nSigma)<Pion_nsigma)&&(abs(bk_nSigma)<Kaon_nsigma) && (abs(bpi1_TOFinvbeta) < Pion_invbeta) && (bk_nHitFit>Kaon_nHits) && (bpi1_nHitFit>Pion_nHits) && (abs(bprimVz-bprimVzVpd)<VzTPCVzVPD)&& (bD_cosThetaStar < cosThetaStarCut)*/) {
                    //  dca_d0 = bD0_decayL * sqrt(1 - cos(bD0_theta) * cos(bD0_theta));

                    hInvMassBackDstar->Fill(bDstar_mass);
                    hInvMassBackDstarw->Fill(bDstar_mass);
                    hbDstar_D0->Fill(bDstar_D0);



                    if ((bDstar_pt > 1) && (bDstar_pt < 2) ) {
                        hInvMassBackDstar12->Fill(bDstar_mass);
                        hInvMassBackDstar12w->Fill(bDstar_mass);
                        hbDstar_D012->Fill(bDstar_D0);
                    }

                    if ((bDstar_pt > 2) && (bDstar_pt < 3) ) {
                        hInvMassBackDstar23->Fill(bDstar_mass);
                        hInvMassBackDstar23w->Fill(bDstar_mass);
                        hbDstar_D023->Fill(bDstar_D0);
                    }

                    if ((bDstar_pt > 3) && (bDstar_pt < 4) ) {
                        hInvMassBackDstar34->Fill(bDstar_mass);
                        hInvMassBackDstar34w->Fill(bDstar_mass);
                        hbDstar_D034->Fill(bDstar_D0);
                    }

                    if (bDstar_pt > 4 ) {
                        hInvMassBackDstar45->Fill(bDstar_mass);
                        hInvMassBackDstar45w->Fill(bDstar_mass);
                        hbDstar_D045->Fill(bDstar_D0);
                    }
          //  }
        }

    }


    /*hInvMassBackPlus -> Clone("background");
    Double_t value, error, valueM, errorM, valueP, errorP;
    for (Int_t j = 1; j < 2001; j++) { //počítá chyby a spojí likesign pozadí plus a minus dokupy
        valueP = hInvMassBackPlus -> GetBinContent(j);
        errorP = hInvMassBackPlus -> GetBinError(j);
        valueM = hInvMassBackMin -> GetBinContent(j);
        errorM = hInvMassBackMin -> GetBinError(j);
        error = sqrt(valueM*errorP*errorP/valueP + valueP*errorM*errorM/valueM);
        value = 2*sqrt(valueP*valueM);

        hInvMassBack -> SetBinContent(j, value);
        hInvMassBack -> SetBinError(j, error);
    }*/


    TH1D *rozdilD0w = new TH1D("D0 rozdilw","D0 rozdilw",200, 0.5, 2.5);
    rozdilD0w->Add(hInvMassSignD0w,hInvMassBackD0w,1,-1);

    TH1D *rozdilD0 = new TH1D("D0 rozdil","D0 rozdil",32, 1.7, 2.05);
    rozdilD0->Add(hInvMassSignD0,hInvMassBackD0,1,-1);

    TH1D *rozdilD012 = new TH1D("D0 rozdil12","D0 rozdil12",32, 1.7, 2.05);
    rozdilD012->Add(hInvMassSignD012,hInvMassBackD012,1,-1);

    TH1D *rozdilD023 = new TH1D("D0 rozdil23","D0 rozdil23",32, 1.7, 2.05);
    rozdilD023->Add(hInvMassSignD023,hInvMassBackD023,1,-1);

    TH1D *rozdilD034 = new TH1D("D0 rozdil34","D0 rozdil34",32, 1.7, 2.05);
    rozdilD034->Add(hInvMassSignD034,hInvMassBackD034,1,-1);

    TH1D *rozdilD045 = new TH1D("D0 rozdil45","D0 rozdil45",32, 1.7, 2.05);
    rozdilD045->Add(hInvMassSignD045,hInvMassBackD045,1,-1);

    TH1D *rozdilD012w = new TH1D("D0 rozdil12w","D0 rozdil12w",200, 0.5, 2.5);
    rozdilD012w->Add(hInvMassSignD012w,hInvMassBackD012w,1,-1);

    TH1D *rozdilD023w = new TH1D("D0 rozdil23w","D0 rozdil23w",200, 0.5, 2.5);
    rozdilD023w->Add(hInvMassSignD023w,hInvMassBackD023w,1,-1);

    TH1D *rozdilD034w = new TH1D("D0 rozdil34w","D0 rozdil34w",200, 0.5, 2.5);
    rozdilD034w->Add(hInvMassSignD034w,hInvMassBackD034w,1,-1);

    TH1D *rozdilD045w = new TH1D("D0 rozdil45w","D0 rozdil45w",200, 0.5, 2.5);
    rozdilD045w->Add(hInvMassSignD045w,hInvMassBackD045w,1,-1);




    TH1D *rozdilDstarw = new TH1D("Dstar rozdilw","Dstar rozdilw",200, 0.5, 2.5);
    rozdilDstarw->Add(hInvMassSignDstarw,hInvMassBackDstarw,1,-0.3333333);

    TH1D *rozdilDstar = new TH1D("Dstar rozdil","Dstar rozdil",32, 1.85, 2.2);
    rozdilDstar->Add(hInvMassSignDstar,hInvMassBackDstar,1,-0.3333333);

    TH1D *rozdilDstar12 = new TH1D("Dstar rozdil12","Dstar rozdil12",32, 1.85, 2.2);
    rozdilDstar12->Add(hInvMassSignDstar12,hInvMassBackDstar12,1,-0.3333333);

    TH1D *rozdilDstar23 = new TH1D("Dstar rozdil23","Dstar rozdil23",32, 1.85, 2.2);
    rozdilDstar23->Add(hInvMassSignDstar23,hInvMassBackDstar23,1,-0.3333333);

    TH1D *rozdilDstar34 = new TH1D("Dstar rozdil34","Dstar rozdil34",32, 1.85, 2.2);
    rozdilDstar34->Add(hInvMassSignDstar34,hInvMassBackDstar34,1,-0.3333333);

    TH1D *rozdilDstar45 = new TH1D("Dstar rozdil45","Dstar rozdil45",32, 1.85, 2.2);
    rozdilDstar45->Add(hInvMassSignDstar45,hInvMassBackDstar45,1,-0.3333333);

    TH1D *rozdilDstar12w = new TH1D("Dstar rozdil12w","Dstar rozdil12w",200, 0.5, 2.5);
    rozdilDstar12w->Add(hInvMassSignDstar12w,hInvMassBackDstar12w,1,-0.3333333);

    TH1D *rozdilDstar23w = new TH1D("Dstar rozdil23w","Dstar rozdil23w",200, 0.5, 2.5);
    rozdilDstar23w->Add(hInvMassSignDstar23w,hInvMassBackDstar23w,1,-0.3333333);

    TH1D *rozdilDstar34w = new TH1D("Dstar rozdil34w","Dstar rozdil34w",200, 0.5, 2.5);
    rozdilDstar34w->Add(hInvMassSignDstar34w,hInvMassBackDstar34w,1,-0.3333333);

    TH1D *rozdilDstar45w = new TH1D("Dstar rozdil45w","Dstar rozdil45w",200, 0.5, 2.5);
    rozdilDstar45w->Add(hInvMassSignDstar45w,hInvMassBackDstar45w,1,-0.3333333);



    TH1D *rozdilDstar_D0 = new TH1D("Dstar-D0 rozdil","Dstar-D0 rozdil",32, 0.14, 0.16);
    rozdilDstar_D0->Add(hDstar_D0,hbDstar_D0,1,-0.3333333);

    TH1D *rozdilDstar_D012 = new TH1D("Dstar-D0 rozdil12","Dstar-D0 rozdil12",32, 0.14, 0.16);
    rozdilDstar_D012->Add(hDstar_D012,hbDstar_D012,1,-0.3333333);

    TH1D *rozdilDstar_D023 = new TH1D("Dstar-D0 rozdil23","Dstar-D0 rozdil23",32, 0.14, 0.16);
    rozdilDstar_D023->Add(hDstar_D023,hbDstar_D023,1,-0.3333333);

    TH1D *rozdilDstar_D034 = new TH1D("Dstar-D0 rozdil34","Dstar-D0 rozdil34",32, 0.14, 0.16);
    rozdilDstar_D034->Add(hDstar_D034,hbDstar_D034,1,-0.3333333);

    TH1D *rozdilDstar_D045 = new TH1D("Dstar-D0 rozdil45","Dstar-D0 rozdil45",32, 0.14, 0.16);
    rozdilDstar_D045->Add(hDstar_D045,hbDstar_D045,1,-0.3333333);





    reject = kTRUE;

    TF1 *fl = new TF1("fl",fline,1.7,2.05,2);
    hInvMassSignD0->Fit(fl,"0");

    TF1 *fq = new TF1("fq",fquadratic,1.7,2.05,3);
    hInvMassSignD0->Fit(fq,"0");

    TF1 *fc = new TF1("fc",fcubic,1.7,2.05,4);
    hInvMassSignD0->Fit(fc,"0");



    Double_t AA=0,BB=0;

    AA = fl->GetParameter(1); //parametr u x
    BB = fl->GetParameter(0);

    Double_t AAA=0, BBB=0, CCC=0;

    AAA = fq->GetParameter(0);
    BBB = fq->GetParameter(1); //parametr u x
    CCC = fq->GetParameter(2); //parametr u x^2

    Double_t AAAA=0,BBBB=0,CCCC=0,DDDD=0;

    AAAA = fc->GetParameter(0);
    BBBB = fc->GetParameter(1); //parametr u x
    CCCC = fc->GetParameter(2); //parametr u x^2
    DDDD = fc->GetParameter(3); //parametr u x^3



    TF1 *bckg_funct_linear = new TF1("bckg_funct_linear","[0]+[1]*x",1.7,2.05);
    bckg_funct_linear->SetParameters(BB,AA);

    TF1 *bckg_funct_quadratic = new TF1("bckg_funct_quadratic","[0]+[1]*x+[2]*x*x",1.7,2.05);
    bckg_funct_quadratic->SetParameters(AAA,BBB,CCC);

    TF1 *bckg_funct_cubic = new TF1("bckg_funct_cubic","[0]+[1]*x+[2]*x*x+[3]*x*x*x",1.7,2.05);
    bckg_funct_cubic->SetParameters(AAAA,BBBB,CCCC,DDDD);




    TH1D *odecet_linear = new TH1D("odecet_linear","odecet_linear",80, 1.7, 2.05);

    Double_t JJ=0, KK=0, LL=0, MM=0;

    for(int k=0; k<=80;k++){
        JJ=hInvMassSignD0->GetBinContent(k);
        MM=hInvMassSignD0->GetBinError(k);
        KK=odecet_linear->GetBinCenter(k);
        LL=bckg_funct_linear->TF1::Eval(KK);
        odecet_linear->SetBinContent(k,JJ-LL);
        odecet_linear->SetBinError(k,MM);
    }


    TH1D *odecet_quadratic = new TH1D("odecet_quadratic","odecet_quadratic",32, 1.7, 2.05);

    Double_t NN=0, OO=0, PP=0, QQ=0;

    for(int l=0; l<=32;l++){
        NN=hInvMassSignD0->GetBinContent(l);
        QQ=hInvMassSignD0->GetBinError(l);
        OO=odecet_quadratic->GetBinCenter(l);
        PP=bckg_funct_quadratic->TF1::Eval(OO);
        odecet_quadratic->SetBinContent(l,NN-PP);
        odecet_quadratic->SetBinError(l,QQ);
    }

    TH1D *odecet_cubic = new TH1D("odecet_cubic","odecet_cubic",32, 1.7, 2.05);

    Double_t RR=0, SS=0, TT=0, UU=0;

    for(int m=0; m<=32;m++){
        RR=hInvMassSignD0->GetBinContent(m);
        UU=hInvMassSignD0->GetBinError(m);
        SS=odecet_cubic->GetBinCenter(m);
        TT=bckg_funct_cubic->TF1::Eval(SS);
        odecet_cubic->SetBinContent(m,RR-TT);
        odecet_cubic->SetBinError(m,UU);

    }



    TFile* dataRes = new TFile("pp500_D0.root","RECREATE");


    hInvMassSignD0 -> Write();
    hInvMassSignD0w -> Write();
    hInvMassSignD0test -> Write();
    hInvMassBackD0 -> Write();
    hInvMassBackD0w -> Write();
    hInvMassBackD012 -> Write();
    hInvMassBackD023 -> Write();
    hInvMassBackD034 -> Write();
    hInvMassBackD045 -> Write();
    hInvMassSignD012 -> Write();
    hInvMassSignD023 -> Write();
    hInvMassSignD034 -> Write();
    hInvMassSignD045 -> Write();
    hInvMassBackD012w -> Write();
    hInvMassBackD023w -> Write();
    hInvMassBackD034w -> Write();
    hInvMassBackD045w -> Write();
    hInvMassSignD012w -> Write();
    hInvMassSignD023w -> Write();
    hInvMassSignD034w -> Write();
    hInvMassSignD045w -> Write();
    hInvMassSignD012test -> Write();
    hInvMassSignD023test -> Write();
    hInvMassSignD034test -> Write();
    hInvMassSignD045test -> Write();
    rozdilD0 -> Write();
    rozdilD012 -> Write();
    rozdilD023 -> Write();
    rozdilD034 -> Write();
    rozdilD045-> Write();
    rozdilD0w -> Write();
    rozdilD012w -> Write();
    rozdilD023w -> Write();
    rozdilD034w -> Write();
    rozdilD045w-> Write();
    /*
    hdecayLength -> Write();
    hpi1_dca -> Write();
    hk_dca -> Write();
    hdcaDaughters -> Write();
    hdcaD0ToPv -> Write();
    hk_pt -> Write();
    hpi1_pt -> Write();
    hbdecayLength -> Write();
    hbpi1_dca -> Write();
    hbk_dca -> Write();
    hbdcaDaughters -> Write();
    hbdcaD0ToPv -> Write();
    hbk_pt -> Write();
    hbpi1_pt -> Write();
    hbD0_pt -> Write();
    */
    hD0_pt -> Write();
    hk_nSigma-> Write();
    hpi1_nSigma-> Write();
    hpi1_TOFinvbeta-> Write();
    hk_TOFinvbeta-> Write();
    hD_cosThetaStar-> Write();
    hprimVzVpd-> Write();
    hprimVz-> Write();
    hk_nHitFit-> Write();
    hpi1_nHitFit-> Write();
    hvzvz->Write();
    hbvzvz->Write();

    hbk_nSigma-> Write();
    hbpi1_nSigma-> Write();
    hbpi1_TOFinvbeta-> Write();
    hbk_TOFinvbeta-> Write();
    hbD_cosThetaStar-> Write();
    hbprimVzVpd-> Write();
    hbprimVz-> Write();
    hbk_nHitFit-> Write();
    hbpi1_nHitFit-> Write();
    costhetastar->Write();
    bcosthetastar->Write();
    kaonnsigma->Write();
    pionnsigma->Write();
    hVzVPD_VzTPC->Write();


    sigmakaon->Write();
    meankaon->Write();

    odecet_linear->Write();
    odecet_quadratic->Write();
    odecet_cubic->Write();

    dataRes->Close();

    TFile* dataRes2 = new TFile("pp500_Dstar.root","RECREATE");



    hDstar_pt -> Write();
    hInvMassSignDstar -> Write();
    hInvMassSignDstarw -> Write();
    hInvMassBackDstar -> Write();
    hInvMassBackDstarw -> Write();
    hInvMassBackDstar12 -> Write();
    hInvMassBackDstar23 -> Write();
    hInvMassBackDstar34 -> Write();
    hInvMassBackDstar45 -> Write();
    hInvMassSignDstar12 -> Write();
    hInvMassSignDstar23 -> Write();
    hInvMassSignDstar34 -> Write();
    hInvMassSignDstar45 -> Write();
    hInvMassBackDstar12w -> Write();
    hInvMassBackDstar23w -> Write();
    hInvMassBackDstar34w -> Write();
    hInvMassBackDstar45w -> Write();
    hInvMassSignDstar12w -> Write();
    hInvMassSignDstar23w -> Write();
    hInvMassSignDstar34w -> Write();
    hInvMassSignDstar45w -> Write();
    rozdilDstar -> Write();
    rozdilDstar12 -> Write();
    rozdilDstar23 -> Write();
    rozdilDstar34 -> Write();
    rozdilDstar45-> Write();
    rozdilDstarw -> Write();
    rozdilDstar12w -> Write();
    rozdilDstar23w -> Write();
    rozdilDstar34w -> Write();
    rozdilDstar45w -> Write();
    hDstar_D0 -> Write();
    hDstar_D012 -> Write();
    hDstar_D023 -> Write();
    hDstar_D034 -> Write();
    hDstar_D045 -> Write();
    hbDstar_D0 -> Write();
    hbDstar_D012 -> Write();
    hbDstar_D023 -> Write();
    hbDstar_D034 -> Write();
    hbDstar_D045 -> Write();
    rozdilDstar_D0 -> Write();
    rozdilDstar_D012 -> Write();
    rozdilDstar_D023 -> Write();
    rozdilDstar_D034 -> Write();
    rozdilDstar_D045 -> Write();



    dataRes2->Close();


    // kaonprojection->Write();
  //  pionprojection->Write();



    cout<<"Michale jsi dobrej"<<endl;
    cout<<"Hotovo, Jarvis"<<endl;

}

