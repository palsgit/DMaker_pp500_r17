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
#include "TMath.h"

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

Double_t nSigmaTOFKaonTrack_upper(Double_t x) { 
	double f_nsigmaTOFKaonTrack_res =  1.21 + (0.10/(pow((x - 0.17), 1.21)));
    double f_nsigmaTOFKaonTrack_pos =  -0.01 + (0.02/(pow((x - 0.17), 1.73))); 

    /*float mSigma = 3.18;
    if (x > 2.0) mSigma = 2.41;*/
    float mSigma = 2.4;//midstep
	double kaon_higher = mSigma*f_nsigmaTOFKaonTrack_res + f_nsigmaTOFKaonTrack_pos;
    return kaon_higher; 
	}
Double_t nSigmaTOFKaonTrack_lower(Double_t x) {
	double f_nsigmaTOFKaonTrack_res =  1.21 + (0.10/(pow((x - 0.17), 1.21)));
    double f_nsigmaTOFKaonTrack_pos =  -0.01 + (0.02/(pow((x - 0.17), 1.73))); 
    /*float mSigma = -2.31;
    if (x > 1.25 && x < 1.50) mSigma = -1.56;
    if (x > 1.50 && x < 1.65) mSigma = -1.18;
    if (x > 1.65) mSigma = -0.8;*/
    float mSigma = -1.2; //midstep
	 double kaon_lower = mSigma*f_nsigmaTOFKaonTrack_res + f_nsigmaTOFKaonTrack_pos;
	 return kaon_lower; 
	 }

     Double_t nSigmaTPCPionTrack_upper(Double_t x) { 
	double f_nsigmaTPCPionTrack_res =  0.9548 + (0.0002/(pow((x + 0.1046), 3.9041)));
    double f_nsigmaTPCPionTrack_pos =  0.0988 + (-0.0012/(pow((x - 0.0706), 3.1974)));

    float mSigma = 2.92;
    //////if (x > 2.0) mSigma = 2.41;
	double kaon_higher = mSigma*f_nsigmaTPCPionTrack_res + f_nsigmaTPCPionTrack_pos;
    return kaon_higher; 
	}

Double_t nSigmaTPCPionTrack_lower(Double_t x) {
	double f_nsigmaTPCPionTrack_res =  0.9548 + (0.0002/(pow((x + 0.1046), 3.9041)));
    double f_nsigmaTPCPionTrack_pos =  0.0988 + (-0.0012/(pow((x - 0.0706), 3.1974))); 
    float mSigma = -2.70;
    /*if (x > 1.25 && x < 1.50) mSigma = -1.56;
    if (x > 1.50 && x < 1.65) mSigma = -1.18;
    if (x > 1.65) mSigma = -0.8;*/
	 double kaon_lower = mSigma*f_nsigmaTPCPionTrack_res + f_nsigmaTPCPionTrack_pos;
	 return kaon_lower; 
	 }



void Danalysis_dstar()
{
    TChain *ntp_RightSign = new TChain("Dntp_RightSign");
    ntp_RightSign->Add("output_all.root");
    TChain *ntp_SideBand = new TChain("Dntp_SideBand");
    ntp_SideBand->Add("output_all.root");
    /////TChain *ntp_ME = new TChain("ntp_signal_ME");
    /////ntp_ME->Add("output_mixedevent_sig.root");
    ////ntp_ME->Add("outputLocal.picoMEtree.sigME.root");
    TChain *ntp_WrongSign = new TChain("Dntp_Wrongsign");
    ntp_WrongSign->Add("output_all.root");


    ///Double_t pTforCos  = 0.0;
    /////Double_t cosThetaStarCutD0 = 1.0;
    //////Double_t cosThetaStarCutD0 = 0.90;
    //////Double_t cosThetaStarCutDstar = 0.77;
    //////Double_t cosThetapointD0cut = 0.1;
    ///Double_t dcaD0ToPvCut = 1.7;

    /*Int_t nSoftPionsTof, nSoftPionsBemc;

    Int_t nD0Candidates = 0;
    Int_t nSoftPionCandidates = 0;

    Int_t nSoftPions = 0;*/


    Float_t triplet_pairdiff, Dstar_rapidity, /*D0_decayL, D0_theta, cosTheta,*/ Dstar_pt, /*pi1_pt, k_pt,  pi1_dca, k_dca,*/ k_nSigmaTPC, pi1_nSigmaTPC, k_nSigmaTOF, pi1_nSigmaTOF, /*pi1_TOFinvbeta, k_TOFinvbeta, pi1_betaBase, k_betaBase,*/ D_cosThetaStar, /*dcaD0ToPv, dcaDaughters, primVz, primVzVpd, k_nHitFit, dcaMax, pi1_eventId, k_eventId,  pi1_nHitFit,*/ pi1_p, k_p, pi1_charge, k_charge, /*VzVPD_VzTPC, D_phi, D_eta,*/ pi2_charge, pi2_nSigmaTPC, pi2_nSigmaTOF,pi2_p;
    ntp_RightSign -> SetBranchAddress("triplet_pairdiff", &triplet_pairdiff);
    ntp_RightSign -> SetBranchAddress("Dstar_rapidity", &Dstar_rapidity);
    /*ntp_RightSign -> SetBranchAddress("D_decayL", &D0_decayL);
    ntp_RightSign -> SetBranchAddress("D_theta", &D0_theta);
    ntp_RightSign -> SetBranchAddress("cosTheta", &cosTheta);*/
    ntp_RightSign -> SetBranchAddress("Dstar_pt", &Dstar_pt);
    /*ntp_RightSign -> SetBranchAddress("pi1_pt", &pi1_pt);
    ntp_RightSign -> SetBranchAddress("k_pt", &k_pt);
    ntp_RightSign -> SetBranchAddress("pi1_dca", &pi1_dca);
    ntp_RightSign -> SetBranchAddress("k_dca", &k_dca);
    ntp_RightSign -> SetBranchAddress("pi1_TOFinvbeta", &pi1_TOFinvbeta);
    ntp_RightSign -> SetBranchAddress("k_TOFinvbeta", &k_TOFinvbeta);
    ntp_RightSign -> SetBranchAddress("pi1_betaBase", &pi1_betaBase);
    ntp_RightSign -> SetBranchAddress("k_betaBase", &k_betaBase);
    ntp_RightSign -> SetBranchAddress("D_cosThetaStar", &D_cosThetaStar);
    ntp_RightSign -> SetBranchAddress("dcaD0ToPv", &dcaD0ToPv);
    ntp_RightSign -> SetBranchAddress("dcaDaughters", &dcaDaughters);
    ntp_RightSign -> SetBranchAddress("primVz", &primVz);
    ntp_RightSign -> SetBranchAddress("primVzVpd", &primVzVpd);
    ntp_RightSign -> SetBranchAddress("k_nHitFit", &k_nHitFit);
    ntp_RightSign -> SetBranchAddress("pi1_nHitFit", &pi1_nHitFit);*/
    ntp_RightSign -> SetBranchAddress("D_cosThetaStar", &D_cosThetaStar);
    ntp_RightSign -> SetBranchAddress("pi1_p", &pi1_p);
    ntp_RightSign -> SetBranchAddress("pi2_p", &pi2_p);
    ntp_RightSign -> SetBranchAddress("k_nSigmaTPC", &k_nSigmaTPC);
    ntp_RightSign -> SetBranchAddress("pi1_nSigmaTPC", &pi1_nSigmaTPC);
    ntp_RightSign -> SetBranchAddress("k_nSigmaTOF", &k_nSigmaTOF);
    ntp_RightSign -> SetBranchAddress("pi1_nSigmaTOF", &pi1_nSigmaTOF);
    ntp_RightSign -> SetBranchAddress("pi2_nSigmaTPC", &pi2_nSigmaTPC);
    ntp_RightSign -> SetBranchAddress("pi2_nSigmaTOF", &pi2_nSigmaTOF);
    ntp_RightSign -> SetBranchAddress("k_p", &k_p);
    ntp_RightSign -> SetBranchAddress("pi1_charge", &pi1_charge);
    ntp_RightSign -> SetBranchAddress("k_charge", &k_charge);
    /*ntp_RightSign -> SetBranchAddress("VzVPD_VzTPC", &VzVPD_VzTPC);
    ntp_RightSign -> SetBranchAddress("D_phi", &D_phi);
    ntp_RightSign -> SetBranchAddress("D_eta", &D_eta);*/
    ntp_RightSign -> SetBranchAddress("pi2_charge", &pi2_charge);
    
    
    
    Float_t striplet_pairdiff, sDstar_rapidity, /*sD0_decayL, sD0_theta, scosTheta,*/  sDstar_pt, /*spi1_pt, sk_pt,  spi1_dca, sk_dca,*/ sk_nSigmaTPC, spi1_nSigmaTPC, sk_nSigmaTOF, spi1_nSigmaTOF,/* spi1_TOFinvbeta, sk_TOFinvbeta, spi1_betaBase, sk_betaBase,*/ sD_cosThetaStar, /*sdcaD0ToPv, sdcaDaughters, sprimVz, sprimVzVpd, sk_nHitFit, sdcaMax, spi1_eventId, sk_eventId,  spi1_nHitFit,*/ spi1_p, sk_p, spi1_charge, sk_charge, /*sVzVPD_VzTPC, sD_phi, sD_eta,*/ spi2_charge, spi2_nSigmaTPC, spi2_nSigmaTOF, spi2_p;
    ntp_SideBand -> SetBranchAddress("triplet_pairdiff", &striplet_pairdiff);
    ntp_SideBand -> SetBranchAddress("Dstar_rapidity", &sDstar_rapidity);
    /*ntp_SideBand -> SetBranchAddress("D_decayL", &sD0_decayL);
    ntp_SideBand -> SetBranchAddress("D_theta", &sD0_theta);
    ntp_SideBand -> SetBranchAddress("cosTheta", &scosTheta);*/
    ntp_SideBand -> SetBranchAddress("Dstar_pt", &sDstar_pt);
    /*ntp_SideBand -> SetBranchAddress("pi1_pt", &spi1_pt);
    ntp_SideBand -> SetBranchAddress("k_pt", &sk_pt);
    ntp_SideBand -> SetBranchAddress("pi1_dca", &spi1_dca);
    ntp_SideBand -> SetBranchAddress("k_dca", &sk_dca);
    ntp_SideBand -> SetBranchAddress("pi1_TOFinvbeta", &spi1_TOFinvbeta);
    ntp_SideBand -> SetBranchAddress("k_TOFinvbeta", &sk_TOFinvbeta);
    ntp_SideBand -> SetBranchAddress("pi1_betaBase", &spi1_betaBase);
    ntp_SideBand -> SetBranchAddress("k_betaBase", &sk_betaBase);
    ntp_SideBand -> SetBranchAddress("dcaD0ToPv", &sdcaD0ToPv);
    ntp_SideBand -> SetBranchAddress("dcaDaughters", &sdcaDaughters);
    ntp_SideBand -> SetBranchAddress("primVz", &sprimVz);
    ntp_SideBand -> SetBranchAddress("primVzVpd", &sprimVzVpd);
    ntp_SideBand -> SetBranchAddress("k_nHitFit", &sk_nHitFit);
    ntp_SideBand -> SetBranchAddress("pi1_nHitFit", &spi1_nHitFit);*/
    ntp_SideBand -> SetBranchAddress("D_cosThetaStar", &sD_cosThetaStar);
    ntp_SideBand -> SetBranchAddress("pi1_p", &spi1_p);
    ntp_SideBand -> SetBranchAddress("pi2_p", &spi2_p);
    ntp_SideBand -> SetBranchAddress("k_p", &sk_p);
    ntp_SideBand -> SetBranchAddress("k_nSigmaTPC", &sk_nSigmaTPC);
    ntp_SideBand -> SetBranchAddress("pi1_nSigmaTPC", &spi1_nSigmaTPC);
    ntp_SideBand -> SetBranchAddress("k_nSigmaTOF", &sk_nSigmaTOF);
    ntp_SideBand -> SetBranchAddress("pi1_nSigmaTOF", &spi1_nSigmaTOF);
    ntp_SideBand -> SetBranchAddress("pi2_nSigmaTPC", &spi2_nSigmaTPC);
    ntp_SideBand -> SetBranchAddress("pi2_nSigmaTOF", &spi2_nSigmaTOF);
    ntp_SideBand -> SetBranchAddress("pi1_charge", &spi1_charge);
    ntp_SideBand -> SetBranchAddress("k_charge", &sk_charge);
    /*ntp_SideBand -> SetBranchAddress("VzVPD_VzTPC", &sVzVPD_VzTPC);
    ntp_SideBand -> SetBranchAddress("D_phi", &sD_phi);
    ntp_SideBand -> SetBranchAddress("D_eta", &sD_eta);*/
    ntp_SideBand -> SetBranchAddress("pi2_charge", &spi2_charge);


    /*Float_t mtriplet_pairdiff, mDstar_rapidity, /*mD0_decayL, mD0_theta, mcosTheta, mDstar_pt, /*mpi1_pt, mk_pt,  mpi1_dca, mk_dca, mk_nSigmaTPC, mpi1_nSigmaTPC, mk_nSigmaTOF, mpi1_nSigmaTOF, mpi1_TOFinvbeta, mk_TOFinvbeta, mpi1_betaBase, mk_betaBase, mD_cosThetaStar, mdcaD0ToPv, mdcaDaughters, mprimVz, mprimVzVpd, mk_nHitFit, mdcaMax, mpi1_eventId, mk_eventId,  mpi1_nHitFit, mpi1_p, mk_p, mpi1_charge, mk_charge, /*mVzVPD_VzTPC, mD_phi, mD_eta, mpi2_charge;
    ntp_ME -> SetBranchAddress("triplet_pairdiff", &mtriplet_pairdiff);
    ntp_ME -> SetBranchAddress("Dstar_rapidity", &mDstar_rapidity);
    /*ntp_ME -> SetBranchAddress("D_decayL", &mD0_decayL);
    ntp_ME -> SetBranchAddress("D_theta", &mD0_theta);
    ntp_ME -> SetBranchAddress("cosTheta", &mcosTheta);
    ntp_ME -> SetBranchAddress("Dstar_pt", &mDstar_pt);
    /*ntp_ME -> SetBranchAddress("pi1_pt", &mpi1_pt);
    ntp_ME -> SetBranchAddress("k_pt", &mk_pt);
    ntp_ME -> SetBranchAddress("pi1_dca", &mpi1_dca);
    ntp_ME -> SetBranchAddress("k_dca", &mk_dca);
    ntp_ME -> SetBranchAddress("k_nSigmaTPC", &mk_nSigmaTPC);
    ntp_ME -> SetBranchAddress("pi1_nSigmaTPC", &mpi1_nSigmaTPC);
    ntp_ME -> SetBranchAddress("k_nSigmaTOF", &mk_nSigmaTOF);
    ntp_ME -> SetBranchAddress("pi1_nSigmaTOF", &mpi1_nSigmaTOF);
    ntp_ME -> SetBranchAddress("pi1_TOFinvbeta", &mpi1_TOFinvbeta);
    ntp_ME -> SetBranchAddress("k_TOFinvbeta", &mk_TOFinvbeta);
    ntp_ME -> SetBranchAddress("pi1_betaBase", &mpi1_betaBase);
    ntp_ME -> SetBranchAddress("k_betaBase", &mk_betaBase);
    ntp_ME -> SetBranchAddress("D_cosThetaStar", &mD_cosThetaStar);
    ntp_ME -> SetBranchAddress("dcaD0ToPv", &mdcaD0ToPv);
    ntp_ME -> SetBranchAddress("dcaDaughters", &mdcaDaughters);
    //////ntp_ME -> SetBranchAddress("primVz", &mprimVz);
    //////ntp_ME -> SetBranchAddress("primVzVpd", &mprimVzVpd);
    ntp_ME -> SetBranchAddress("k_nHitFit", &mk_nHitFit);
    ntp_ME -> SetBranchAddress("pi1_nHitFit", &mpi1_nHitFit);
    ntp_ME -> SetBranchAddress("pi1_p", &mpi1_p);
    ntp_ME -> SetBranchAddress("k_p", &mk_p);
    ntp_ME -> SetBranchAddress("pi1_charge", &mpi1_charge);
    ntp_ME -> SetBranchAddress("k_charge", &mk_charge);
    ///////ntp_ME -> SetBranchAddress("VzVPD_VzTPC", &mVzVPD_VzTPC);
    ntp_ME -> SetBranchAddress("D_phi", &mD_phi);
    ntp_ME -> SetBranchAddress("D_eta", &mD_eta);
    ntp_ME -> SetBranchAddress("pi2_charge", &mpi2_charge);*/

    Float_t wtriplet_pairdiff, wDstar_rapidity, /*wD0_decayL, wD0_theta, wcosTheta,*/ wDstar_pt, /*wpi1_pt, wk_pt,  wpi1_dca, wk_dca,*/ wk_nSigmaTPC, wpi1_nSigmaTPC, wk_nSigmaTOF, wpi1_nSigmaTOF,/* wpi1_TOFinvbeta, wk_TOFinvbeta, wpi1_betaBase, wk_betaBase,*/ wD_cosThetaStar, /*wdcaD0ToPv, wdcaDaughters, wprimVz, wprimVzVpd, wk_nHitFit, wdcaMax, wpi1_eventId, wk_eventId,  wpi1_nHitFit,*/ wpi1_p, wk_p, wpi1_charge, wk_charge, /*wVzVPD_VzTPC, wD_phi, wD_eta,*/ wpi2_charge, wpi2_nSigmaTPC, wpi2_nSigmaTOF, wpi2_p;
    ntp_WrongSign -> SetBranchAddress("triplet_pairdiff", &wtriplet_pairdiff);
    ntp_WrongSign -> SetBranchAddress("Dstar_rapidity", &wDstar_rapidity);
    /*ntp_WrongSign -> SetBranchAddress("D_decayL", &wD0_decayL);
    ntp_WrongSign -> SetBranchAddress("D_theta", &wD0_theta);
    ntp_WrongSign -> SetBranchAddress("cosTheta", &wcosTheta);*/
    ntp_WrongSign -> SetBranchAddress("Dstar_pt", &wDstar_pt);
    /*ntp_WrongSign -> SetBranchAddress("pi1_pt", &wpi1_pt);
    ntp_WrongSign -> SetBranchAddress("k_pt", &wk_pt);
    ntp_WrongSign -> SetBranchAddress("pi1_dca", &wpi1_dca);
    ntp_WrongSign -> SetBranchAddress("k_dca", &wk_dca);
    ntp_WrongSign -> SetBranchAddress("pi1_TOFinvbeta", &wpi1_TOFinvbeta);
    ntp_WrongSign -> SetBranchAddress("k_TOFinvbeta", &wk_TOFinvbeta);
    ntp_WrongSign -> SetBranchAddress("pi1_betaBase", &wpi1_betaBase);
    ntp_WrongSign -> SetBranchAddress("k_betaBase", &wk_betaBase);
    ntp_WrongSign -> SetBranchAddress("D_cosThetaStar", &wD_cosThetaStar);
    ntp_WrongSign -> SetBranchAddress("dcaD0ToPv", &wdcaD0ToPv);
    ntp_WrongSign -> SetBranchAddress("dcaDaughters", &wdcaDaughters);
    ntp_WrongSign -> SetBranchAddress("primVz", &wprimVz);
    ntp_WrongSign -> SetBranchAddress("primVzVpd", &wprimVzVpd);
    ntp_WrongSign -> SetBranchAddress("k_nHitFit", &wk_nHitFit);
    ntp_WrongSign -> SetBranchAddress("pi1_nHitFit", &wpi1_nHitFit);*/
    ntp_WrongSign -> SetBranchAddress("D_cosThetaStar", &wD_cosThetaStar);
    ntp_WrongSign -> SetBranchAddress("pi1_p", &wpi1_p);
     ntp_WrongSign -> SetBranchAddress("pi2_p", &wpi2_p);
    ntp_WrongSign -> SetBranchAddress("k_nSigmaTPC", &wk_nSigmaTPC);
    ntp_WrongSign -> SetBranchAddress("pi1_nSigmaTPC", &wpi1_nSigmaTPC);
    ntp_WrongSign -> SetBranchAddress("k_nSigmaTOF", &wk_nSigmaTOF);
    ntp_WrongSign -> SetBranchAddress("pi1_nSigmaTOF", &wpi1_nSigmaTOF);
     ntp_SideBand -> SetBranchAddress("pi2_nSigmaTPC", &wpi2_nSigmaTPC);
    ntp_SideBand -> SetBranchAddress("pi2_nSigmaTOF", &wpi2_nSigmaTOF);
    ntp_WrongSign -> SetBranchAddress("k_p", &wk_p);
    ntp_WrongSign -> SetBranchAddress("pi1_charge", &wpi1_charge);
    ntp_WrongSign -> SetBranchAddress("k_charge", &wk_charge);
    /*ntp_WrongSign -> SetBranchAddress("VzVPD_VzTPC", &wVzVPD_VzTPC);
    ntp_WrongSign -> SetBranchAddress("D_phi", &wD_phi);
    ntp_WrongSign -> SetBranchAddress("D_eta", &wD_eta);*/
    ntp_WrongSign -> SetBranchAddress("pi2_charge", &wpi2_charge);
    

    TH2F *DstarPlusRightSign = new TH2F("DstarPlusRightSign","DeltaM of DstarPlus RightSign",100,0,10.,4000,0.05,2.05);
    TH2F *DstarMinusRightSign = new TH2F("DstarMinusRightSign","DeltaM of DstarMinus RightSign",100,0,10.,4000,0.05,2.05); 
    TH2F *DstarPlusSideBand = new TH2F("DstarPlusSideBand","DeltaM of DstarPlus SideBand",100,0,10.,4000,0.05,2.05);
    TH2F *DstarMinusSideBand = new TH2F("DstarMinusSideBand","DeltaM of DstarMinus SideBand",100,0,10.,4000,0.05,2.05);
    TH2F *MixedDstarPlusRightSign = new TH2F("MixedDstarPlusRightSign","Mixed event background",100,0,10.,4000,0.05,2.05);
    TH2F *MixedDstarMinusRightSign = new TH2F("MixedDstarMinusRightSign","Mixed event background",100,0,10.,4000,0.05,2.05);
    TH2F *DstarPlusWrongSign = new TH2F("DstarPlusWrongSign","DeltaM of DstarPlus WrongSign",100,0,10.,4000,0.05,2.05); 
    TH2F *DstarMinusWrongSign = new TH2F("DstarMinusWrongSign","DeltaM of DstarMinus WrongSign",100,0,10.,4000,0.05,2.05);


    /*
    TH1D *dcaDaughters_unlike = new TH1D("dcaDaughters_unlike","dcaDaughters; dca[cm]; Entries",300,-1.0, 5.0);
    TH1D *D_theta_unlike = new TH1D("D_theta_unlike","D_theta; #theta_{point} [rad]; Entries",300,-1.0, 4.0);
    TH1D *cosTheta_unlike = new TH1D("cosTheta_unlike","cosTheta; cos #theta_{point}; Entries",300,-2.0, 2.0);
    TH1D *D_decayL_unlike = new TH1D("D_decayL_unlike","decay length; decay length [cm]; Entries",2000,-1.0, 1000.0);
    TH1D *dcaD0ToPv_unlike = new TH1D("dcaD0ToPv_unlike","dcaD0ToPv; dcaD0ToPv [cm]; Entries",3000,-1.0, 900.0);
    TH1D *D_cosThetaStar_unlike = new TH1D("D_cosThetaStar_unlike","D_cosThetaStar; cos #theta^{*}; Entries",300,-2.0, 2.0);
    TH1D *Dstar_pt_unlike = new TH1D("Dstar_pt_unlike","Dstar_pt; p_{T} [GeV]; Entries",100,0,10.);
    TH1D *triplet_pairdiff_unlike = new TH1D("triplet_pairdiff_unlike","triplet_pairdiff; invariant mass [GeV]; Entries",4000,-0.0,2.0);
    TH1D *Dstar_rapidity_unlike = new TH1D("Dstar_rapidity_unlike","Dstar_rapidity; #pi K pair rapidity; Entries",3000,-3.0, 3.0);
    TH1D *D_phi_unlike = new TH1D("D_phi_unlike","D_phi; #pi K pair phi; Entries",100,-5.0, 5.0);
    TH1D *D_eta_unlike = new TH1D("D_eta_unlike","D_eta; #pi K pair eta; Entries",4000,-20.0, 20.0);
    TH1D *pi_DCA_unlike = new TH1D("pi_DCA_unlike","pi_DCA; DCA #pi [cm] ; Entries",100,-1.0, 5.0);
    TH1D *k_DCA_unlike = new TH1D("k_DCA_unlike","k_DCA;  DCA k [cm]; Entries",100,-1.0, 5.0);
    TH2D *D_y_D_eta_unlike = new TH2D("D_y_D_eta_unlike","D_{y}_D_{#eta}_unlike; D_{y}; D_{#eta}",2000,-10.,10.,2000,-10.,10.);


    TH1D *dcaDaughters_like = new TH1D("dcaDaughters_like","dcaDaughters; dca[cm]; Entries",300,-1.0, 5.0);
    TH1D *D_theta_like = new TH1D("D_theta_like","D_theta; #theta_{point} [rad]; Entries",300,-1.0, 4.0);
    TH1D *cosTheta_like = new TH1D("cosTheta_like","cosTheta; cos #theta_{point}; Entries",300,-2.0, 2.0);
    TH1D *D_decayL_like = new TH1D("D_decayL_like","decay length; decay length [cm]; Entries",2000,-1.0, 1000.0);
    TH1D *dcaD0ToPv_like = new TH1D("dcaD0ToPv_like","dcaD0ToPv; dcaD0ToPv [cm]; Entries",3000,-1.0, 900.0);
    TH1D *D_cosThetaStar_like = new TH1D("D_cosThetaStar_like","D_cosThetaStar; cos #theta^{*}; Entries",300,-2.0, 2.0);
    TH1D *Dstar_pt_like = new TH1D("Dstar_pt_like","Dstar_pt; p_{T} [GeV]; Entries",100,0,10.);
    TH1D *triplet_pairdiff_like = new TH1D("triplet_pairdiff_like","triplet_pairdiff; invariant mass [GeV]; Entries",4000,-0.0,2.0);
    TH1D *Dstar_rapidity_like = new TH1D("Dstar_rapidity_like","Dstar_rapidity; #pi K pair rapidity; Entries",3000,-3.0, 3.0);
    TH1D *D_phi_like = new TH1D("D_phi_like","D_phi; #pi K pair phi; Entries",100,-5.0, 5.0);
    TH1D *D_eta_like = new TH1D("D_eta_like","D_eta; #pi K pair eta; Entries",4000,-20.0, 20.0);
    TH1D *pi_DCA_like = new TH1D("pi_DCA_like","pi_DCA; DCA #pi [cm] ; Entries",100,-1.0, 5.0);
    TH1D *k_DCA_like = new TH1D("k_DCA_like","k_DCA;  DCA k [cm]; Entries",100,-1.0, 5.0);
    TH2D *D_y_D_eta_like = new TH2D("D_y_D_eta_like","D_{y}_D_{#eta}_like; D_{y}; D_{#eta}",2000,-10.,10.,2000,-10.,10.);


    TH1D *dcaDaughters_rot = new TH1D("dcaDaughters_rot","dcaDaughters; dca[cm]; Entries",300,-1.0, 5.0);
    TH1D *D_theta_rot = new TH1D("D_theta_rot","D_theta; #theta_{point} [rad]; Entries",300,-1.0, 4.0);
    TH1D *cosTheta_rot = new TH1D("cosTheta_rot","cosTheta; cos #theta_{point}; Entries",300,-2.0, 2.0);
    TH1D *D_decayL_rot = new TH1D("D_decayL_rot","decay length; decay length [cm]; Entries",2000,-1.0, 1000.0);
    TH1D *dcaD0ToPv_rot = new TH1D("dcaD0ToPv_rot","dcaD0ToPv; dcaD0ToPv [cm]; Entries",3000,-1.0, 900.0);
    TH1D *D_cosThetaStar_rot = new TH1D("D_cosThetaStar_rot","D_cosThetaStar; cos #theta^{*}; Entries",300,-2.0, 2.0);
    TH1D *Dstar_pt_rot = new TH1D("Dstar_pt_rot","Dstar_pt; p_{T} [GeV]; Entries",100,0,10.);
    TH1D *triplet_pairdiff_rot = new TH1D("triplet_pairdiff_rot","triplet_pairdiff; invariant mass [GeV]; Entries",4000,-0.0,2.0);
    TH1D *Dstar_rapidity_rot = new TH1D("Dstar_rapidity_rot","Dstar_rapidity; #pi K pair rapidity; Entries",3000,-3.0, 3.0);
    TH1D *D_phi_rot = new TH1D("D_phi_rot","D_phi; #pi K pair phi; Entries",100,-5.0, 5.0);
    TH1D *D_eta_rot = new TH1D("D_eta_rot","D_eta; #pi K pair eta; Entries",4000,-20.0, 20.0);
    TH1D *pi_DCA_rot = new TH1D("pi_DCA_rot","pi_DCA; DCA #pi [cm] ; Entries",100,-1.0, 5.0);
    TH1D *k_DCA_rot = new TH1D("k_DCA_rot","k_DCA;  DCA k [cm]; Entries",100,-1.0, 5.0);
    TH2D *D_y_D_eta_rot = new TH2D("D_y_D_eta_rot","D_{y}_D_{#eta}_rot; D_{y}; D_{#eta}",2000,-10.,10.,2000,-10.,10.);

    TH1D *dcaDaughters_mixed = new TH1D("dcaDaughters_mixed","dcaDaughters; dca[cm]; Entries",300,-1.0, 5.0);
    TH1D *D_theta_mixed = new TH1D("D_theta_mixed","D_theta; #theta_{point} [rad]; Entries",300,-1.0, 4.0);
    TH1D *cosTheta_mixed = new TH1D("cosTheta_mixed","cosTheta; cos #theta_{point}; Entries",300,-2.0, 2.0);
    TH1D *D_decayL_mixed = new TH1D("D_decayL_mixed","decay length; decay length [cm]; Entries",2000,-1.0, 1000.0);
    TH1D *dcaD0ToPv_mixed = new TH1D("dcaD0ToPv_mixed","dcaD0ToPv; dcaD0ToPv [cm]; Entries",3000,-1.0, 900.0);
    TH1D *D_cosThetaStar_mixed = new TH1D("D_cosThetaStar_mixed","D_cosThetaStar; cos #theta^{*}; Entries",300,-2.0, 2.0);
    TH1D *Dstar_pt_mixed = new TH1D("Dstar_pt_mixed","Dstar_pt; p_{T} [GeV]; Entries",100,0,10.);
    TH1D *triplet_pairdiff_mixed = new TH1D("triplet_pairdiff_mixed","triplet_pairdiff; invariant mass [GeV]; Entries",4000,-0.0,2.0);
    TH1D *Dstar_rapidity_mixed = new TH1D("Dstar_rapidity_mixed","Dstar_rapidity; #pi K pair rapidity; Entries",3000,-3.0, 3.0);
    TH1D *D_phi_mixed = new TH1D("D_phi_mixed","D_phi; #pi K pair phi; Entries",100,-5.0, 5.0);
    TH1D *D_eta_mixed = new TH1D("D_eta_mixed","D_eta; #pi K pair eta; Entries",4000,-20.0, 20.0);
    TH1D *pi_DCA_mixed = new TH1D("pi_DCA_mixed","pi_DCA; DCA #pi [cm] ; Entries",100,-1.0, 5.0);
    TH1D *k_DCA_mixed = new TH1D("k_DCA_mixed","k_DCA;  DCA k [cm]; Entries",100,-1.0, 5.0);
    TH2D *D_y_D_eta_mixed = new TH2D("D_y_D_eta_mixed","D_{y}_D_{#eta}_mixed; D_{y}; D_{#eta}",2000,-10.,10.,2000,-10.,10.);

    /*TH2F *hRapidityD = new TH2F("RapidityD","#pi K pair rapidity",100,0.,5.,100,-5.,5.);
    TH2F *hDptPiPt = new TH2F("hDptPiPt","X:D^{0} p_{T}, Y:#pi p_{T}",200,0.,10.,200,0.,10.);
    TH2F *hDptKPt = new TH2F("hDptKPt","X:D^{0} p_{T}, Y:K p_{T}",200,0.,10.,200,0.,10.);
    TH2F *hDptPiP = new TH2F("hDptPiP","X:D^{0} p_{T}, Y:#pi momentum",200,0.,10.,200,0.,10.);
    TH2F *hDptKP = new TH2F("hDptKP","X:D^{0} p_{T}, Y:K momentum",200,0.,10.,200,0.,10.);
    TH2F *hPtPionVsKaon = new TH2F("hPtPionVsKaon","hPtPionVsKaon",500,0.,5.,500,0.,5.);
    TH2F *hNHitsKvsD0pt = new TH2F("NHitsKvsD0pt","Number of hits of Kaon Track vs D0 p_{#perp}",200,0.,10.,50,0,50);
    TH2F *hNHitsPivsD0pt = new TH2F("NHitsPivsD0pt","Number of hits of Pion Track vs D0 p_{#perp}",200,0.,10.,50,0,50);
    TH2F *hDstarRapidity = new TH2F("DstarRapidity","D* rapidity",500,0.,5.,100,-2.,2.);
    TH2F *hSoftPionPlus = new TH2F("SoftPionPlus","Invariant mass of SoftPion plus ",200,0,20.,4000,0.135,0.175); 
    TH2F *hSoftPionMinus = new TH2F("SoftPionMinus","Invariant mass of SoftPion minus ",200,0,20.,4000,0.135,0.175);
    TH2F *hD0PtVsSoftPionPt = new TH2F("D0PtVsSoftPionPt","p_{T}(D^{0}) vs p_{T}(soft #pi)",40,0.,2.,4000,0.,20.);
    TH2F *hSideBand = new TH2F("SideBand","Side band ",200,0,20.,4000,0.135,0.175);
    TH2F *hWrongSign = new TH2F("WrongSign","Wrong Sign",200,0,20.,4000,0.135,0.175);
    TH2I *hD0SoftPionCandidatesCorr = new TH2I("D0SoftPionCandidatesCorr","(Soft #pi candidates, D^{0} candidates)",10,0,10,10,0,10);
    TH1F *hVertexZ = new TH1F("VertexZ","Vertex Z",500,-250.,250);*/
    
    Long64_t numberEntr = ntp_RightSign -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<<"RightSign Dstar "<<i<<endl;}
        ntp_RightSign -> GetEntry(i);
        ///if (Dstar_pt<=0) continue;
        ///////if ((Dstar_pt<1.1) || (Dstar_pt>2.1))  continue;
        if(TMath::Abs(Dstar_rapidity)>=1.0) continue;
        if(k_nSigmaTPC <= -2.5 || k_nSigmaTPC >= 3.0) continue; 
        if(k_nSigmaTOF <= nSigmaTOFKaonTrack_lower(k_p) || k_nSigmaTOF >= nSigmaTOFKaonTrack_upper(k_p)) continue;

        if(pi1_nSigmaTPC <= -3.0 || pi1_nSigmaTPC >= 3.0) continue;
        if(pi1_nSigmaTOF <= -2.5 || pi1_nSigmaTOF >= 2.5) continue;

        if(pi2_nSigmaTPC <= -3.0 || pi2_nSigmaTPC >= 3.0) continue;
        if(pi2_nSigmaTOF <= -2.5 || pi2_nSigmaTOF >= 2.5) continue;
        
        //////if(D0_decayL < 1.0) continue;  ///could be strict
        /*if(k_nSigmaTPC < -2.5 || k_nSigmaTPC > 2.5) continue; 
        if(pi1_nSigmaTOF < -2.5 || pi1_nSigmaTOF > 2.5) continue;
        if (k_nSigmaTOF < nSigmaTOFKaonTrack_lower(k_p) || k_nSigmaTOF > nSigmaTOFKaonTrack_upper(k_p)) continue;
        */
        //////if(dcaD0ToPv > 99.9) continue;  ///could be strict
        ///////if(TMath::Abs(cosTheta) > 0.90) continue; ////come back here
        ///////if(TMath::Abs(D_eta) > 1.0) continue;
        
        if(D_cosThetaStar >= 0.7) continue;  ////come back here
        /////if(pi1_dca > 0.62) continue;  ///????
        /////if(k_dca < 0.38) continue;
        
        if(pi1_charge == 1 && k_charge == -1 && pi2_charge == 1) {
            ////if(((triplet_pairdiff>1.75) && (triplet_pairdiff<1.85)) || ((triplet_pairdiff>1.9) && (triplet_pairdiff<2.0))){
            /////if(((triplet_pairdiff>1.78) && (triplet_pairdiff<1.83))){
               /*nD0Candidates++;
               hNHitsKvsD0pt->Fill(Dstar_pt,k_nHitFit);
			   hNHitsPivsD0pt->Fill(Dstar_pt,pi1_nHitFit);
               hDptPiPt->Fill(Dstar_pt,pi1_pt);
               hDptKPt->Fill(Dstar_pt,k_pt);
               hDptPiP->Fill(Dstar_pt,pi1_p);
			   hDptKP->Fill(Dstar_pt,k_p);
               hPtPionVsKaon->Fill(pi1_pt,k_pt);*/
               /*dcaDaughters_unlike->Fill(dcaDaughters);
               D_theta_unlike->Fill(D0_theta);
               cosTheta_unlike->Fill(cosTheta);
               D_decayL_unlike->Fill(D0_decayL);
               dcaD0ToPv_unlike->Fill(dcaD0ToPv);
               D_cosThetaStar_unlike->Fill(D_cosThetaStar);
               Dstar_pt_unlike->Fill(Dstar_pt);
               triplet_pairdiff_unlike->Fill(triplet_pairdiff);
               Dstar_rapidity_unlike->Fill(Dstar_rapidity);
               D_phi_unlike->Fill(D_phi);
               D_eta_unlike->Fill(D_eta);
               pi_DCA_unlike->Fill(pi1_dca);
               k_DCA_unlike->Fill(k_dca);
               D_y_D_eta_unlike->Fill(Dstar_rapidity, D_eta);*/
            /////}
            /////if(Dstar_pt>pTforCos) {
            //////if((D_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(D0_theta)) > cosThetapointD0cut)) {
                ///////if((D_cosThetaStar < cosThetaStarCutD0)){
                    DstarPlusRightSign->Fill(Dstar_pt,triplet_pairdiff);
                    //////hRapidityD->Fill(Dstar_pt,Dstar_rapidity);
            //////}               
            //////}
            /*else{
					DstarPlusRightSign->Fill(Dstar_pt,triplet_pairdiff);	
                    hRapidityD->Fill(Dstar_pt,Dstar_rapidity);
				}
                */
        }

        if(pi1_charge == -1 && k_charge == 1 && pi2_charge == -1) {
             ////if(((triplet_pairdiff>1.75) && (triplet_pairdiff<1.85)) || ((triplet_pairdiff>1.9) && (triplet_pairdiff<2.0))){
             ///////if(((triplet_pairdiff>1.78) && (triplet_pairdiff<1.83))){
               /*hNHitsKvsD0pt->Fill(Dstar_pt,k_nHitFit);
			   hNHitsPivsD0pt->Fill(Dstar_pt,pi1_nHitFit);
               hDptPiPt->Fill(Dstar_pt,pi1_pt);
               hDptKPt->Fill(Dstar_pt,k_pt);
               hDptPiP->Fill(Dstar_pt,pi1_p);
			   hDptKP->Fill(Dstar_pt,k_p);
               hPtPionVsKaon->Fill(pi1_pt,k_pt);
               */

              /*dcaDaughters_unlike->Fill(dcaDaughters);
               D_theta_unlike->Fill(D0_theta);
               cosTheta_unlike->Fill(cosTheta);
               D_decayL_unlike->Fill(D0_decayL);
               dcaD0ToPv_unlike->Fill(dcaD0ToPv);
               D_cosThetaStar_unlike->Fill(D_cosThetaStar);
               Dstar_pt_unlike->Fill(Dstar_pt);
               triplet_pairdiff_unlike->Fill(triplet_pairdiff);
               Dstar_rapidity_unlike->Fill(Dstar_rapidity);
               D_phi_unlike->Fill(D_phi);
               D_eta_unlike->Fill(D_eta);
               pi_DCA_unlike->Fill(pi1_dca);
               k_DCA_unlike->Fill(k_dca);
               D_y_D_eta_unlike->Fill(Dstar_rapidity, D_eta);*/
            //////}
            //////if(Dstar_pt>pTforCos) {
            //////if((D_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(D0_theta)) > cosThetapointD0cut)) {
            ///////if((D_cosThetaStar < cosThetaStarCutD0)){
                    DstarMinusRightSign->Fill(Dstar_pt,triplet_pairdiff);
                    ////hRapidityD->Fill(Dstar_pt,Dstar_rapidity);
            ///////}               
            /////}
            /*else{
					DstarMinusRightSign->Fill(Dstar_pt,triplet_pairdiff);	
                    hRapidityD->Fill(Dstar_pt,Dstar_rapidity);
				}
                */
        }

    }


    numberEntr = ntp_SideBand -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<< "Sideband Dstar "<<i<<endl;}
        ntp_SideBand -> GetEntry(i);
        ///if (Dstar_pt<=0) continue;
        ////if ((sDstar_pt<1.1) || (sDstar_pt>2.1))  continue;
        if((TMath::Abs(sDstar_rapidity)>= 1.0)) continue;

        if(sk_nSigmaTPC <= -2.5 || sk_nSigmaTPC >= 3.0) continue; 
        if(sk_nSigmaTOF <= nSigmaTOFKaonTrack_lower(sk_p) || sk_nSigmaTOF >= nSigmaTOFKaonTrack_upper(sk_p)) continue;

        if(spi1_nSigmaTPC <= -3.0 || spi1_nSigmaTPC >= 3.0) continue;
        if(spi1_nSigmaTOF <= -2.5 || spi1_nSigmaTOF >= 2.5) continue;

        if(spi2_nSigmaTPC <= -3.0 || spi2_nSigmaTPC >= 3.0) continue;
        if(spi2_nSigmaTOF <= -2.5 || spi2_nSigmaTOF >= 2.5) continue;
        
        //////if(sD0_decayL < 0.0 || sD0_decayL > 9999999.) continue;  ///could be strict
        //////if(sdcaDaughters > 10.0) continue;
        /*if(sk_nSigmaTPC < -2.5 || sk_nSigmaTPC > 2.5) continue; 
        if(spi1_nSigmaTOF < -2.5 || spi1_nSigmaTOF > 2.5) continue;
        if (sk_nSigmaTOF < nSigmaTOFKaonTrack_lower(sk_p) || sk_nSigmaTOF > nSigmaTOFKaonTrack_upper(sk_p)) continue;
        */
        /////if(sdcaD0ToPv > 99.9) continue;  ///could be strict
        //////if(TMath::Abs(scosTheta) > 0.9) continue;  ////come back here
        //////if(TMath::Abs(sD_eta) > 1.0) continue;
        
        if(sD_cosThetaStar >= 0.7) continue;   ////come back here
        ////if(spi1_dca > 0.62) continue;
        ////if(sk_dca < 0.38) continue;
        //////if(striplet_pairdiff>3.5 || striplet_pairdiff<0.1) continue;
        //////cout << "after cosThetaStarCutD0" << endl;

        
        if(spi1_charge == 1 && sk_charge == -1 && spi2_charge == 1) {
            ////if(((striplet_pairdiff>1.75) && (striplet_pairdiff<1.85)) || ((striplet_pairdiff>1.9) && (striplet_pairdiff<2.0))){
            /////if(((striplet_pairdiff>1.78) && (striplet_pairdiff<1.83))){
               /*dcaDaughters_rot->Fill(sdcaDaughters);
               D_theta_rot->Fill(sD0_theta);
               cosTheta_rot->Fill(scosTheta);
               D_decayL_rot->Fill(sD0_decayL);
               dcaD0ToPv_rot->Fill(sdcaD0ToPv);
               D_cosThetaStar_rot->Fill(sD_cosThetaStar);
               Dstar_pt_rot->Fill(sDstar_pt);
               triplet_pairdiff_rot->Fill(striplet_pairdiff);
               Dstar_rapidity_rot->Fill(sDstar_rapidity);
               D_phi_rot->Fill(sD_phi);
               D_eta_rot->Fill(sD_eta);
               pi_DCA_rot->Fill(spi1_dca);
               k_DCA_rot->Fill(sk_dca);
               D_y_D_eta_rot->Fill(sDstar_rapidity, sD_eta);*/
            /////}
            /////if((bD_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(bD0_theta)) > cosThetapointD0cut)) {
          ////////if((bD_cosThetaStar < cosThetaStarCutD0)){
            /////cout << "in charge loop" << endl;
            DstarPlusSideBand->Fill(sDstar_pt,striplet_pairdiff);
            ///////}
        }
        if(spi1_charge == -1 && sk_charge == 1 && spi2_charge == -1) {
             ////if(((striplet_pairdiff>1.75) && (striplet_pairdiff<1.85)) || ((striplet_pairdiff>1.9) && (striplet_pairdiff<2.0))){
            ////if(((striplet_pairdiff>1.78) && (striplet_pairdiff<1.83))){
               /*dcaDaughters_rot->Fill(sdcaDaughters);
               D_theta_rot->Fill(sD0_theta);
               cosTheta_rot->Fill(scosTheta);
               D_decayL_rot->Fill(sD0_decayL);
               dcaD0ToPv_rot->Fill(sdcaD0ToPv);
               D_cosThetaStar_rot->Fill(sD_cosThetaStar);
               Dstar_pt_rot->Fill(sDstar_pt);
               triplet_pairdiff_rot->Fill(striplet_pairdiff);
               Dstar_rapidity_rot->Fill(sDstar_rapidity);
               D_phi_rot->Fill(sD_phi);
               D_eta_rot->Fill(sD_eta);
               pi_DCA_rot->Fill(spi1_dca);
               k_DCA_rot->Fill(sk_dca);
               D_y_D_eta_rot->Fill(sDstar_rapidity, sD_eta);*/
            /////}
            /////if((bD_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(bD0_theta)) > cosThetapointD0cut)) {
           ///////if((bD_cosThetaStar < cosThetaStarCutD0)){
            //////cout << "in charge loop" << endl;
            DstarMinusSideBand->Fill(sDstar_pt,striplet_pairdiff);
            ///////}
        }
   
    }
  
    
    
    /*numberEntr = ntp_ME -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<< "Mixed Event D0 "<<i<<endl;}
        ntp_ME -> GetEntry(i);
        ///if (Dstar_pt<=0) continue;
        //////if ((sDstar_pt<1.1) || (sDstar_pt>2.1))  continue;
        if((TMath::Abs(mDstar_rapidity)>1.0)) continue;
        
        //////if(sD0_decayL < 1.0) continue;  ///could be strict
        //////if(sdcaDaughters > 0.48) continue;
        /////if(mdcaD0ToPv > 99.9) continue;  ///could be strict
        /*if(mk_nSigmaTPC < -2.5 || mk_nSigmaTPC > 2.5) continue;  
        if(mpi1_nSigmaTOF < -2.5 || mpi1_nSigmaTOF > 2.5) continue;
        if (mk_nSigmaTOF < nSigmaTOFKaonTrack_lower(mk_p) || mk_nSigmaTOF > nSigmaTOFKaonTrack_upper(mk_p)) continue;
        
        //////if(TMath::Abs(mcosTheta) > 0.90) continue;  ////come back here
        //////if(TMath::Abs(mD_eta) > 1.0) continue;
        
        ////////if(mD_cosThetaStar > 0.77) continue;   ////come back here
        ////if(spi1_dca > 0.62) continue;
        ////if(sk_dca < 0.38) continue;
        if(mpi1_charge == 1 && mk_charge == -1 && mpi2_charge == 1) {
            ////if(((mtriplet_pairdiff>1.75) && (mtriplet_pairdiff<1.85)) || ((mtriplet_pairdiff>1.9) && (mtriplet_pairdiff<2.0))){
            /////if(((mtriplet_pairdiff>1.78) && (mtriplet_pairdiff<1.83))){
               dcaDaughters_mixed->Fill(mdcaDaughters);
               D_theta_mixed->Fill(mD0_theta);
               cosTheta_mixed->Fill(mcosTheta);
               D_decayL_mixed->Fill(mD0_decayL);
               dcaD0ToPv_mixed->Fill(mdcaD0ToPv);
               D_cosThetaStar_mixed->Fill(mD_cosThetaStar);
               Dstar_pt_mixed->Fill(mDstar_pt);
               triplet_pairdiff_mixed->Fill(mtriplet_pairdiff);
               Dstar_rapidity_mixed->Fill(mDstar_rapidity);
               D_phi_mixed->Fill(mD_phi);
               D_eta_mixed->Fill(mD_eta);
               pi_DCA_mixed->Fill(mpi1_dca);
               k_DCA_mixed->Fill(mk_dca);
               D_y_D_eta_mixed->Fill(mDstar_rapidity, mD_eta);
            ////}
            ///////if(mDstar_pt>pTforCos) {
            ///////if((mD_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(mD0_theta)) > cosThetapointD0cut)) {
            /////if((mD_cosThetaStar < cosThetaStarCutD0)){
                    MixedDstarPlusRightSign->Fill(mDstar_pt,mtriplet_pairdiff);
            /////}               
            ///////}
            /*else{
					MixedDstarPlusRightSign->Fill(mDstar_pt,mtriplet_pairdiff);	
				}
                
        }
        if(mpi1_charge == -1 && mk_charge == 1 && mpi2_charge == -1) {
                  ////if(((mtriplet_pairdiff>1.75) && (mtriplet_pairdiff<1.85)) || ((mtriplet_pairdiff>1.9) && (mtriplet_pairdiff<2.0))){
            /////if(((mtriplet_pairdiff>1.78) && (mtriplet_pairdiff<1.83))){
               dcaDaughters_mixed->Fill(mdcaDaughters);
               D_theta_mixed->Fill(mD0_theta);
               cosTheta_mixed->Fill(mcosTheta);
               D_decayL_mixed->Fill(mD0_decayL);
               dcaD0ToPv_mixed->Fill(mdcaD0ToPv);
               D_cosThetaStar_mixed->Fill(mD_cosThetaStar);
               Dstar_pt_mixed->Fill(mDstar_pt);
               triplet_pairdiff_mixed->Fill(mtriplet_pairdiff);
               Dstar_rapidity_mixed->Fill(mDstar_rapidity);
               D_phi_mixed->Fill(mD_phi);
               D_eta_mixed->Fill(mD_eta);
               pi_DCA_mixed->Fill(mpi1_dca);
               k_DCA_mixed->Fill(mk_dca);
               D_y_D_eta_mixed->Fill(mDstar_rapidity, mD_eta);
            ////}
            ///////if(mDstar_pt>pTforCos) {
            ///////if((mD_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(mD0_theta)) > cosThetapointD0cut)) {
            ///////if((mD_cosThetaStar < cosThetaStarCutD0)){
                    MixedDstarMinusRightSign->Fill(mDstar_pt,mtriplet_pairdiff);
            //////}               
            ///////}
            /*else{
					MixedDstarMinusRightSign->Fill(mDstar_pt,mtriplet_pairdiff);	
				}
                
        }
   
    }*/
    
    


    numberEntr = ntp_WrongSign -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<< "WrongSign Dstar "<<i<<endl;}
        ntp_WrongSign -> GetEntry(i);
        ///if (Dstar_pt<=0) continue;
        //////if ((sDstar_pt<1.1) || (sDstar_pt>2.1))  continue;
        if((TMath::Abs(wDstar_rapidity)>=1.0)) continue;

        if(wk_nSigmaTPC <= -2.5 || wk_nSigmaTPC >= 3.0) continue; 
        if(wk_nSigmaTOF <= nSigmaTOFKaonTrack_lower(wk_p) || wk_nSigmaTOF >= nSigmaTOFKaonTrack_upper(wk_p)) continue;

        if(wpi1_nSigmaTPC < -3.0 || wpi1_nSigmaTPC > 3.0) continue;
        if(wpi1_nSigmaTOF < -2.5 || wpi1_nSigmaTOF > 2.5) continue;

        if(wpi2_nSigmaTPC < -3.0 || wpi2_nSigmaTPC > 3.0) continue;
        if(wpi2_nSigmaTOF < -2.5 || wpi2_nSigmaTOF > 2.5) continue;
        
        //////if(sD0_decayL < 1.0) continue;  ///could be strict
        //////if(sdcaDaughters > 0.48) continue;
        /////if(wdcaD0ToPv > 99.9) continue;  ///could be strict
        /*if(wk_nSigmaTPC < -2.5 || wk_nSigmaTPC > 2.5) continue;
        if(wpi1_nSigmaTOF < -2.5 || wpi1_nSigmaTOF > 2.5) continue;
        if (wk_nSigmaTOF < nSigmaTOFKaonTrack_lower(wk_p) || wk_nSigmaTOF > nSigmaTOFKaonTrack_upper(wk_p)) continue; 
        */
        ///////if(TMath::Abs(wcosTheta) > 0.90) continue;  ////come back here
        ///////if(TMath::Abs(wD_eta) > 1.0) continue;
        
        if(wD_cosThetaStar >= 0.7) continue;  ////come back here
        ////if(spi1_dca > 0.62) continue;
        ////if(sk_dca < 0.38) continue;

        if(wpi1_charge == 1 && wk_charge == -1 && wpi2_charge == -1) {
               ////if(((wtriplet_pairdiff>1.75) && (wtriplet_pairdiff<1.85)) || ((wtriplet_pairdiff>1.9) && (wtriplet_pairdiff<2.0))){
            /////if(((wtriplet_pairdiff>1.78) && (wtriplet_pairdiff<1.83))){
              /*dcaDaughters_like->Fill(wdcaDaughters);
               D_theta_like->Fill(wD0_theta);
               cosTheta_like->Fill(wcosTheta);
               D_decayL_like->Fill(wD0_decayL);
               dcaD0ToPv_like->Fill(wdcaD0ToPv);
               D_cosThetaStar_like->Fill(wD_cosThetaStar);
               Dstar_pt_like->Fill(wDstar_pt);
               triplet_pairdiff_like->Fill(wtriplet_pairdiff);
               Dstar_rapidity_like->Fill(wDstar_rapidity);
               D_phi_like->Fill(wD_phi);
               D_eta_like->Fill(wD_eta);
               pi_DCA_like->Fill(wpi1_dca);
               k_DCA_like->Fill(wk_dca);
               D_y_D_eta_like->Fill(wDstar_rapidity, wD_eta);*/
           ///// }
            ///////if(wDstar_pt>pTforCos) {
            ///////if((wD_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(wD0_theta)) > cosThetapointD0cut)) {
            /////if((wD_cosThetaStar < cosThetaStarCutD0)){
                    DstarPlusWrongSign->Fill(wDstar_pt,wtriplet_pairdiff);
            /////}               
            ///////}
            /*else{
					DstarPlusWrongSign->Fill(wDstar_pt,wtriplet_pairdiff);	
				}*/
                
        }
        if(wpi1_charge == -1 && wk_charge == 1 && wpi2_charge == 1) {

            ////if(((wtriplet_pairdiff>1.75) && (wtriplet_pairdiff<1.85)) || ((wtriplet_pairdiff>1.9) && (wtriplet_pairdiff<2.0))){
            ////if(((wtriplet_pairdiff>1.78) && (wtriplet_pairdiff<1.83))){
              /* dcaDaughters_like->Fill(wdcaDaughters);
               D_theta_like->Fill(wD0_theta);
               cosTheta_like->Fill(wcosTheta);
               D_decayL_like->Fill(wD0_decayL);
               dcaD0ToPv_like->Fill(wdcaD0ToPv);
               D_cosThetaStar_like->Fill(wD_cosThetaStar);
               Dstar_pt_like->Fill(wDstar_pt);
               triplet_pairdiff_like->Fill(wtriplet_pairdiff);
               Dstar_rapidity_like->Fill(wDstar_rapidity);
               D_phi_like->Fill(wD_phi);
               D_eta_like->Fill(wD_eta);
               pi_DCA_like->Fill(wpi1_dca);
               k_DCA_like->Fill(wpi1_dca);
               D_y_D_eta_like->Fill(wDstar_rapidity, wD_eta);*/
            /////}
            ///////if(wDstar_pt>pTforCos) {
            ///////if((wD_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(wD0_theta)) > cosThetapointD0cut)) {
            ///////if((wD_cosThetaStar < cosThetaStarCutD0)){
                    DstarMinusWrongSign->Fill(wDstar_pt,wtriplet_pairdiff);
            //////}               
            ///////}
            /*else{
					DstarMinusWrongSign->Fill(wDstar_pt,wtriplet_pairdiff);	
				}*/
                
        }
   
    }
    



   
    TFile* dataRes = new TFile("pp500_prereq_Dstarmidstep_costhetastar0.7.root","RECREATE");


    DstarPlusRightSign->Write();
    DstarMinusRightSign->Write();
    DstarPlusSideBand->Write();
    DstarMinusSideBand->Write();
    MixedDstarPlusRightSign->Write();
    MixedDstarMinusRightSign->Write();
    DstarPlusWrongSign->Write();
    DstarMinusWrongSign->Write();
    /*hRapidityD->Write();
    hDptPiPt->Write();
    hDptKPt->Write();
    hDptPiP->Write();
	hDptKP->Write();
    hPtPionVsKaon->Write();
    hNHitsKvsD0pt->Write();
	hNHitsPivsD0pt->Write();
    hDstarRapidity->Write();
    hSoftPionPlus->Write();
	hSoftPionMinus->Write();
    hD0PtVsSoftPionPt->Write();
    hSideBand->Write();
    hWrongSign->Write();
    hD0SoftPionCandidatesCorr->Write();*/

    /*dcaDaughters_unlike->Scale(1.0/dcaDaughters_unlike->GetEntries());
    D_theta_unlike->Scale(1.0/D_theta_unlike->GetEntries());
    cosTheta_unlike->Scale(1.0/cosTheta_unlike->GetEntries());
    D_decayL_unlike->Scale(1.0/D_decayL_unlike->GetEntries());
    dcaD0ToPv_unlike->Scale(1.0/dcaD0ToPv_unlike->GetEntries());
    D_cosThetaStar_unlike->Scale(1.0/D_cosThetaStar_unlike->GetEntries());
    Dstar_pt_unlike->Scale(1.0/Dstar_pt_unlike->GetEntries());
    triplet_pairdiff_unlike->Scale(1.0/triplet_pairdiff_unlike->GetEntries());
    Dstar_rapidity_unlike->Scale(1.0/Dstar_rapidity_unlike->GetEntries());
    D_phi_unlike->Scale(1.0/D_phi_unlike->GetEntries());
    D_eta_unlike->Scale(1.0/D_eta_unlike->GetEntries());*/

    /*dcaDaughters_unlike->Write();
    D_theta_unlike->Write();
    cosTheta_unlike->Write();
    D_decayL_unlike->Write();
    dcaD0ToPv_unlike->Write();
    D_cosThetaStar_unlike->Write();
    Dstar_pt_unlike->Write();
    triplet_pairdiff_unlike->Write();
    Dstar_rapidity_unlike->Write();
    D_phi_unlike->Write();
    D_eta_unlike->Write();
    pi_DCA_unlike->Write();
    k_DCA_unlike->Write();
    D_y_D_eta_unlike->Write();*/


    /*dcaDaughters_like->Scale(1.0/dcaDaughters_like->GetEntries());
    D_theta_like->Scale(1.0/D_theta_like->GetEntries());
    cosTheta_like->Scale(1.0/cosTheta_like->GetEntries());
    D_decayL_like->Scale(1.0/D_decayL_like->GetEntries());
    dcaD0ToPv_like->Scale(1.0/dcaD0ToPv_like->GetEntries());
    D_cosThetaStar_like->Scale(1.0/D_cosThetaStar_like->GetEntries());
    Dstar_pt_like->Scale(1.0/Dstar_pt_like->GetEntries());
    triplet_pairdiff_like->Scale(1.0/triplet_pairdiff_like->GetEntries());
    Dstar_rapidity_like->Scale(1.0/Dstar_rapidity_like->GetEntries());
    D_phi_like->Scale(1.0/D_phi_like->GetEntries());
    D_eta_like->Scale(1.0/D_eta_like->GetEntries());*/

    /*dcaDaughters_like->Write();
    D_theta_like->Write();
    cosTheta_like->Write();
    D_decayL_like->Write();
    dcaD0ToPv_like->Write();
    D_cosThetaStar_like->Write();
    Dstar_pt_like->Write();
    triplet_pairdiff_like->Write();
    Dstar_rapidity_like->Write();
    D_phi_like->Write();
    D_eta_like->Write();
    pi_DCA_like->Write();
    k_DCA_like->Write();
    D_y_D_eta_like->Write();*/

    /*dcaDaughters_rot->Scale(1.0/dcaDaughters_rot->GetEntries());
    D_theta_rot->Scale(1.0/D_theta_rot->GetEntries());
    cosTheta_rot->Scale(1.0/cosTheta_rot->GetEntries());
    D_decayL_rot->Scale(1.0/D_decayL_rot->GetEntries());
    dcaD0ToPv_rot->Scale(1.0/dcaD0ToPv_rot->GetEntries());
    D_cosThetaStar_rot->Scale(1.0/D_cosThetaStar_rot->GetEntries());
    Dstar_pt_rot->Scale(1.0/Dstar_pt_rot->GetEntries());
    triplet_pairdiff_rot->Scale(1.0/triplet_pairdiff_rot->GetEntries());
    Dstar_rapidity_rot->Scale(1.0/Dstar_rapidity_rot->GetEntries());
    D_phi_rot->Scale(1.0/D_phi_rot->GetEntries());
    D_eta_rot->Scale(1.0/D_eta_rot->GetEntries());*/

    /*dcaDaughters_rot->Write();
    D_theta_rot->Write();
    cosTheta_rot->Write();
    D_decayL_rot->Write();
    dcaD0ToPv_rot->Write();
    D_cosThetaStar_rot->Write();
    Dstar_pt_rot->Write();
    triplet_pairdiff_rot->Write();
    Dstar_rapidity_rot->Write();
    D_phi_rot->Write();
    D_eta_rot->Write();
    pi_DCA_rot->Write();
    k_DCA_rot->Write();
    D_y_D_eta_rot->Write();*/


    /*dcaDaughters_mixed->Scale(1.0/dcaDaughters_mixed->GetEntries());
    D_theta_mixed->Scale(1.0/D_theta_mixed->GetEntries());
    cosTheta_mixed->Scale(1.0/cosTheta_mixed->GetEntries());
    D_decayL_mixed->Scale(1.0/D_decayL_mixed->GetEntries());
    dcaD0ToPv_mixed->Scale(1.0/dcaD0ToPv_mixed->GetEntries());
    D_cosThetaStar_mixed->Scale(1.0/D_cosThetaStar_mixed->GetEntries());
    Dstar_pt_mixed->Scale(1.0/Dstar_pt_mixed->GetEntries());
    triplet_pairdiff_mixed->Scale(1.0/triplet_pairdiff_mixed->GetEntries());
    Dstar_rapidity_mixed->Scale(1.0/Dstar_rapidity_mixed->GetEntries());
    D_phi_mixed->Scale(1.0/D_phi_mixed->GetEntries());
    D_eta_mixed->Scale(1.0/D_eta_mixed->GetEntries());*/

    /*dcaDaughters_mixed->Write();
    D_theta_mixed->Write();
    cosTheta_mixed->Write();
    D_decayL_mixed->Write();
    dcaD0ToPv_mixed->Write();
    D_cosThetaStar_mixed->Write();
    Dstar_pt_mixed->Write();
    triplet_pairdiff_mixed->Write();
    Dstar_rapidity_mixed->Write();
    D_phi_mixed->Write();
    D_eta_mixed->Write();
    pi_DCA_mixed->Write();
    k_DCA_mixed->Write();
    D_y_D_eta_mixed->Write();*/


    /*TCanvas *c1 = new TCanvas("c1");
    c1->cd();
    dcaDaughters_unlike->SetMarkerStyle(7);
    dcaDaughters_unlike->SetMarkerColor(1);
    dcaDaughters_unlike->SetLineColor(1);
	dcaDaughters_unlike->SetLineWidth(2);
	dcaDaughters_unlike->GetYaxis()->SetTitleFont(42);
	dcaDaughters_unlike->GetYaxis()->SetLabelFont(42);
	dcaDaughters_unlike->Draw("");
    dcaDaughters_mixed->SetMarkerStyle(7);
    dcaDaughters_mixed->SetMarkerColor(2);
	dcaDaughters_mixed->SetLineColor(2);
	dcaDaughters_mixed->SetLineWidth(2);
    dcaDaughters_mixed->Draw("Csame");
	dcaDaughters_like->SetMarkerStyle(7);
	dcaDaughters_like->SetMarkerColor(4);
	dcaDaughters_like->SetLineColor(4);
	dcaDaughters_like->SetLineWidth(2);
	dcaDaughters_like->Draw("hist C same");
	dcaDaughters_rot->SetMarkerStyle(7);
	dcaDaughters_rot->SetMarkerColor(8);
	dcaDaughters_rot->SetLineColor(8);
	dcaDaughters_rot->SetLineWidth(2);
	dcaDaughters_rot->Draw("hist C same");

    TCanvas *c2 = new TCanvas("c2");
    c2->cd();
    D_theta_unlike->SetMarkerStyle(7);
    D_theta_unlike->SetMarkerColor(1);
    D_theta_unlike->SetLineColor(1);
	D_theta_unlike->SetLineWidth(2);
	D_theta_unlike->GetYaxis()->SetTitleFont(42);
	D_theta_unlike->GetYaxis()->SetLabelFont(42);
    D_theta_unlike->GetXaxis()->SetTitleFont(42);
	D_theta_unlike->GetXaxis()->SetLabelFont(42);
	D_theta_unlike->Draw("");
    D_theta_mixed->SetMarkerStyle(7);
    D_theta_mixed->SetMarkerColor(2);
	D_theta_mixed->SetLineColor(2);
	D_theta_mixed->SetLineWidth(2);
    D_theta_mixed->Draw("Csame");
	D_theta_like->SetMarkerStyle(7);
	D_theta_like->SetMarkerColor(4);
	D_theta_like->SetLineColor(4);
	D_theta_like->SetLineWidth(2);
	D_theta_like->Draw("hist C same");
	D_theta_rot->SetMarkerStyle(7);
	D_theta_rot->SetMarkerColor(8);
	D_theta_rot->SetLineColor(8);
	D_theta_rot->SetLineWidth(2);
	D_theta_rot->Draw("hist C same");

    TCanvas *c3 = new TCanvas("c3");
    c3->cd();
    cosTheta_unlike->SetMarkerStyle(7);
    cosTheta_unlike->SetMarkerColor(1);
    cosTheta_unlike->SetLineColor(1);
	cosTheta_unlike->SetLineWidth(2);
	cosTheta_unlike->GetYaxis()->SetTitleFont(42);
	cosTheta_unlike->GetYaxis()->SetLabelFont(42);
    cosTheta_unlike->GetXaxis()->SetTitleFont(42);
	cosTheta_unlike->GetXaxis()->SetLabelFont(42);
	cosTheta_unlike->Draw("");
    cosTheta_mixed->SetMarkerStyle(7);
    cosTheta_mixed->SetMarkerColor(2);
	cosTheta_mixed->SetLineColor(2);
	cosTheta_mixed->SetLineWidth(2);
    cosTheta_mixed->Draw("Csame");
	cosTheta_like->SetMarkerStyle(7);
	cosTheta_like->SetMarkerColor(4);
	cosTheta_like->SetLineColor(4);
	cosTheta_like->SetLineWidth(2);
	cosTheta_like->Draw("hist C same");
	cosTheta_rot->SetMarkerStyle(7);
	cosTheta_rot->SetMarkerColor(8);
	cosTheta_rot->SetLineColor(8);
	cosTheta_rot->SetLineWidth(2);
	cosTheta_rot->Draw("hist C same");

    TCanvas *c4 = new TCanvas("c4");
    c4->cd();
    D_decayL_unlike->SetMarkerStyle(7);
    D_decayL_unlike->SetMarkerColor(1);
    D_decayL_unlike->SetLineColor(1);
	D_decayL_unlike->SetLineWidth(2);
	D_decayL_unlike->GetYaxis()->SetTitleFont(42);
	D_decayL_unlike->GetYaxis()->SetLabelFont(42);
    D_decayL_unlike->GetXaxis()->SetTitleFont(42);
	D_decayL_unlike->GetXaxis()->SetLabelFont(42);
	D_decayL_unlike->Draw("");
    D_decayL_mixed->SetMarkerStyle(7);
    D_decayL_mixed->SetMarkerColor(2);
	D_decayL_mixed->SetLineColor(2);
	D_decayL_mixed->SetLineWidth(2);
    D_decayL_mixed->Draw("Csame");
	D_decayL_like->SetMarkerStyle(7);
	D_decayL_like->SetMarkerColor(4);
	D_decayL_like->SetLineColor(4);
	D_decayL_like->SetLineWidth(2);
	D_decayL_like->Draw("hist C same");
	D_decayL_rot->SetMarkerStyle(7);
	D_decayL_rot->SetMarkerColor(8);
	D_decayL_rot->SetLineColor(8);
	D_decayL_rot->SetLineWidth(2);
	D_decayL_rot->Draw("hist C same");

    TCanvas *c5 = new TCanvas("c5");
    c5->cd();
    dcaD0ToPv_unlike->SetMarkerStyle(7);
    dcaD0ToPv_unlike->SetMarkerColor(1);
    dcaD0ToPv_unlike->SetLineColor(1);
	dcaD0ToPv_unlike->SetLineWidth(2);
	dcaD0ToPv_unlike->GetYaxis()->SetTitleFont(42);
	dcaD0ToPv_unlike->GetYaxis()->SetLabelFont(42);
    dcaD0ToPv_unlike->GetXaxis()->SetTitleFont(42);
	dcaD0ToPv_unlike->GetXaxis()->SetLabelFont(42);
	dcaD0ToPv_unlike->Draw("");
    dcaD0ToPv_mixed->SetMarkerStyle(7);
    dcaD0ToPv_mixed->SetMarkerColor(2);
	dcaD0ToPv_mixed->SetLineColor(2);
	dcaD0ToPv_mixed->SetLineWidth(2);
    dcaD0ToPv_mixed->Draw("Csame");
	dcaD0ToPv_like->SetMarkerStyle(7);
	dcaD0ToPv_like->SetMarkerColor(4);
	dcaD0ToPv_like->SetLineColor(4);
	dcaD0ToPv_like->SetLineWidth(2);
	dcaD0ToPv_like->Draw("hist C same");
	dcaD0ToPv_rot->SetMarkerStyle(7);
	dcaD0ToPv_rot->SetMarkerColor(8);
	dcaD0ToPv_rot->SetLineColor(8);
	dcaD0ToPv_rot->SetLineWidth(2);
	dcaD0ToPv_rot->Draw("hist C same");

    TCanvas *c6 = new TCanvas("c6");
    c6->cd();
    D_cosThetaStar_unlike->SetMarkerStyle(7);
    D_cosThetaStar_unlike->SetMarkerColor(1);
    D_cosThetaStar_unlike->SetLineColor(1);
	D_cosThetaStar_unlike->SetLineWidth(2);
	D_cosThetaStar_unlike->GetYaxis()->SetTitleFont(42);
	D_cosThetaStar_unlike->GetYaxis()->SetLabelFont(42);
    D_cosThetaStar_unlike->GetXaxis()->SetTitleFont(42);
	D_cosThetaStar_unlike->GetXaxis()->SetLabelFont(42);
	D_cosThetaStar_unlike->Draw("");
    D_cosThetaStar_mixed->SetMarkerStyle(7);
    D_cosThetaStar_mixed->SetMarkerColor(2);
	D_cosThetaStar_mixed->SetLineColor(2);
	D_cosThetaStar_mixed->SetLineWidth(2);
    D_cosThetaStar_mixed->Draw("Csame");
	D_cosThetaStar_like->SetMarkerStyle(7);
	D_cosThetaStar_like->SetMarkerColor(4);
	D_cosThetaStar_like->SetLineColor(4);
	D_cosThetaStar_like->SetLineWidth(2);
	D_cosThetaStar_like->Draw("hist C same");
	D_cosThetaStar_rot->SetMarkerStyle(7);
	D_cosThetaStar_rot->SetMarkerColor(8);
	D_cosThetaStar_rot->SetLineColor(8);
	D_cosThetaStar_rot->SetLineWidth(2);
	D_cosThetaStar_rot->Draw("hist C same");

    TCanvas *c7 = new TCanvas("c7");
    c7->cd();
    Dstar_pt_unlike->SetMarkerStyle(7);
    Dstar_pt_unlike->SetMarkerColor(1);
    Dstar_pt_unlike->SetLineColor(1);
	Dstar_pt_unlike->SetLineWidth(2);
	Dstar_pt_unlike->GetYaxis()->SetTitleFont(42);
	Dstar_pt_unlike->GetYaxis()->SetLabelFont(42);
    Dstar_pt_unlike->GetXaxis()->SetTitleFont(42);
	Dstar_pt_unlike->GetXaxis()->SetLabelFont(42);
	Dstar_pt_unlike->Draw("");
    Dstar_pt_mixed->SetMarkerStyle(7);
    Dstar_pt_mixed->SetMarkerColor(2);
	Dstar_pt_mixed->SetLineColor(2);
	Dstar_pt_mixed->SetLineWidth(2);
    Dstar_pt_mixed->Draw("Csame");
	Dstar_pt_like->SetMarkerStyle(7);
	Dstar_pt_like->SetMarkerColor(4);
	Dstar_pt_like->SetLineColor(4);
	Dstar_pt_like->SetLineWidth(2);
	Dstar_pt_like->Draw("hist C same");
	Dstar_pt_rot->SetMarkerStyle(7);
	Dstar_pt_rot->SetMarkerColor(8);
	Dstar_pt_rot->SetLineColor(8);
	Dstar_pt_rot->SetLineWidth(2);
	Dstar_pt_rot->Draw("hist C same");

    TCanvas *c8 = new TCanvas("c8");
    c8->cd();
    triplet_pairdiff_unlike->SetMarkerStyle(7);
    triplet_pairdiff_unlike->SetMarkerColor(1);
    triplet_pairdiff_unlike->SetLineColor(1);
	triplet_pairdiff_unlike->SetLineWidth(2);
	triplet_pairdiff_unlike->GetYaxis()->SetTitleFont(42);
	triplet_pairdiff_unlike->GetYaxis()->SetLabelFont(42);
    triplet_pairdiff_unlike->GetXaxis()->SetTitleFont(42);
	triplet_pairdiff_unlike->GetXaxis()->SetLabelFont(42);
	triplet_pairdiff_unlike->Draw("");
    triplet_pairdiff_mixed->SetMarkerStyle(7);
    triplet_pairdiff_mixed->SetMarkerColor(2);
	triplet_pairdiff_mixed->SetLineColor(2);
	triplet_pairdiff_mixed->SetLineWidth(2);
    triplet_pairdiff_mixed->Draw("Csame");
	triplet_pairdiff_like->SetMarkerStyle(7);
	triplet_pairdiff_like->SetMarkerColor(4);
	triplet_pairdiff_like->SetLineColor(4);
	triplet_pairdiff_like->SetLineWidth(2);
	triplet_pairdiff_like->Draw("hist C same");
	triplet_pairdiff_rot->SetMarkerStyle(7);
	triplet_pairdiff_rot->SetMarkerColor(8);
	triplet_pairdiff_rot->SetLineColor(8);
	triplet_pairdiff_rot->SetLineWidth(2);
	triplet_pairdiff_rot->Draw("hist C same");


    TCanvas *c9 = new TCanvas("c9");
    c9->cd();
    Dstar_rapidity_unlike->SetMarkerStyle(7);
    Dstar_rapidity_unlike->SetMarkerColor(1);
    Dstar_rapidity_unlike->SetLineColor(1);
	Dstar_rapidity_unlike->SetLineWidth(2);
	Dstar_rapidity_unlike->GetYaxis()->SetTitleFont(42);
	Dstar_rapidity_unlike->GetYaxis()->SetLabelFont(42);
    Dstar_rapidity_unlike->GetXaxis()->SetTitleFont(42);
	Dstar_rapidity_unlike->GetXaxis()->SetLabelFont(42);
	Dstar_rapidity_unlike->Draw("");
    Dstar_rapidity_mixed->SetMarkerStyle(7);
    Dstar_rapidity_mixed->SetMarkerColor(2);
	Dstar_rapidity_mixed->SetLineColor(2);
	Dstar_rapidity_mixed->SetLineWidth(2);
    Dstar_rapidity_mixed->Draw("Csame");
	Dstar_rapidity_like->SetMarkerStyle(7);
	Dstar_rapidity_like->SetMarkerColor(4);
	Dstar_rapidity_like->SetLineColor(4);
	Dstar_rapidity_like->SetLineWidth(2);
	Dstar_rapidity_like->Draw("hist C same");
	Dstar_rapidity_rot->SetMarkerStyle(7);
	Dstar_rapidity_rot->SetMarkerColor(8);
	Dstar_rapidity_rot->SetLineColor(8);
	Dstar_rapidity_rot->SetLineWidth(2);
	Dstar_rapidity_rot->Draw("hist C same");

    TCanvas *c10 = new TCanvas("c10");
    c10->cd();
    D_phi_unlike->SetMarkerStyle(7);
    D_phi_unlike->SetMarkerColor(1);
    D_phi_unlike->SetLineColor(1);
	D_phi_unlike->SetLineWidth(2);
	D_phi_unlike->GetYaxis()->SetTitleFont(42);
	D_phi_unlike->GetYaxis()->SetLabelFont(42);
    D_phi_unlike->GetXaxis()->SetTitleFont(42);
	D_phi_unlike->GetXaxis()->SetLabelFont(42);
	D_phi_unlike->Draw("");
    D_phi_mixed->SetMarkerStyle(7);
    D_phi_mixed->SetMarkerColor(2);
	D_phi_mixed->SetLineColor(2);
	D_phi_mixed->SetLineWidth(2);
    D_phi_mixed->Draw("Csame");
	D_phi_like->SetMarkerStyle(7);
	D_phi_like->SetMarkerColor(4);
	D_phi_like->SetLineColor(4);
	D_phi_like->SetLineWidth(2);
	D_phi_like->Draw("hist C same");
	D_phi_rot->SetMarkerStyle(7);
	D_phi_rot->SetMarkerColor(8);
	D_phi_rot->SetLineColor(8);
	D_phi_rot->SetLineWidth(2);
	D_phi_rot->Draw("hist C same");

    TCanvas *c11 = new TCanvas("c11");
    c11->cd();
    D_eta_unlike->SetMarkerStyle(7);
    D_eta_unlike->SetMarkerColor(1);
    D_eta_unlike->SetLineColor(1);
	D_eta_unlike->SetLineWidth(2);
	D_eta_unlike->GetYaxis()->SetTitleFont(42);
	D_eta_unlike->GetYaxis()->SetLabelFont(42);
    D_eta_unlike->GetXaxis()->SetTitleFont(42);
	D_eta_unlike->GetXaxis()->SetLabelFont(42);
	D_eta_unlike->Draw("");
    D_eta_mixed->SetMarkerStyle(7);
    D_eta_mixed->SetMarkerColor(2);
	D_eta_mixed->SetLineColor(2);
	D_eta_mixed->SetLineWidth(2);
    D_eta_mixed->Draw("Csame");
	D_eta_like->SetMarkerStyle(7);
	D_eta_like->SetMarkerColor(4);
	D_eta_like->SetLineColor(4);
	D_eta_like->SetLineWidth(2);
	D_eta_like->Draw("hist C same");
	D_eta_rot->SetMarkerStyle(7);
	D_eta_rot->SetMarkerColor(8);
	D_eta_rot->SetLineColor(8);
	D_eta_rot->SetLineWidth(2);
	D_eta_rot->Draw("hist C same");*/

    dataRes->Close();

    
    cout<<"Finished"<<endl;


}

