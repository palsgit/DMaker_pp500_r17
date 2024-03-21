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
	double f_nsigmaTOFKaonTrack_res =  1.35 + (0.02/(pow((x + 0.05), 3.76)));
    double f_nsigmaTOFKaonTrack_pos =  -0.01 + (0.02/(pow((x - 0.17), 1.81))); 

    /*float mSigmahigher = 3.2;
    ///if (x > 2.0) mSigmahigher = 2.4;*/
    float mSigma = 2.22;
	double kaon_higher = mSigma*f_nsigmaTOFKaonTrack_res + f_nsigmaTOFKaonTrack_pos;
    //////double kaon_higher = 3.0;
    //////////if (x < 2.0) kaon_higher = (37/9.0) - (5.0/9.0)*x;
    return kaon_higher; 
	}
Double_t nSigmaTOFKaonTrack_lower(Double_t x) {
	double f_nsigmaTOFKaonTrack_res =  1.35 + (0.02/(pow((x + 0.05), 3.76)));
    double f_nsigmaTOFKaonTrack_pos =  -0.01 + (0.02/(pow((x - 0.17), 1.81))); 
    /*float mSigma = -2.30;
    /*if (x > 1.25 && x < 1.50) mSigmalower = -1.56;
    if (x > 1.5 && x < 1.65) mSigmalower = -1.18;
    if (x > 1.65) mSigmalower = -0.79;*/
    float mSigma = -1.11;
    //////if (x < 1.0) mSigma = -2.94;
	double kaon_lower = mSigma*f_nsigmaTOFKaonTrack_res + f_nsigmaTOFKaonTrack_pos;
	 //////double kaon_lower = -1.5;
    if (x < 1.5) kaon_lower = 2.05*x - 4.575;
	return kaon_lower;
	 }


Double_t DeltaTOFKaonTrack_upper(Double_t x) { 
	double f_DeltaTOFKaonTrack_res =  0.0107886 + (0.00139274/(pow((x - 0.182004), 1.01763)));
    double f_DeltaTOFKaonTrack_pos =  -9.27541e-05 + (0.000215893/(pow((x - 0.17385), 1.67005))); 

    /*float mSigma = 3.18;
    if (x > 2.0) mSigma = 2.41;*/
    float mSigma = 2.5;
	double kaon_higher = mSigma*f_DeltaTOFKaonTrack_res + f_DeltaTOFKaonTrack_pos;
    return kaon_higher; 
	}
Double_t DeltaTOFKaonTrack_lower(Double_t x) {
	double f_DeltaTOFKaonTrack_res =  0.0107886 + (0.00139274/(pow((x - 0.183716), 1.01763)));
    double f_DeltaTOFKaonTrack_pos =  -9.27541e-05 + (0.000215893/(pow((x -0.17385), 1.67005))); 
    /*float mSigma = -2.31;
    if (x > 1.25 && x < 1.50) mSigma = -1.56;
    if (x > 1.50 && x < 1.65) mSigma = -1.18;
    if (x > 1.65) mSigma = -0.8;*/
    float mSigma = -1.0;
	 double kaon_lower = mSigma*f_DeltaTOFKaonTrack_res + f_DeltaTOFKaonTrack_pos;
	 return kaon_lower; 
	 }

    Double_t nSigmaTPCPionTrack_upper(Double_t x) { 
	double f_nSigmaTPCPionTrack_res =  0.9548 + (0.0002/(pow((x + 0.1046), 3.9041)));
    double f_nSigmaTPCPionTrack_pos =  0.0988 + (-0.0012/(pow((x - 0.0706), 3.1974))); 

    float mSigma = 2.92;
    //////if (x > 2.0) mSigma = 2.41;
	double pion_higher = mSigma*f_nSigmaTPCPionTrack_res + f_nSigmaTPCPionTrack_pos;
    return pion_higher; 
	}
     Double_t nSigmaTPCPionTrack_lower(Double_t x) {
	double f_nSigmaTPCPionTrack_res =  0.9548 + (0.0002/(pow((x + 0.1046), 3.9041)));
    double f_nSigmaTPCPionTrack_pos =  0.0988 + (-0.0012/(pow((x - 0.0706), 3.1974))); 
    float mSigma = -2.70;
    /*if (x > 1.25 && x < 1.50) mSigma = -1.56;
    if (x > 1.50 && x < 1.65) mSigma = -1.18;
    if (x > 1.65) mSigma = -0.8;*/
	 double pion_lower = mSigma*f_nSigmaTPCPionTrack_res + f_nSigmaTPCPionTrack_pos;
	 return pion_lower; 
	 }

    Double_t nSigmaTOFPionTrack_upper(Double_t x) { 
	double f_nSigmaTOFPionTrack_res =  0.907 + (0.002/(pow((x - 0.164), 1.990)));
    double f_nSigmaTOFPionTrack_pos =  0.016 + (0.004/(pow((x - 0.064), 2.828))); 

    float mSigma = 2.74;

    double Pion_higher = mSigma*f_nSigmaTOFPionTrack_res + f_nSigmaTOFPionTrack_pos;
    /*double Pion_higher = 2.5;
    if (x < 2.0) Pion_higher = (-5.0/18.0)*x + (55.0/18.0);*/
    return Pion_higher; 
	}
Double_t nSigmaTOFPionTrack_lower(Double_t x) {
	double f_nSigmaTOFPionTrack_res =  0.907 + (0.002/(pow((x - 0.164), 1.990)));
    double f_nSigmaTOFPionTrack_pos = 0.016 + (0.004/(pow((x - 0.064), 2.828))); 
    float mSigma = -4.43;
    /////float mSigma = -3.32;


    double Pion_lower = mSigma*f_nSigmaTOFPionTrack_res + f_nSigmaTOFPionTrack_pos;
    /////double Pion_lower = -4.0;
    ///////double Pion_lower = -3.0;
    /////if (x < 2.0) Pion_lower = (5.0/9.0)*x - (46.0/9.0);

	 return Pion_lower; 
	 }






void Danalysis_test_SDCC_finalise(TString Filename, TString OutputFilename)
{
    TChain *ntp_signal = new TChain("ntp_UnlikeSign");
    ntp_signal->Add(Filename);
    TChain *ntp_Rotated = new TChain("ntp_Rotated");
    ntp_Rotated->Add(Filename);
    TChain *ntp_ME = new TChain("ntp_signal_ME");
    ntp_ME->Add(Filename);
    ////ntp_ME->Add("outputLocal.picoMEtree.sigME.root");
    TChain *ntp_LikeSign = new TChain("ntp_LikeSign");
    ntp_LikeSign->Add(Filename);


    ///Double_t pTforCos  = 0.0;
    /////Double_t cosThetaStarCutD0 = 1.0;
    //////Double_t cosThetaStarCutD0 = 0.90;
    //////Double_t cosThetaStarCutDstar = 0.77;
    //////Double_t cosThetapointD0cut = 0.1;
    ///Double_t dcaD0ToPvCut = 1.7;

    Int_t nSoftPionsTof, nSoftPionsBemc;

    Int_t nD0Candidates = 0;
    Int_t nSoftPionCandidates = 0;

    Int_t nSoftPions = 0;


    Float_t D0_mass, D0_rapidity, D0_decayL, D0_theta, cosTheta, D0_pt, pi1_pt, k_pt,  pi1_dca, k_dca, k_nSigmaTPC, pi1_nSigmaTPC, k_nSigmaTOF, pi1_nSigmaTOF, pi1_TOFinvbeta, k_TOFinvbeta, pi1_betaBase, k_betaBase, D_cosThetaStar, dcaD0ToPv, dcaDaughters, primVz, primVzVpd, k_nHitFit, dcaMax, pi1_eventId, k_eventId,  pi1_nHitFit, pi1_p, k_p, pi1_charge, k_charge, VzVPD_VzTPC, D_phi, D_eta, k_isTOFmatched, pi1_isTOFmatched, k_NHitsDedx, pi1_NHitsDedx, k_btofYLocal, pi1_btofYLocal, k_isBEMCmatched, pi1_isBEMCmatched;
    ntp_signal -> SetBranchAddress("D_mass", &D0_mass);
    ntp_signal -> SetBranchAddress("D_rapidity", &D0_rapidity);
    ntp_signal -> SetBranchAddress("D_decayL", &D0_decayL);
    ntp_signal -> SetBranchAddress("D_theta", &D0_theta);
    //////////ntp_signal -> SetBranchAddress("cosTheta", &cosTheta);
    ntp_signal -> SetBranchAddress("D_pt", &D0_pt);
    ntp_signal -> SetBranchAddress("pi1_pt", &pi1_pt);
    ntp_signal -> SetBranchAddress("k_pt", &k_pt);
    ntp_signal -> SetBranchAddress("pi1_gDCA", &pi1_dca);
    ntp_signal -> SetBranchAddress("k_gDCA", &k_dca);
    ntp_signal -> SetBranchAddress("k_nSigmaTPC", &k_nSigmaTPC);
    ntp_signal -> SetBranchAddress("pi1_nSigmaTPC", &pi1_nSigmaTPC);
    ntp_signal -> SetBranchAddress("k_nSigmaTOF", &k_nSigmaTOF);
    ntp_signal -> SetBranchAddress("pi1_nSigmaTOF", &pi1_nSigmaTOF);
    //////////ntp_signal -> SetBranchAddress("pi1_TOFinvbeta", &pi1_TOFinvbeta);
    //////////ntp_signal -> SetBranchAddress("k_TOFinvbeta", &k_TOFinvbeta);
    //////////ntp_signal -> SetBranchAddress("pi1_betaBase", &pi1_betaBase);
    //////////ntp_signal -> SetBranchAddress("k_betaBase", &k_betaBase);
    ntp_signal -> SetBranchAddress("D_cosThetaStar", &D_cosThetaStar);
    ntp_signal -> SetBranchAddress("dcaD0ToPv", &dcaD0ToPv);
    ntp_signal -> SetBranchAddress("dcaDaughters", &dcaDaughters);
    ntp_signal -> SetBranchAddress("VzTPC", &primVz);
    ntp_signal -> SetBranchAddress("VzVPD", &primVzVpd);
    ntp_signal -> SetBranchAddress("k_nHitsFit", &k_nHitFit);
    ntp_signal -> SetBranchAddress("pi1_nHitsFit", &pi1_nHitFit);
    ntp_signal -> SetBranchAddress("pi1_p", &pi1_p);
    ntp_signal -> SetBranchAddress("k_p", &k_p);
    ntp_signal -> SetBranchAddress("pi1_charge", &pi1_charge);
    ntp_signal -> SetBranchAddress("k_charge", &k_charge);
    //////////ntp_signal -> SetBranchAddress("VzVPD_VzTPC", &VzVPD_VzTPC);
    ntp_signal -> SetBranchAddress("D_phi", &D_phi);
    ntp_signal -> SetBranchAddress("D_eta", &D_eta);
    ntp_signal -> SetBranchAddress("k_isTOFmatched", &k_isTOFmatched);
    ntp_signal -> SetBranchAddress("pi1_isTOFmatched", &pi1_isTOFmatched);
    ntp_signal -> SetBranchAddress("k_NHitsDedx", &k_NHitsDedx);
    ntp_signal -> SetBranchAddress("pi1_NHitsDedx", &pi1_NHitsDedx);
    ntp_signal -> SetBranchAddress("k_btofYLocal", &k_btofYLocal);
    ntp_signal -> SetBranchAddress("pi1_btofYLocal", &pi1_btofYLocal);
    ntp_signal -> SetBranchAddress("k_isBEMCmatched", &k_isBEMCmatched);
    ntp_signal -> SetBranchAddress("pi1_isBEMCmatched", &pi1_isBEMCmatched);
    
    
    
    Float_t rD0_mass, rD0_rapidity, rD0_decayL, rD0_theta, rcosTheta,  rD0_pt, rpi1_pt, rk_pt,  rpi1_dca, rk_dca, rk_nSigmaTPC, rpi1_nSigmaTPC, rk_nSigmaTOF, rpi1_nSigmaTOF, rpi1_TOFinvbeta, rk_TOFinvbeta, rpi1_betaBase, rk_betaBase, rD_cosThetaStar, rdcaD0ToPv, rdcaDaughters, rprimVz, rprimVzVpd, rk_nHitFit, rdcaMax, rpi1_eventId, rk_eventId,  rpi1_nHitFit, rpi1_p, rk_p, rpi1_charge, rk_charge, rVzVPD_VzTPC, rD_phi, rD_eta, rk_isTOFmatched, rpi1_isTOFmatched, rpi1_NHitsDedx, rk_NHitsDedx, rpi1_btofYLocal, rk_btofYLocal, rk_isBEMCmatched, rpi1_isBEMCmatched;
    ntp_Rotated -> SetBranchAddress("D_mass", &rD0_mass);
    ntp_Rotated -> SetBranchAddress("D_rapidity", &rD0_rapidity);
    ntp_Rotated -> SetBranchAddress("D_decayL", &rD0_decayL);
    ntp_Rotated -> SetBranchAddress("D_theta", &rD0_theta);
    //////////ntp_Rotated -> SetBranchAddress("cosTheta", &rcosTheta);
    ntp_Rotated -> SetBranchAddress("D_pt", &rD0_pt);
    ntp_Rotated -> SetBranchAddress("pi1_pt", &rpi1_pt);
    ntp_Rotated -> SetBranchAddress("k_pt", &rk_pt);
    ntp_Rotated -> SetBranchAddress("pi1_gDCA", &rpi1_dca);
    ntp_Rotated -> SetBranchAddress("k_gDCA", &rk_dca);
    ntp_Rotated -> SetBranchAddress("k_nSigmaTPC", &rk_nSigmaTPC);
    ntp_Rotated -> SetBranchAddress("pi1_nSigmaTPC", &rpi1_nSigmaTPC);
    ntp_Rotated -> SetBranchAddress("k_nSigmaTOF", &rk_nSigmaTOF);
    ntp_Rotated -> SetBranchAddress("pi1_nSigmaTOF", &rpi1_nSigmaTOF);
    //////////ntp_Rotated -> SetBranchAddress("pi1_TOFinvbeta", &rpi1_TOFinvbeta);
    //////////ntp_Rotated -> SetBranchAddress("k_TOFinvbeta", &rk_TOFinvbeta);
    //////////ntp_Rotated -> SetBranchAddress("pi1_betaBase", &rpi1_betaBase);
    //////////ntp_Rotated -> SetBranchAddress("k_betaBase", &rk_betaBase);
    ntp_Rotated -> SetBranchAddress("D_cosThetaStar", &rD_cosThetaStar);
    ntp_Rotated -> SetBranchAddress("dcaD0ToPv", &rdcaD0ToPv);
    ntp_Rotated -> SetBranchAddress("dcaDaughters", &rdcaDaughters);
    ntp_Rotated -> SetBranchAddress("VzTPC", &rprimVz);
    ntp_Rotated -> SetBranchAddress("VzVPD", &rprimVzVpd);
    ntp_Rotated -> SetBranchAddress("k_nHitsFit", &rk_nHitFit);
    ntp_Rotated -> SetBranchAddress("pi1_nHitsFit", &rpi1_nHitFit);
    ntp_Rotated -> SetBranchAddress("pi1_p", &rpi1_p);
    ntp_Rotated -> SetBranchAddress("k_p", &rk_p);
    ntp_Rotated -> SetBranchAddress("pi1_charge", &rpi1_charge);
    ntp_Rotated -> SetBranchAddress("k_charge", &rk_charge);
    //////////ntp_Rotated -> SetBranchAddress("VzVPD_VzTPC", &rVzVPD_VzTPC);
    ntp_Rotated -> SetBranchAddress("D_phi", &rD_phi);
    ntp_Rotated -> SetBranchAddress("D_eta", &rD_eta);
    ntp_Rotated -> SetBranchAddress("k_isTOFmatched", &rk_isTOFmatched);
    ntp_Rotated -> SetBranchAddress("pi1_isTOFmatched", &rpi1_isTOFmatched);
    ntp_Rotated -> SetBranchAddress("k_NHitsDedx", &rk_NHitsDedx);
    ntp_Rotated -> SetBranchAddress("pi1_NHitsDedx", &rpi1_NHitsDedx);
    ntp_Rotated -> SetBranchAddress("k_btofYLocal", &rk_btofYLocal);
    ntp_Rotated -> SetBranchAddress("pi1_btofYLocal", &rpi1_btofYLocal);
    ntp_Rotated -> SetBranchAddress("k_isBEMCmatched", &rk_isBEMCmatched);
    ntp_Rotated -> SetBranchAddress("pi1_isBEMCmatched", &rpi1_isBEMCmatched);


    Float_t mD0_mass, mD0_rapidity, mD0_decayL, mD0_theta, mcosTheta, mD0_pt, mpi1_pt, mk_pt,  mpi1_dca, mk_dca, mk_nSigmaTPC, mpi1_nSigmaTPC, mk_nSigmaTOF, mpi1_nSigmaTOF, mpi1_TOFinvbeta, mk_TOFinvbeta, mpi1_betaBase, mk_betaBase, mD_cosThetaStar, mdcaD0ToPv, mdcaDaughters, mprimVz, mprimVzVpd, mk_nHitFit, mdcaMax, mpi1_eventId, mk_eventId,  mpi1_nHitFit, mpi1_p, mk_p, mpi1_charge, mk_charge, mVzVPD_VzTPC, mD_phi, mD_eta, mk_isTOFmatched, mpi1_isTOFmatched;
    ntp_ME -> SetBranchAddress("D_mass", &mD0_mass);
    ntp_ME -> SetBranchAddress("D_rapidity", &mD0_rapidity);
    ntp_ME -> SetBranchAddress("D_decayL", &mD0_decayL);
    ntp_ME -> SetBranchAddress("D_theta", &mD0_theta);
    //////////ntp_ME -> SetBranchAddress("cosTheta", &mcosTheta);
    ntp_ME -> SetBranchAddress("D_pt", &mD0_pt);
    ntp_ME -> SetBranchAddress("pi1_pt", &mpi1_pt);
    ntp_ME -> SetBranchAddress("k_pt", &mk_pt);
    //////////ntp_ME -> SetBranchAddress("pi1_dca", &mpi1_dca);
    //////////ntp_ME -> SetBranchAddress("k_dca", &mk_dca);
    ntp_ME -> SetBranchAddress("k_nSigmaTPC", &mk_nSigmaTPC);
    ntp_ME -> SetBranchAddress("pi1_nSigmaTPC", &mpi1_nSigmaTPC);
    ntp_ME -> SetBranchAddress("k_nSigmaTOF", &mk_nSigmaTOF);
    ntp_ME -> SetBranchAddress("pi1_nSigmaTOF", &mpi1_nSigmaTOF);
    //////////ntp_ME -> SetBranchAddress("pi1_TOFinvbeta", &mpi1_TOFinvbeta);
    //////////ntp_ME -> SetBranchAddress("k_TOFinvbeta", &mk_TOFinvbeta);
    //////////ntp_ME -> SetBranchAddress("pi1_betaBase", &mpi1_betaBase);
    //////////ntp_ME -> SetBranchAddress("k_betaBase", &mk_betaBase);
    ntp_ME -> SetBranchAddress("D_cosThetaStar", &mD_cosThetaStar);
    ntp_ME -> SetBranchAddress("dcaD0ToPv", &mdcaD0ToPv);
    ntp_ME -> SetBranchAddress("dcaDaughters", &mdcaDaughters);
    //////ntp_ME -> SetBranchAddress("primVz", &mprimVz);
    //////ntp_ME -> SetBranchAddress("primVzVpd", &mprimVzVpd);
    //////////ntp_ME -> SetBranchAddress("k_nHitFit", &mk_nHitFit);
    //////////ntp_ME -> SetBranchAddress("pi1_nHitFit", &mpi1_nHitFit);
    ntp_ME -> SetBranchAddress("pi1_p", &mpi1_p);
    ntp_ME -> SetBranchAddress("k_p", &mk_p);
    ntp_ME -> SetBranchAddress("pi1_charge", &mpi1_charge);
    ntp_ME -> SetBranchAddress("k_charge", &mk_charge);
    ///////ntp_ME -> SetBranchAddress("VzVPD_VzTPC", &mVzVPD_VzTPC);
    ntp_ME -> SetBranchAddress("D_phi", &mD_phi);
    ntp_ME -> SetBranchAddress("D_eta", &mD_eta);
    ntp_ME -> SetBranchAddress("k_isTOFmatched", &mk_isTOFmatched);
    ntp_ME -> SetBranchAddress("pi1_isTOFmatched", &mpi1_isTOFmatched);

    Float_t lD0_mass, lD0_rapidity, lD0_decayL, lD0_theta, lcosTheta, lD0_pt, lpi1_pt, lk_pt,  lpi1_dca, lk_dca, lk_nSigmaTPC, lpi1_nSigmaTPC, lk_nSigmaTOF, lpi1_nSigmaTOF, lpi1_TOFinvbeta, lk_TOFinvbeta, lpi1_betaBase, lk_betaBase, lD_cosThetaStar, ldcaD0ToPv, ldcaDaughters, lprimVz, lprimVzVpd, lk_nHitFit, ldcaMax, lpi1_eventId, lk_eventId,  lpi1_nHitFit, lpi1_p, lk_p, lpi1_charge, lk_charge, lVzVPD_VzTPC, lD_phi, lD_eta, lk_isTOFmatched, lpi1_isTOFmatched, lpi1_NHitsDedx, lk_NHitsDedx, lpi1_btofYLocal, lk_btofYLocal, lk_isBEMCmatched, lpi1_isBEMCmatched;
    ntp_LikeSign -> SetBranchAddress("D_mass", &lD0_mass);
    ntp_LikeSign -> SetBranchAddress("D_rapidity", &lD0_rapidity);
    ntp_LikeSign -> SetBranchAddress("D_decayL", &lD0_decayL);
    ntp_LikeSign -> SetBranchAddress("D_theta", &lD0_theta);
    //////////ntp_LikeSign -> SetBranchAddress("cosTheta", &lcosTheta);
    ntp_LikeSign -> SetBranchAddress("D_pt", &lD0_pt);
    ntp_LikeSign -> SetBranchAddress("pi1_pt", &lpi1_pt);
    ntp_LikeSign -> SetBranchAddress("k_pt", &lk_pt);
    ntp_LikeSign -> SetBranchAddress("pi1_gDCA", &lpi1_dca);
    ntp_LikeSign -> SetBranchAddress("k_gDCA", &lk_dca);
    ntp_LikeSign -> SetBranchAddress("k_nSigmaTPC", &lk_nSigmaTPC);
    ntp_LikeSign -> SetBranchAddress("pi1_nSigmaTPC", &lpi1_nSigmaTPC);
    ntp_LikeSign -> SetBranchAddress("k_nSigmaTOF", &lk_nSigmaTOF);
    ntp_LikeSign -> SetBranchAddress("pi1_nSigmaTOF", &lpi1_nSigmaTOF);
    //////////ntp_LikeSign -> SetBranchAddress("pi1_TOFinvbeta", &lpi1_TOFinvbeta);
    //////////ntp_LikeSign -> SetBranchAddress("k_TOFinvbeta", &lk_TOFinvbeta);
    //////////ntp_LikeSign -> SetBranchAddress("pi1_betaBase", &lpi1_betaBase);
    //////////ntp_LikeSign -> SetBranchAddress("k_betaBase", &lk_betaBase);
    ntp_LikeSign -> SetBranchAddress("D_cosThetaStar", &lD_cosThetaStar);
    ntp_LikeSign -> SetBranchAddress("dcaD0ToPv", &ldcaD0ToPv);
    ntp_LikeSign -> SetBranchAddress("dcaDaughters", &ldcaDaughters);
    ntp_LikeSign -> SetBranchAddress("VzTPC", &lprimVz);
    ntp_LikeSign -> SetBranchAddress("VzVPD", &lprimVzVpd);
    ntp_LikeSign -> SetBranchAddress("k_nHitsFit", &lk_nHitFit);
    ntp_LikeSign -> SetBranchAddress("pi1_nHitsFit", &lpi1_nHitFit);
    ntp_LikeSign -> SetBranchAddress("pi1_p", &lpi1_p);
    ntp_LikeSign -> SetBranchAddress("k_p", &lk_p);
    ntp_LikeSign -> SetBranchAddress("pi1_charge", &lpi1_charge);
    ntp_LikeSign -> SetBranchAddress("k_charge", &lk_charge);
    //////////ntp_LikeSign -> SetBranchAddress("VzVPD_VzTPC", &lVzVPD_VzTPC);
    ntp_LikeSign -> SetBranchAddress("D_phi", &lD_phi);
    ntp_LikeSign -> SetBranchAddress("D_eta", &lD_eta);
    ntp_LikeSign -> SetBranchAddress("k_isTOFmatched", &lk_isTOFmatched);
    ntp_LikeSign -> SetBranchAddress("pi1_isTOFmatched", &lpi1_isTOFmatched);
    ntp_LikeSign -> SetBranchAddress("k_NHitsDedx", &lk_NHitsDedx);
    ntp_LikeSign -> SetBranchAddress("pi1_NHitsDedx", &lpi1_NHitsDedx);
    ntp_LikeSign -> SetBranchAddress("k_btofYLocal", &lk_btofYLocal);
    ntp_LikeSign -> SetBranchAddress("pi1_btofYLocal", &lpi1_btofYLocal);
    ntp_LikeSign -> SetBranchAddress("k_isBEMCmatched", &lk_isBEMCmatched);
    ntp_LikeSign -> SetBranchAddress("pi1_isBEMCmatched", &lpi1_isBEMCmatched);

    

    TH2F *D = new TH2F("D","Invariant mass of unlike pairs ",100,0,10.,400,0.5,2.5);   ////////400,0.5,2.5
    TH2F *Dbar = new TH2F("Dbar","Invariant mass of unlike pairs",100,0,10.,400,0.5,2.5); 
    TH2F *DRotate = new TH2F("DRotate","D Rotational background",100,0,10.,400,0.5,2.5);
    TH2F *DbarRotate = new TH2F("DbarRotate","D bar Rotational background",100,0,10.,400,0.5,2.5);
    TH2F *MixedD = new TH2F("MixedD","Mixed event background",100,0,10.,400,0.5,2.5);
    TH2F *MixedDbar = new TH2F("MixedDbar","Mixed event background",100,0,10.,400,0.5,2.5);
    TH2F *LikeBgD = new TH2F("LikeBgD","Combinatorial background",100,0,10.,400,0.5,2.5); 
    TH2F *LikeBgDbar = new TH2F("LikeBgDbar","Combinatorial background",100,0,10.,400,0.5,2.5);


    TH1D *dcaDaughters_unlike = new TH1D("dcaDaughters_unlike","dcaDaughters; dca[cm]; Entries",300,-1.0, 5.0);
    TH1D *D_theta_unlike = new TH1D("D_theta_unlike","D_theta; #theta_{point} [rad]; Entries",300,-1.0, 4.0);
    TH1D *cosTheta_unlike = new TH1D("cosTheta_unlike","cosTheta; cos #theta_{point}; Entries",300,-2.0, 2.0);
    TH1D *D_decayL_unlike = new TH1D("D_decayL_unlike","decay length; decay length [cm]; Entries",2000,-1.0, 1000.0);
    TH1D *dcaD0ToPv_unlike = new TH1D("dcaD0ToPv_unlike","dcaD0ToPv; dcaD0ToPv [cm]; Entries",3000,-1.0, 900.0);
    TH1D *D_cosThetaStar_unlike = new TH1D("D_cosThetaStar_unlike","D_cosThetaStar; cos #theta^{*}; Entries",300,-2.0, 2.0);
    TH1D *D_pt_unlike = new TH1D("D_pt_unlike","D_pt; p_{T} [GeV]; Entries",100,0,10.);
    TH1D *D_mass_unlike = new TH1D("D_mass_unlike","D_mass; invariant mass [GeV]; Entries",400,0.5,2.5);
    TH1D *D_rapidity_unlike = new TH1D("D_rapidity_unlike","D_rapidity; #pi K pair rapidity; Entries",3000,-3.0, 3.0);
    TH1D *D_phi_unlike = new TH1D("D_phi_unlike","D_phi; #pi K pair phi; Entries",100,-5.0, 5.0);
    TH1D *D_eta_unlike = new TH1D("D_eta_unlike","D_eta; #pi K pair eta; Entries",400,-20.0, 20.0);
    TH1D *pi_DCA_unlike = new TH1D("pi_DCA_unlike","pi_DCA; DCA #pi [cm] ; Entries",100,-1.0, 5.0);
    TH1D *k_DCA_unlike = new TH1D("k_DCA_unlike","k_DCA;  DCA k [cm]; Entries",100,-1.0, 5.0);
    TH2D *D_y_D_eta_unlike = new TH2D("D_y_D_eta_unlike","D_{y}_D_{#eta}_unlike; D_{y}; D_{#eta}",2000,-10.,10.,2000,-10.,10.);
    TH2F *h_nSigmaOneOverBetaPi_TOF_TPC_unlike = new TH2F("h_nSigmaOneOverBetaPi_TOF_TPC_unlike", "n#sigma^{#pi}_{1/#beta} after TOF and TPC cut;p[GeV/c] ;n#sigma^{#pi}_{TOF}", 250, 0,10,400,-40,40);
    TH2F *h_nSigmadEdxPi_TOF_TPC_unlike = new TH2F("h_nSigmadEdxPi_TOF_TPC_unlike", "n#sigma^{#pi}_{dE/dx} after TOF and TPC cut; p [GeV/c]; n#sigma^{#pi}_{TPC}", 250, 0,10,400,-40,40);
    TH2F *h_nSigmaOneOverBetaK_TOF_TPC_unlike = new TH2F("h_nSigmaOneOverBetaK_TOF_TPC_unlike", "n#sigma^{K}_{1/#beta} after TOF and TPC cut; p[GeV/c]; n#sigma^{K}_{TOF}", 250, 0,10,400,-40,40);
    TH2F *h_nSigmadEdxK_TOF_TPC_unlike = new TH2F("h_nSigmadEdxK_TOF_TPC_unlike", "n#sigma^{K}_{dE/dx} after TOF and TPC cut; p [GeV/c]; n#sigma^{K}_{TPC}", 250, 0,10,400,-40,40);



    TH1D *dcaDaughters_like = new TH1D("dcaDaughters_like","dcaDaughters; dca[cm]; Entries",300,-1.0, 5.0);
    TH1D *D_theta_like = new TH1D("D_theta_like","D_theta; #theta_{point} [rad]; Entries",300,-1.0, 4.0);
    TH1D *cosTheta_like = new TH1D("cosTheta_like","cosTheta; cos #theta_{point}; Entries",300,-2.0, 2.0);
    TH1D *D_decayL_like = new TH1D("D_decayL_like","decay length; decay length [cm]; Entries",2000,-1.0, 1000.0);
    TH1D *dcaD0ToPv_like = new TH1D("dcaD0ToPv_like","dcaD0ToPv; dcaD0ToPv [cm]; Entries",3000,-1.0, 900.0);
    TH1D *D_cosThetaStar_like = new TH1D("D_cosThetaStar_like","D_cosThetaStar; cos #theta^{*}; Entries",300,-2.0, 2.0);
    TH1D *D_pt_like = new TH1D("D_pt_like","D_pt; p_{T} [GeV]; Entries",100,0,10.);
    TH1D *D_mass_like = new TH1D("D_mass_like","D_mass; invariant mass [GeV]; Entries",400,0.5,2.5);
    TH1D *D_rapidity_like = new TH1D("D_rapidity_like","D_rapidity; #pi K pair rapidity; Entries",3000,-3.0, 3.0);
    TH1D *D_phi_like = new TH1D("D_phi_like","D_phi; #pi K pair phi; Entries",100,-5.0, 5.0);
    TH1D *D_eta_like = new TH1D("D_eta_like","D_eta; #pi K pair eta; Entries",400,-20.0, 20.0);
    TH1D *pi_DCA_like = new TH1D("pi_DCA_like","pi_DCA; DCA #pi [cm] ; Entries",100,-1.0, 5.0);
    TH1D *k_DCA_like = new TH1D("k_DCA_like","k_DCA;  DCA k [cm]; Entries",100,-1.0, 5.0);
    TH2D *D_y_D_eta_like = new TH2D("D_y_D_eta_like","D_{y}_D_{#eta}_like; D_{y}; D_{#eta}",2000,-10.,10.,2000,-10.,10.);
    TH2F *h_nSigmaOneOverBetaPi_TOF_TPC_like = new TH2F("h_nSigmaOneOverBetaPi_TOF_TPC_like", "n#sigma^{#pi}_{1/#beta} after TOF and TPC cut;p[GeV/c] ;n#sigma^{#pi}_{TOF}", 250, 0,10,400,-40,40);
    TH2F *h_nSigmadEdxPi_TOF_TPC_like = new TH2F("h_nSigmadEdxPi_TOF_TPC_like", "n#sigma^{#pi}_{dE/dx} after TOF and TPC cut; p [GeV/c]; n#sigma^{#pi}_{TPC}", 250, 0,10,400,-40,40);
    TH2F *h_nSigmaOneOverBetaK_TOF_TPC_like = new TH2F("h_nSigmaOneOverBetaK_TOF_TPC_like", "n#sigma^{K}_{1/#beta} after TOF and TPC cut; p[GeV/c]; n#sigma^{K}_{TOF}", 250, 0,10,400,-40,40);
    TH2F *h_nSigmadEdxK_TOF_TPC_like = new TH2F("h_nSigmadEdxK_TOF_TPC_like", "n#sigma^{K}_{dE/dx} after TOF and TPC cut; p [GeV/c]; n#sigma^{K}_{TPC}", 250, 0,10,400,-40,40);


    TH1D *dcaDaughters_rot = new TH1D("dcaDaughters_rot","dcaDaughters; dca[cm]; Entries",300,-1.0, 5.0);
    TH1D *D_theta_rot = new TH1D("D_theta_rot","D_theta; #theta_{point} [rad]; Entries",300,-1.0, 4.0);
    TH1D *cosTheta_rot = new TH1D("cosTheta_rot","cosTheta; cos #theta_{point}; Entries",300,-2.0, 2.0);
    TH1D *D_decayL_rot = new TH1D("D_decayL_rot","decay length; decay length [cm]; Entries",2000,-1.0, 1000.0);
    TH1D *dcaD0ToPv_rot = new TH1D("dcaD0ToPv_rot","dcaD0ToPv; dcaD0ToPv [cm]; Entries",3000,-1.0, 900.0);
    TH1D *D_cosThetaStar_rot = new TH1D("D_cosThetaStar_rot","D_cosThetaStar; cos #theta^{*}; Entries",300,-2.0, 2.0);
    TH1D *D_pt_rot = new TH1D("D_pt_rot","D_pt; p_{T} [GeV]; Entries",100,0,10.);
    TH1D *D_mass_rot = new TH1D("D_mass_rot","D_mass; invariant mass [GeV]; Entries",400,0.5,2.5);
    TH1D *D_rapidity_rot = new TH1D("D_rapidity_rot","D_rapidity; #pi K pair rapidity; Entries",3000,-3.0, 3.0);
    TH1D *D_phi_rot = new TH1D("D_phi_rot","D_phi; #pi K pair phi; Entries",100,-5.0, 5.0);
    TH1D *D_eta_rot = new TH1D("D_eta_rot","D_eta; #pi K pair eta; Entries",400,-20.0, 20.0);
    TH1D *pi_DCA_rot = new TH1D("pi_DCA_rot","pi_DCA; DCA #pi [cm] ; Entries",100,-1.0, 5.0);
    TH1D *k_DCA_rot = new TH1D("k_DCA_rot","k_DCA;  DCA k [cm]; Entries",100,-1.0, 5.0);
    TH2D *D_y_D_eta_rot = new TH2D("D_y_D_eta_rot","D_{y}_D_{#eta}_rot; D_{y}; D_{#eta}",2000,-10.,10.,2000,-10.,10.);
    TH2F *h_nSigmaOneOverBetaPi_TOF_TPC_rot = new TH2F("h_nSigmaOneOverBetaPi_TOF_TPC_rot", "n#sigma^{#pi}_{1/#beta} after TOF and TPC cut;p[GeV/c] ;n#sigma^{#pi}_{TOF}", 250, 0,10,400,-40,40);
    TH2F *h_nSigmadEdxPi_TOF_TPC_rot = new TH2F("h_nSigmadEdxPi_TOF_TPC_rot", "n#sigma^{#pi}_{dE/dx} after TOF and TPC cut; p [GeV/c]; n#sigma^{#pi}_{TPC}", 250, 0,10,400,-40,40);
    TH2F *h_nSigmaOneOverBetaK_TOF_TPC_rot = new TH2F("h_nSigmaOneOverBetaK_TOF_TPC_rot", "n#sigma^{K}_{1/#beta} after TOF and TPC cut; p[GeV/c]; n#sigma^{K}_{TOF}", 250, 0,10,400,-40,40);
    TH2F *h_nSigmadEdxK_TOF_TPC_rot = new TH2F("h_nSigmadEdxK_TOF_TPC_rot", "n#sigma^{K}_{dE/dx} after TOF and TPC cut; p [GeV/c]; n#sigma^{K}_{TPC}", 250, 0,10,400,-40,40);

    TH1D *dcaDaughters_mixed = new TH1D("dcaDaughters_mixed","dcaDaughters; dca[cm]; Entries",300,-1.0, 5.0);
    TH1D *D_theta_mixed = new TH1D("D_theta_mixed","D_theta; #theta_{point} [rad]; Entries",300,-1.0, 4.0);
    TH1D *cosTheta_mixed = new TH1D("cosTheta_mixed","cosTheta; cos #theta_{point}; Entries",300,-2.0, 2.0);
    TH1D *D_decayL_mixed = new TH1D("D_decayL_mixed","decay length; decay length [cm]; Entries",2000,-1.0, 1000.0);
    TH1D *dcaD0ToPv_mixed = new TH1D("dcaD0ToPv_mixed","dcaD0ToPv; dcaD0ToPv [cm]; Entries",3000,-1.0, 900.0);
    TH1D *D_cosThetaStar_mixed = new TH1D("D_cosThetaStar_mixed","D_cosThetaStar; cos #theta^{*}; Entries",300,-2.0, 2.0);
    TH1D *D_pt_mixed = new TH1D("D_pt_mixed","D_pt; p_{T} [GeV]; Entries",100,0,10.);
    TH1D *D_mass_mixed = new TH1D("D_mass_mixed","D_mass; invariant mass [GeV]; Entries",400,0.5,2.5);
    TH1D *D_rapidity_mixed = new TH1D("D_rapidity_mixed","D_rapidity; #pi K pair rapidity; Entries",3000,-3.0, 3.0);
    TH1D *D_phi_mixed = new TH1D("D_phi_mixed","D_phi; #pi K pair phi; Entries",100,-5.0, 5.0);
    TH1D *D_eta_mixed = new TH1D("D_eta_mixed","D_eta; #pi K pair eta; Entries",400,-20.0, 20.0);
    TH1D *pi_DCA_mixed = new TH1D("pi_DCA_mixed","pi_DCA; DCA #pi [cm] ; Entries",100,-1.0, 5.0);
    TH1D *k_DCA_mixed = new TH1D("k_DCA_mixed","k_DCA;  DCA k [cm]; Entries",100,-1.0, 5.0);
    TH2D *D_y_D_eta_mixed = new TH2D("D_y_D_eta_mixed","D_{y}_D_{#eta}_mixed; D_{y}; D_{#eta}",2000,-10.,10.,2000,-10.,10.);

    TH2F *hRapidityD = new TH2F("RapidityD","#pi K pair rapidity",100,0.,5.,100,-5.,5.);
    TH2F *hDptPiPt = new TH2F("hDptPiPt","X:D^{0} p_{T}, Y:#pi p_{T}",200,0.,10.,2000,0.,10.);
    TH2F *hDptKPt = new TH2F("hDptKPt","X:D^{0} p_{T}, Y:K p_{T}",200,0.,10.,2000,0.,10.);
    TH2F *hDptPiP = new TH2F("hDptPiP","X:D^{0} p_{T}, Y:#pi momentum",200,0.,10.,200,0.,10.);
    TH2F *hDptKP = new TH2F("hDptKP","X:D^{0} p_{T}, Y:K momentum",200,0.,10.,200,0.,10.);
    TH2F *hPtPionVsKaon = new TH2F("hPtPionVsKaon","hPtPionVsKaon",500,0.,5.,500,0.,5.);
    TH2F *hNHitsdEdxKvsD0pt = new TH2F("hNHitsdEdxKvsD0pt","Number of hits dEdx of Kaon Track vs D0 p_{#perp}",200,0.,10.,50,0,50);
    TH2F *hNHitsdEdxPivsD0pt = new TH2F("hNHitsdEdxPivsD0pt","Number of hits dEdx of Pion Track vs D0 p_{#perp}",200,0.,10.,50,0,50);
    TH2F *hDstarRapidity = new TH2F("DstarRapidity","D* rapidity",500,0.,5.,100,-2.,2.);
    TH2F *hSoftPionPlus = new TH2F("SoftPionPlus","Invariant mass of SoftPion plus ",200,0,20.,400,0.135,0.175); 
    TH2F *hSoftPionMinus = new TH2F("SoftPionMinus","Invariant mass of SoftPion minus ",200,0,20.,400,0.135,0.175);
    TH2F *hD0PtVsSoftPionPt = new TH2F("D0PtVsSoftPionPt","p_{T}(D^{0}) vs p_{T}(soft #pi)",40,0.,2.,400,0.,20.);
    TH2F *hSideBand = new TH2F("SideBand","Side band ",200,0,20.,400,0.135,0.175);
    TH2F *hWrongSign = new TH2F("WrongSign","Wrong Sign",200,0,20.,400,0.135,0.175);
    TH2I *hD0SoftPionCandidatesCorr = new TH2I("D0SoftPionCandidatesCorr","(Soft #pi candidates, D^{0} candidates)",10,0,10,10,0,10);
    TH1F *hVertexZ = new TH1F("VertexZ","Vertex Z",500,-250.,250);
    
    
    Long64_t numberEntr = ntp_signal -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<<"Unlike D0 "<<i<<endl;}
        ntp_signal -> GetEntry(i);
        ///if (D0_pt<=0) continue;
        ///////if ((D0_pt<1.1) || (D0_pt>2.1))  continue;
        hDptPiPt->Fill(D0_pt,pi1_pt);
        hDptKPt->Fill(D0_pt,k_pt);
        hNHitsdEdxKvsD0pt->Fill(D0_pt,k_NHitsDedx);
        hNHitsdEdxPivsD0pt->Fill(D0_pt,pi1_NHitsDedx);

        /////////if (TMath::Abs(primVz) >= 60.0) continue;
        if (TMath::Abs(primVzVpd - primVz) >= 4.0) continue;
        
        if (k_nHitFit < 18 || pi1_nHitFit < 18) continue;
        if (k_dca >= 1.5 || pi1_dca >= 1.5) continue;
        if (k_pt <= 0.2 || pi1_pt <= 0.2) continue;
        /////////if (k_pt >= 1.6 || pi1_pt >= 1.6) continue;
        ///////if (k_NHitsDedx < 10 || pi1_NHitsDedx < 10) continue;
        ///////if (TMath::Abs(k_btofYLocal) > 3.0 || TMath::Abs(pi1_btofYLocal) > 3.0) continue;


        //////if(D0_decayL < 1.0) continue;  ///could be strict
        if(k_nSigmaTPC <= -2.5 || k_nSigmaTPC >= 3.0) continue;
        if(pi1_nSigmaTPC <= -3.0 || pi1_nSigmaTPC >= 3.0) continue;
        
        /*if(k_nSigmaTPC <= -2.5 || k_nSigmaTPC >= 3.0) continue;
        if(pi1_nSigmaTPC <= -2.5 || pi1_nSigmaTPC >= 2.5) continue;*/


        ////if (k_pt < 1.6) {
        if (k_isTOFmatched == 1) {
                if (k_nSigmaTOF <= nSigmaTOFKaonTrack_lower(k_p) || k_nSigmaTOF >= nSigmaTOFKaonTrack_upper(k_p)) continue;
            }
        if (k_isTOFmatched == 0) {
            ///////if (TMath::Abs(k_nSigmaTPC) >= 1.0) continue;
            //////if (k_isBEMCmatched == 0 || TMath::Abs(k_nSigmaTPC) >= 1.0) continue; 
            continue;
            }

        ////if (pi1_pt < 1.6) {
        if (pi1_isTOFmatched == 1) {
                if (pi1_nSigmaTOF <= nSigmaTOFPionTrack_lower(pi1_p) || pi1_nSigmaTOF >= nSigmaTOFPionTrack_upper(pi1_p)) continue;
                }
        if (pi1_isTOFmatched == 0) {
            if (TMath::Abs(pi1_nSigmaTPC) >= 2.5) continue;
            //////if (pi1_isBEMCmatched == 0) continue;
            ///////continue;
            }
        /*   }
        if (pi1_pt >= 1.6) {
            //////if (TMath::Abs(pi1_nSigmaTPC) >= 2.5 && pi1_pt > 1.6) continue;
            if (pi1_isBEMCmatched == 0 || pi1_isTOFmatched == 0) continue;
            }*/


        /*if (pi1_isTOFmatched == 1) {
            if (pi1_nSigmaTOF <= nSigmaTOFPionTrack_lower(pi1_p) || pi1_nSigmaTOF >= nSigmaTOFPionTrack_upper(pi1_p)) continue;
            //////if (pi1_nSigmaTOF <= -2.5 || pi1_nSigmaTOF >= 2.5) continue;
            }
        */

        h_nSigmaOneOverBetaPi_TOF_TPC_unlike->Fill(pi1_p,pi1_nSigmaTOF);
        h_nSigmadEdxPi_TOF_TPC_unlike->Fill(pi1_p,pi1_nSigmaTPC);
        h_nSigmaOneOverBetaK_TOF_TPC_unlike->Fill(k_p,k_nSigmaTOF);
        h_nSigmadEdxK_TOF_TPC_unlike->Fill(k_p,k_nSigmaTPC);
        //////if(pi1_TOFinvbeta < -0.020 || pi1_TOFinvbeta > 0.022) continue;
        //////if(pi1_nSigmaTOF < -2.5 || pi1_nSigmaTOF > 2.5) continue;
        
        /////if (k_TOFinvbeta < DeltaTOFKaonTrack_lower(k_p) || k_TOFinvbeta > DeltaTOFKaonTrack_upper(k_p)) continue;
        /////if (pi1_nSigmaTPC < nSigmaTPCPionTrack_lower(pi1_p) || pi1_nSigmaTPC > nSigmaTPCPionTrack_upper(pi1_p)) continue;
        
        //////if(dcaD0ToPv > 99.9) continue;  ///could be strict
        ///////if(TMath::Abs(cosTheta) > 0.90) continue; ////come back here
        ///////if(TMath::Abs(D_eta) > 1.0) continue;
        
        //////if(D_cosThetaStar > 0.77) continue;  ////come back here
        /////if(pi1_dca > 0.62) continue;  ///????
        /////if(k_dca < 0.38) continue;*/

        if(TMath::Abs(D0_rapidity)>=1.0) continue;
        
        if(pi1_charge == 1 && k_charge == -1 ) {
            ////if(((D0_mass>1.75) && (D0_mass<1.85)) || ((D0_mass>1.9) && (D0_mass<2.0))){
            /////if(((D0_mass>1.78) && (D0_mass<1.83))){
               /*nD0Candidates++;
               hNHitsKvsD0pt->Fill(D0_pt,k_nHitFit);
			   hNHitsPivsD0pt->Fill(D0_pt,pi1_nHitFit);
               hDptPiPt->Fill(D0_pt,pi1_pt);
               hDptKPt->Fill(D0_pt,k_pt);
               hDptPiP->Fill(D0_pt,pi1_p);
			   hDptKP->Fill(D0_pt,k_p);
               hPtPionVsKaon->Fill(pi1_pt,k_pt);*/
               dcaDaughters_unlike->Fill(dcaDaughters);
               D_theta_unlike->Fill(D0_theta);
               cosTheta_unlike->Fill(cosTheta);
               D_decayL_unlike->Fill(D0_decayL);
               dcaD0ToPv_unlike->Fill(dcaD0ToPv);
               D_cosThetaStar_unlike->Fill(D_cosThetaStar);
               D_pt_unlike->Fill(D0_pt);
               D_mass_unlike->Fill(D0_mass);
               D_rapidity_unlike->Fill(D0_rapidity);
               D_phi_unlike->Fill(D_phi);
               D_eta_unlike->Fill(D_eta);
               pi_DCA_unlike->Fill(pi1_dca);
               k_DCA_unlike->Fill(k_dca);
               D_y_D_eta_unlike->Fill(D0_rapidity, D_eta);
            /////}
            /////if(D0_pt>pTforCos) {
            //////if((D_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(D0_theta)) > cosThetapointD0cut)) {
                ///////if((D_cosThetaStar < cosThetaStarCutD0)){
                    D->Fill(D0_pt,D0_mass);
                    //////hRapidityD->Fill(D0_pt,D0_rapidity);
            //////}               
            //////}
            /*else{
					D->Fill(D0_pt,D0_mass);	
                    hRapidityD->Fill(D0_pt,D0_rapidity);
				}
                */
        }

        if(pi1_charge == -1 && k_charge == 1 ) {
             ////if(((D0_mass>1.75) && (D0_mass<1.85)) || ((D0_mass>1.9) && (D0_mass<2.0))){
             ///////if(((D0_mass>1.78) && (D0_mass<1.83))){
               /*hNHitsKvsD0pt->Fill(D0_pt,k_nHitFit);
			   hNHitsPivsD0pt->Fill(D0_pt,pi1_nHitFit);
               hDptPiPt->Fill(D0_pt,pi1_pt);
               hDptKPt->Fill(D0_pt,k_pt);
               hDptPiP->Fill(D0_pt,pi1_p);
			   hDptKP->Fill(D0_pt,k_p);
               hPtPionVsKaon->Fill(pi1_pt,k_pt);
               */

              dcaDaughters_unlike->Fill(dcaDaughters);
               D_theta_unlike->Fill(D0_theta);
               cosTheta_unlike->Fill(cosTheta);
               D_decayL_unlike->Fill(D0_decayL);
               dcaD0ToPv_unlike->Fill(dcaD0ToPv);
               D_cosThetaStar_unlike->Fill(D_cosThetaStar);
               D_pt_unlike->Fill(D0_pt);
               D_mass_unlike->Fill(D0_mass);
               D_rapidity_unlike->Fill(D0_rapidity);
               D_phi_unlike->Fill(D_phi);
               D_eta_unlike->Fill(D_eta);
               pi_DCA_unlike->Fill(pi1_dca);
               k_DCA_unlike->Fill(k_dca);
               D_y_D_eta_unlike->Fill(D0_rapidity, D_eta);
            //////}
            //////if(D0_pt>pTforCos) {
            //////if((D_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(D0_theta)) > cosThetapointD0cut)) {
            ///////if((D_cosThetaStar < cosThetaStarCutD0)){
                    Dbar->Fill(D0_pt,D0_mass);
                    ////hRapidityD->Fill(D0_pt,D0_rapidity);
            ///////}               
            /////}
            /*else{
					Dbar->Fill(D0_pt,D0_mass);	
                    hRapidityD->Fill(D0_pt,D0_rapidity);
				}
                */
        }

    }


    numberEntr = ntp_Rotated -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<< "Rotated D0 "<<i<<endl;}
        ntp_Rotated -> GetEntry(i);
        ///if (D0_pt<=0) continue;
        ////if ((rD0_pt<1.1) || (rD0_pt>2.1))  continue;
        
        ////////if (TMath::Abs(rprimVz) >= 60.0) continue;
        if (TMath::Abs(rprimVzVpd - rprimVz) >= 4.0) continue;

        if (rk_nHitFit < 18 || rpi1_nHitFit < 18) continue;
        if (rk_dca >= 1.5 || rpi1_dca >= 1.5) continue;
        if (rk_pt <= 0.2 || rpi1_pt <= 0.2) continue;
        /////if (rk_pt >= 1.6 || rpi1_pt >= 1.6) continue;
        ////////if (rk_NHitsDedx < 10 || rpi1_NHitsDedx < 10) continue;
        ////////if (TMath::Abs(rk_btofYLocal) > 2.0 || TMath::Abs(rpi1_btofYLocal) > 2.0) continue;

        
        //////if(rD0_decayL < 0.0 || rD0_decayL > 9999999.) continue;  ///could be strict
        //////if(rdcaDaughters > 10.0) continue;
        if(rk_nSigmaTPC <= -2.5 || rk_nSigmaTPC >= 3.0) continue; 
        if (rpi1_nSigmaTPC <= -3.0 || rpi1_nSigmaTPC >= 3.0) continue;
        /*if(rk_nSigmaTPC <= -2.5 || rk_nSigmaTPC >= 3.0) continue;
        if (rpi1_nSigmaTPC <= -2.5 || rpi1_nSigmaTPC >= 2.5) continue;*/
        
        if (rk_isTOFmatched == 1) {
            if (rk_nSigmaTOF <= nSigmaTOFKaonTrack_lower(rk_p) || rk_nSigmaTOF >= nSigmaTOFKaonTrack_upper(rk_p)) continue;
            } 
        if (rk_isTOFmatched == 0) {
            //////if (TMath::Abs(rk_nSigmaTPC) >= 1.0) continue;
            /////if (rk_isBEMCmatched == 0 || TMath::Abs(rk_nSigmaTPC) >= 1.0) continue;
            continue; 
            }

        if (rpi1_isTOFmatched == 1) {
            if (rpi1_nSigmaTOF <= nSigmaTOFPionTrack_lower(rpi1_p) || rpi1_nSigmaTOF >= nSigmaTOFPionTrack_upper(rpi1_p)) continue;
            //////if (rpi1_nSigmaTOF <= -2.5 || rpi1_nSigmaTOF >= 2.5) continue;
            }
        if (rpi1_isTOFmatched == 0) {
            if (TMath::Abs(rpi1_nSigmaTPC) >= 2.5) continue;
            /////continue;
            }

         /*   if (rk_pt < 1.6) {
            if (rk_isTOFmatched == 1) {
                if (rk_nSigmaTOF <= nSigmaTOFKaonTrack_lower(rk_p) || rk_nSigmaTOF >= nSigmaTOFKaonTrack_upper(rk_p)) continue;
                }
            if (rk_isTOFmatched == 0) continue;
            }
           if (rk_pt >= 1.6) {
            ///////if (TMath::Abs(rk_nSigmaTPC) >= 1.0 && rk_pt > 1.6) continue;
            if (rk_isBEMCmatched == 0 || rk_isTOFmatched == 0) continue; 
            }

        if (rpi1_pt < 1.6) {
            if (rpi1_isTOFmatched == 1) {
                if (rpi1_nSigmaTOF <= nSigmaTOFPionTrack_lower(rpi1_p) || rpi1_nSigmaTOF >= nSigmaTOFPionTrack_upper(rpi1_p)) continue;
                }
            if (rpi1_isTOFmatched == 0) continue;
            }
        if (rpi1_pt >= 1.6) {
            //////if (TMath::Abs(rpi1_nSigmaTPC) >= 2.5 && rpi1_pt > 1.6) continue;
            if (rpi1_isBEMCmatched == 0 || rpi1_isTOFmatched == 0) continue;
            }*/


        h_nSigmaOneOverBetaPi_TOF_TPC_rot->Fill(rpi1_p,rpi1_nSigmaTOF);
        h_nSigmadEdxPi_TOF_TPC_rot->Fill(rpi1_p,rpi1_nSigmaTPC);
        h_nSigmaOneOverBetaK_TOF_TPC_rot->Fill(rk_p,rk_nSigmaTOF);
        h_nSigmadEdxK_TOF_TPC_rot->Fill(rk_p,rk_nSigmaTPC);
        //////if(rpi1_nSigmaTOF < -2.5 || rpi1_nSigmaTOF > 2.5) continue;
        //////if(rpi1_TOFinvbeta < -0.020 || rpi1_TOFinvbeta > 0.022) continue;
        ///////if (rk_TOFinvbeta < DeltaTOFKaonTrack_lower(rk_p) || rk_TOFinvbeta > DeltaTOFKaonTrack_upper(rk_p)) continue;
        ///////if (rk_nSigmaTOF < nSigmaTOFKaonTrack_lower(rk_p) || rk_nSigmaTOF > nSigmaTOFKaonTrack_upper(rk_p)) continue;
        /////if (rpi1_nSigmaTPC < nSigmaTPCPionTrack_lower(rpi1_p) || rpi1_nSigmaTPC > nSigmaTPCPionTrack_upper(rpi1_p)) continue;
        /////if(rdcaD0ToPv > 99.9) continue;  ///could be strict
        //////if(TMath::Abs(rcosTheta) > 0.9) continue;  ////come back here
        //////if(TMath::Abs(rD_eta) > 1.0) continue;
        
        /////if(rD_cosThetaStar > 0.77) continue;   ////come back here
        ////if(rpi1_dca > 0.62) continue;
        ////if(rk_dca < 0.38) continue;
        //////if(rD0_mass>3.5 || rD0_mass<0.1) continue;*/

        if((TMath::Abs(rD0_rapidity)>=1.0)) continue;
        
        if(rpi1_charge == 1 && rk_charge == -1) {
            ////if(((rD0_mass>1.75) && (rD0_mass<1.85)) || ((rD0_mass>1.9) && (rD0_mass<2.0))){
            /////if(((rD0_mass>1.78) && (rD0_mass<1.83))){
               dcaDaughters_rot->Fill(rdcaDaughters);
               D_theta_rot->Fill(rD0_theta);
               cosTheta_rot->Fill(rcosTheta);
               D_decayL_rot->Fill(rD0_decayL);
               dcaD0ToPv_rot->Fill(rdcaD0ToPv);
               D_cosThetaStar_rot->Fill(rD_cosThetaStar);
               D_pt_rot->Fill(rD0_pt);
               D_mass_rot->Fill(rD0_mass);
               D_rapidity_rot->Fill(rD0_rapidity);
               D_phi_rot->Fill(rD_phi);
               D_eta_rot->Fill(rD_eta);
               pi_DCA_rot->Fill(rpi1_dca);
               k_DCA_rot->Fill(rk_dca);
               D_y_D_eta_rot->Fill(rD0_rapidity, rD_eta);
            /////}
            /////if((bD_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(bD0_theta)) > cosThetapointD0cut)) {
          ////////if((bD_cosThetaStar < cosThetaStarCutD0)){
            DRotate->Fill(rD0_pt,rD0_mass);
            ///////}
        }
        if(rpi1_charge == -1 && rk_charge == 1) {
             ////if(((rD0_mass>1.75) && (rD0_mass<1.85)) || ((rD0_mass>1.9) && (rD0_mass<2.0))){
            ////if(((rD0_mass>1.78) && (rD0_mass<1.83))){
               dcaDaughters_rot->Fill(rdcaDaughters);
               D_theta_rot->Fill(rD0_theta);
               cosTheta_rot->Fill(rcosTheta);
               D_decayL_rot->Fill(rD0_decayL);
               dcaD0ToPv_rot->Fill(rdcaD0ToPv);
               D_cosThetaStar_rot->Fill(rD_cosThetaStar);
               D_pt_rot->Fill(rD0_pt);
               D_mass_rot->Fill(rD0_mass);
               D_rapidity_rot->Fill(rD0_rapidity);
               D_phi_rot->Fill(rD_phi);
               D_eta_rot->Fill(rD_eta);
               pi_DCA_rot->Fill(rpi1_dca);
               k_DCA_rot->Fill(rk_dca);
               D_y_D_eta_rot->Fill(rD0_rapidity, rD_eta);
            /////}
            /////if((bD_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(bD0_theta)) > cosThetapointD0cut)) {
           ///////if((bD_cosThetaStar < cosThetaStarCutD0)){
            DbarRotate->Fill(rD0_pt,rD0_mass);
            ///////}
        }
   
    }
  
    
    
    numberEntr = ntp_ME -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<< "Mixed Event D0 "<<i<<endl;}
        ntp_ME -> GetEntry(i);
        ///if (D0_pt<=0) continue;
        //////if ((rD0_pt<1.1) || (rD0_pt>2.1))  continue;
        if((TMath::Abs(mD0_rapidity)>=1.0)) continue;
        
        //////if(rD0_decayL < 1.0) continue;  ///could be strict
        //////if(rdcaDaughters > 0.48) continue;
        /////if(mdcaD0ToPv > 99.9) continue;  ///could be strict
        /*if(mk_nSigmaTPC < -2.55 || mk_nSigmaTPC > 2.5) continue;  
        if(mpi1_nSigmaTOF < -2.5 || mpi1_nSigmaTOF > 2.5) continue;
        if (mk_nSigmaTOF < nSigmaTOFKaonTrack_lower(mk_p) || mk_nSigmaTOF > nSigmaTOFKaonTrack_upper(mk_p)) continue;
        
        //////if(TMath::Abs(mcosTheta) > 0.90) continue;  ////come back here
        //////if(TMath::Abs(mD_eta) > 1.0) continue;
        
        ////////if(mD_cosThetaStar > 0.77) continue;   ////come back here
        ////if(rpi1_dca > 0.62) continue;
        ////if(rk_dca < 0.38) continue;*/
        if(mpi1_charge == 1 && mk_charge == -1 ) {
            ////if(((mD0_mass>1.75) && (mD0_mass<1.85)) || ((mD0_mass>1.9) && (mD0_mass<2.0))){
            /////if(((mD0_mass>1.78) && (mD0_mass<1.83))){
               dcaDaughters_mixed->Fill(mdcaDaughters);
               D_theta_mixed->Fill(mD0_theta);
               cosTheta_mixed->Fill(mcosTheta);
               D_decayL_mixed->Fill(mD0_decayL);
               dcaD0ToPv_mixed->Fill(mdcaD0ToPv);
               D_cosThetaStar_mixed->Fill(mD_cosThetaStar);
               D_pt_mixed->Fill(mD0_pt);
               D_mass_mixed->Fill(mD0_mass);
               D_rapidity_mixed->Fill(mD0_rapidity);
               D_phi_mixed->Fill(mD_phi);
               D_eta_mixed->Fill(mD_eta);
               pi_DCA_mixed->Fill(mpi1_dca);
               k_DCA_mixed->Fill(mk_dca);
               D_y_D_eta_mixed->Fill(mD0_rapidity, mD_eta);
            ////}
            ///////if(mD0_pt>pTforCos) {
            ///////if((mD_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(mD0_theta)) > cosThetapointD0cut)) {
            /////if((mD_cosThetaStar < cosThetaStarCutD0)){
                    MixedD->Fill(mD0_pt,mD0_mass);
            /////}               
            ///////}
            /*else{
					MixedD->Fill(mD0_pt,mD0_mass);	
				}*/
                
        }
        if(mpi1_charge == -1 && mk_charge == 1 ) {
                  ////if(((mD0_mass>1.75) && (mD0_mass<1.85)) || ((mD0_mass>1.9) && (mD0_mass<2.0))){
            /////if(((mD0_mass>1.78) && (mD0_mass<1.83))){
               dcaDaughters_mixed->Fill(mdcaDaughters);
               D_theta_mixed->Fill(mD0_theta);
               cosTheta_mixed->Fill(mcosTheta);
               D_decayL_mixed->Fill(mD0_decayL);
               dcaD0ToPv_mixed->Fill(mdcaD0ToPv);
               D_cosThetaStar_mixed->Fill(mD_cosThetaStar);
               D_pt_mixed->Fill(mD0_pt);
               D_mass_mixed->Fill(mD0_mass);
               D_rapidity_mixed->Fill(mD0_rapidity);
               D_phi_mixed->Fill(mD_phi);
               D_eta_mixed->Fill(mD_eta);
               pi_DCA_mixed->Fill(mpi1_dca);
               k_DCA_mixed->Fill(mk_dca);
               D_y_D_eta_mixed->Fill(mD0_rapidity, mD_eta);
            ////}
            ///////if(mD0_pt>pTforCos) {
            ///////if((mD_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(mD0_theta)) > cosThetapointD0cut)) {
            ///////if((mD_cosThetaStar < cosThetaStarCutD0)){
                    MixedDbar->Fill(mD0_pt,mD0_mass);
            //////}               
            ///////}
            /*else{
					MixedDbar->Fill(mD0_pt,mD0_mass);	
				}*/
                
        }
   
    }
    
    


    numberEntr = ntp_LikeSign -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<< "Like D0 "<<i<<endl;}
        ntp_LikeSign -> GetEntry(i);
        ///if (D0_pt<=0) continue;
        //////if ((rD0_pt<1.1) || (rD0_pt>2.1))  continue;

        ///////if (TMath::Abs(lprimVz) >= 60.0) continue;
        if (TMath::Abs(lprimVzVpd - lprimVz) >= 4.0) continue;

        if (lk_nHitFit < 18 || lpi1_nHitFit < 18) continue;
        if (lk_dca >= 1.5 || lpi1_dca >= 1.5) continue;
        if (lk_pt <= 0.2 || lpi1_pt <= 0.2) continue;
        /////////if (lk_pt >= 1.6 || lpi1_pt >= 1.6) continue;
        /////////if (lk_NHitsDedx < 10 || lpi1_NHitsDedx < 10) continue;
        ////////if (TMath::Abs(lk_btofYLocal) > 2.0 || TMath::Abs(lpi1_btofYLocal) > 2.0) continue;




        
        
        //////if(rD0_decayL < 1.0) continue;  ///could be strict
        //////if(rdcaDaughters > 0.48) continue;
        /////if(ldcaD0ToPv > 99.9) continue;  ///could be strict
        if(lk_nSigmaTPC <= -2.5 || lk_nSigmaTPC >= 3.0) continue;
        if (lpi1_nSigmaTPC <= -3.0 || lpi1_nSigmaTPC >= 3.0) continue;
        /*if(lk_nSigmaTPC <= -2.5 || lk_nSigmaTPC >= 3.0) continue;
        if (lpi1_nSigmaTPC <= -2.5 || lpi1_nSigmaTPC >= 2.5) continue;*/

         if (lk_isTOFmatched == 1) {
            if (lk_nSigmaTOF <= nSigmaTOFKaonTrack_lower(lk_p) || lk_nSigmaTOF >= nSigmaTOFKaonTrack_upper(lk_p)) continue;
            } 
        if (lk_isTOFmatched == 0) {
            //////if (TMath::Abs(lk_nSigmaTPC) >= 1.0 && lk_pt > 1.6) continue;
            /////if (lk_isBEMCmatched == 0 || TMath::Abs(lk_nSigmaTPC) >= 1.0) continue;
            continue;
            }

        if (lpi1_isTOFmatched == 1) {
            if (lpi1_nSigmaTOF <= nSigmaTOFPionTrack_lower(lpi1_p) || lpi1_nSigmaTOF >= nSigmaTOFPionTrack_upper(lpi1_p)) continue;
            /////if (lpi1_nSigmaTOF <= -2.5 || lpi1_nSigmaTOF >= 2.5) continue;
            }
        if (lpi1_isTOFmatched == 0) {
            ////if (lpi1_isBEMCmatched == 0) continue;
            if (TMath::Abs(lpi1_nSigmaTPC) >= 2.5) continue;
            /////continue;
            }

        /*if (lk_pt < 1.6) {
            if (lk_isTOFmatched == 1) {
                if (lk_nSigmaTOF <= nSigmaTOFKaonTrack_lower(lk_p) || lk_nSigmaTOF >= nSigmaTOFKaonTrack_upper(lk_p)) continue;
                }
            if (lk_isTOFmatched == 0) continue;
            }
           if (lk_pt >= 1.6) {
            ///////if (TMath::Abs(lk_nSigmaTPC) >= 1.0 && lk_pt > 1.6) continue;
            if (lk_isBEMCmatched == 0 || lk_isTOFmatched == 0) continue; 
            }

        if (lpi1_pt < 1.6) {
            if (lpi1_isTOFmatched == 1) {
                if (lpi1_nSigmaTOF <= nSigmaTOFPionTrack_lower(lpi1_p) || lpi1_nSigmaTOF >= nSigmaTOFPionTrack_upper(lpi1_p)) continue;
                }
            if (lpi1_isTOFmatched == 0) continue;
            }
        if (lpi1_pt >= 1.6) {
            //////if (TMath::Abs(lpi1_nSigmaTPC) >= 2.5 && lpi1_pt > 1.6) continue;
            if (lpi1_isBEMCmatched == 0 || lpi1_isTOFmatched == 0) continue;
            }*/

        h_nSigmaOneOverBetaPi_TOF_TPC_like->Fill(lpi1_p,lpi1_nSigmaTOF);
        h_nSigmadEdxPi_TOF_TPC_like->Fill(lpi1_p,lpi1_nSigmaTPC);
        h_nSigmaOneOverBetaK_TOF_TPC_like->Fill(lk_p,lk_nSigmaTOF);
        h_nSigmadEdxK_TOF_TPC_like->Fill(lk_p,lk_nSigmaTPC);
        ///////if(lpi1_nSigmaTOF < -2.5 || lpi1_nSigmaTOF > 2.5) continue;
        //////if(lpi1_TOFinvbeta < -0.020 || lpi1_TOFinvbeta > 0.022) continue;
        //////if (lk_TOFinvbeta < DeltaTOFKaonTrack_lower(lk_p) || lk_TOFinvbeta > DeltaTOFKaonTrack_upper(lk_p)) continue;
        /////if (lk_nSigmaTOF < nSigmaTOFKaonTrack_lower(lk_p) || lk_nSigmaTOF > nSigmaTOFKaonTrack_upper(lk_p)) continue; 
        //////if (lpi1_nSigmaTPC < nSigmaTPCPionTrack_lower(lpi1_p) || lpi1_nSigmaTPC > nSigmaTPCPionTrack_upper(lpi1_p)) continue;
        ///////if(TMath::Abs(lcosTheta) > 0.90) continue;  ////come back here
        ///////if(TMath::Abs(lD_eta) > 1.0) continue;
        
        /////if(lD_cosThetaStar > 0.77) continue;  ////come back here
        ////if(rpi1_dca > 0.62) continue;
        ////if(rk_dca < 0.38) continue;*/

        if((TMath::Abs(lD0_rapidity)>=1.0)) continue;

        if(lpi1_charge == 1 && lk_charge == 1 ) {
               ////if(((lD0_mass>1.75) && (lD0_mass<1.85)) || ((lD0_mass>1.9) && (lD0_mass<2.0))){
            /////if(((lD0_mass>1.78) && (lD0_mass<1.83))){
               dcaDaughters_like->Fill(ldcaDaughters);
               D_theta_like->Fill(lD0_theta);
               cosTheta_like->Fill(lcosTheta);
               D_decayL_like->Fill(lD0_decayL);
               dcaD0ToPv_like->Fill(ldcaD0ToPv);
               D_cosThetaStar_like->Fill(lD_cosThetaStar);
               D_pt_like->Fill(lD0_pt);
               D_mass_like->Fill(lD0_mass);
               D_rapidity_like->Fill(lD0_rapidity);
               D_phi_like->Fill(lD_phi);
               D_eta_like->Fill(lD_eta);
               pi_DCA_like->Fill(lpi1_dca);
               k_DCA_like->Fill(lk_dca);
               D_y_D_eta_like->Fill(lD0_rapidity, lD_eta);
           ///// }
            ///////if(lD0_pt>pTforCos) {
            ///////if((lD_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(lD0_theta)) > cosThetapointD0cut)) {
            /////if((lD_cosThetaStar < cosThetaStarCutD0)){
                    LikeBgD->Fill(lD0_pt,lD0_mass);
            /////}               
            ///////}
            /*else{
					LikeBgD->Fill(lD0_pt,lD0_mass);	
				}*/
                
        }
        if(lpi1_charge == -1 && lk_charge == -1 ) {

            ////if(((lD0_mass>1.75) && (lD0_mass<1.85)) || ((lD0_mass>1.9) && (lD0_mass<2.0))){
            ////if(((lD0_mass>1.78) && (lD0_mass<1.83))){
               dcaDaughters_like->Fill(ldcaDaughters);
               D_theta_like->Fill(lD0_theta);
               cosTheta_like->Fill(lcosTheta);
               D_decayL_like->Fill(lD0_decayL);
               dcaD0ToPv_like->Fill(ldcaD0ToPv);
               D_cosThetaStar_like->Fill(lD_cosThetaStar);
               D_pt_like->Fill(lD0_pt);
               D_mass_like->Fill(lD0_mass);
               D_rapidity_like->Fill(lD0_rapidity);
               D_phi_like->Fill(lD_phi);
               D_eta_like->Fill(lD_eta);
               pi_DCA_like->Fill(lpi1_dca);
               k_DCA_like->Fill(lk_dca);
               D_y_D_eta_like->Fill(lD0_rapidity, lD_eta);
            /////}
            ///////if(lD0_pt>pTforCos) {
            ///////if((lD_cosThetaStar < cosThetaStarCutD0) && (TMath::Abs(cos(lD0_theta)) > cosThetapointD0cut)) {
            ///////if((lD_cosThetaStar < cosThetaStarCutD0)){
                    LikeBgDbar->Fill(lD0_pt,lD0_mass);
            //////}               
            ///////}
            /*else{
					LikeBgDbar->Fill(lD0_pt,lD0_mass);	
				}*/
                
        }
   
    }
    



   
    TFile* dataRes = new TFile("pp500_D0_prereq_"+OutputFilename,"RECREATE");


    D->Write();
    Dbar->Write();
    DRotate->Write();
    DbarRotate->Write();
    MixedD->Write();
    MixedDbar->Write();
    LikeBgD->Write();
    LikeBgDbar->Write();


    h_nSigmaOneOverBetaPi_TOF_TPC_unlike->Write();
    h_nSigmadEdxPi_TOF_TPC_unlike->Write();
    h_nSigmaOneOverBetaK_TOF_TPC_unlike->Write();
    h_nSigmadEdxK_TOF_TPC_unlike->Write();
    h_nSigmaOneOverBetaPi_TOF_TPC_rot->Write();
    h_nSigmadEdxPi_TOF_TPC_rot->Write();
    h_nSigmaOneOverBetaK_TOF_TPC_rot->Write();
    h_nSigmadEdxK_TOF_TPC_rot->Write();
    h_nSigmaOneOverBetaPi_TOF_TPC_like->Write();
    h_nSigmadEdxPi_TOF_TPC_like->Write();
    h_nSigmaOneOverBetaK_TOF_TPC_like->Write();
    h_nSigmadEdxK_TOF_TPC_like->Write();

    hRapidityD->Write();
    hDptPiPt->Write();
    hDptKPt->Write();
    hDptPiP->Write();
	hDptKP->Write();
    hPtPionVsKaon->Write();
    hNHitsdEdxKvsD0pt->Write();
	hNHitsdEdxPivsD0pt->Write();
    hDstarRapidity->Write();
    hSoftPionPlus->Write();
	hSoftPionMinus->Write();
    hD0PtVsSoftPionPt->Write();
    hSideBand->Write();
    hWrongSign->Write();
    hD0SoftPionCandidatesCorr->Write();

    /*dcaDaughters_unlike->Scale(1.0/dcaDaughters_unlike->GetEntries());
    D_theta_unlike->Scale(1.0/D_theta_unlike->GetEntries());
    cosTheta_unlike->Scale(1.0/cosTheta_unlike->GetEntries());
    D_decayL_unlike->Scale(1.0/D_decayL_unlike->GetEntries());
    dcaD0ToPv_unlike->Scale(1.0/dcaD0ToPv_unlike->GetEntries());
    D_cosThetaStar_unlike->Scale(1.0/D_cosThetaStar_unlike->GetEntries());
    D_pt_unlike->Scale(1.0/D_pt_unlike->GetEntries());
    D_mass_unlike->Scale(1.0/D_mass_unlike->GetEntries());
    D_rapidity_unlike->Scale(1.0/D_rapidity_unlike->GetEntries());
    D_phi_unlike->Scale(1.0/D_phi_unlike->GetEntries());
    D_eta_unlike->Scale(1.0/D_eta_unlike->GetEntries());

    

    dcaDaughters_unlike->Write();
    D_theta_unlike->Write();
    cosTheta_unlike->Write();
    D_decayL_unlike->Write();
    dcaD0ToPv_unlike->Write();
    D_cosThetaStar_unlike->Write();
    D_pt_unlike->Write();
    D_mass_unlike->Write();
    D_rapidity_unlike->Write();
    D_phi_unlike->Write();
    D_eta_unlike->Write();
    pi_DCA_unlike->Write();
    k_DCA_unlike->Write();
    D_y_D_eta_unlike->Write();


    /*dcaDaughters_like->Scale(1.0/dcaDaughters_like->GetEntries());
    D_theta_like->Scale(1.0/D_theta_like->GetEntries());
    cosTheta_like->Scale(1.0/cosTheta_like->GetEntries());
    D_decayL_like->Scale(1.0/D_decayL_like->GetEntries());
    dcaD0ToPv_like->Scale(1.0/dcaD0ToPv_like->GetEntries());
    D_cosThetaStar_like->Scale(1.0/D_cosThetaStar_like->GetEntries());
    D_pt_like->Scale(1.0/D_pt_like->GetEntries());
    D_mass_like->Scale(1.0/D_mass_like->GetEntries());
    D_rapidity_like->Scale(1.0/D_rapidity_like->GetEntries());
    D_phi_like->Scale(1.0/D_phi_like->GetEntries());
    D_eta_like->Scale(1.0/D_eta_like->GetEntries());

    dcaDaughters_like->Write();
    D_theta_like->Write();
    cosTheta_like->Write();
    D_decayL_like->Write();
    dcaD0ToPv_like->Write();
    D_cosThetaStar_like->Write();
    D_pt_like->Write();
    D_mass_like->Write();
    D_rapidity_like->Write();
    D_phi_like->Write();
    D_eta_like->Write();
    pi_DCA_like->Write();
    k_DCA_like->Write();
    D_y_D_eta_like->Write();

    /*dcaDaughters_rot->Scale(1.0/dcaDaughters_rot->GetEntries());
    D_theta_rot->Scale(1.0/D_theta_rot->GetEntries());
    cosTheta_rot->Scale(1.0/cosTheta_rot->GetEntries());
    D_decayL_rot->Scale(1.0/D_decayL_rot->GetEntries());
    dcaD0ToPv_rot->Scale(1.0/dcaD0ToPv_rot->GetEntries());
    D_cosThetaStar_rot->Scale(1.0/D_cosThetaStar_rot->GetEntries());
    D_pt_rot->Scale(1.0/D_pt_rot->GetEntries());
    D_mass_rot->Scale(1.0/D_mass_rot->GetEntries());
    D_rapidity_rot->Scale(1.0/D_rapidity_rot->GetEntries());
    D_phi_rot->Scale(1.0/D_phi_rot->GetEntries());
    D_eta_rot->Scale(1.0/D_eta_rot->GetEntries());

    dcaDaughters_rot->Write();
    D_theta_rot->Write();
    cosTheta_rot->Write();
    D_decayL_rot->Write();
    dcaD0ToPv_rot->Write();
    D_cosThetaStar_rot->Write();
    D_pt_rot->Write();
    D_mass_rot->Write();
    D_rapidity_rot->Write();
    D_phi_rot->Write();
    D_eta_rot->Write();
    pi_DCA_rot->Write();
    k_DCA_rot->Write();
    D_y_D_eta_rot->Write();


    /*dcaDaughters_mixed->Scale(1.0/dcaDaughters_mixed->GetEntries());
    D_theta_mixed->Scale(1.0/D_theta_mixed->GetEntries());
    cosTheta_mixed->Scale(1.0/cosTheta_mixed->GetEntries());
    D_decayL_mixed->Scale(1.0/D_decayL_mixed->GetEntries());
    dcaD0ToPv_mixed->Scale(1.0/dcaD0ToPv_mixed->GetEntries());
    D_cosThetaStar_mixed->Scale(1.0/D_cosThetaStar_mixed->GetEntries());
    D_pt_mixed->Scale(1.0/D_pt_mixed->GetEntries());
    D_mass_mixed->Scale(1.0/D_mass_mixed->GetEntries());
    D_rapidity_mixed->Scale(1.0/D_rapidity_mixed->GetEntries());
    D_phi_mixed->Scale(1.0/D_phi_mixed->GetEntries());
    D_eta_mixed->Scale(1.0/D_eta_mixed->GetEntries());

    dcaDaughters_mixed->Write();
    D_theta_mixed->Write();
    cosTheta_mixed->Write();
    D_decayL_mixed->Write();
    dcaD0ToPv_mixed->Write();
    D_cosThetaStar_mixed->Write();
    D_pt_mixed->Write();
    D_mass_mixed->Write();
    D_rapidity_mixed->Write();
    D_phi_mixed->Write();
    D_eta_mixed->Write();
    pi_DCA_mixed->Write();
    k_DCA_mixed->Write();
    D_y_D_eta_mixed->Write();


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
    D_pt_unlike->SetMarkerStyle(7);
    D_pt_unlike->SetMarkerColor(1);
    D_pt_unlike->SetLineColor(1);
	D_pt_unlike->SetLineWidth(2);
	D_pt_unlike->GetYaxis()->SetTitleFont(42);
	D_pt_unlike->GetYaxis()->SetLabelFont(42);
    D_pt_unlike->GetXaxis()->SetTitleFont(42);
	D_pt_unlike->GetXaxis()->SetLabelFont(42);
	D_pt_unlike->Draw("");
    D_pt_mixed->SetMarkerStyle(7);
    D_pt_mixed->SetMarkerColor(2);
	D_pt_mixed->SetLineColor(2);
	D_pt_mixed->SetLineWidth(2);
    D_pt_mixed->Draw("Csame");
	D_pt_like->SetMarkerStyle(7);
	D_pt_like->SetMarkerColor(4);
	D_pt_like->SetLineColor(4);
	D_pt_like->SetLineWidth(2);
	D_pt_like->Draw("hist C same");
	D_pt_rot->SetMarkerStyle(7);
	D_pt_rot->SetMarkerColor(8);
	D_pt_rot->SetLineColor(8);
	D_pt_rot->SetLineWidth(2);
	D_pt_rot->Draw("hist C same");

    TCanvas *c8 = new TCanvas("c8");
    c8->cd();
    D_mass_unlike->SetMarkerStyle(7);
    D_mass_unlike->SetMarkerColor(1);
    D_mass_unlike->SetLineColor(1);
	D_mass_unlike->SetLineWidth(2);
	D_mass_unlike->GetYaxis()->SetTitleFont(42);
	D_mass_unlike->GetYaxis()->SetLabelFont(42);
    D_mass_unlike->GetXaxis()->SetTitleFont(42);
	D_mass_unlike->GetXaxis()->SetLabelFont(42);
	D_mass_unlike->Draw("");
    D_mass_mixed->SetMarkerStyle(7);
    D_mass_mixed->SetMarkerColor(2);
	D_mass_mixed->SetLineColor(2);
	D_mass_mixed->SetLineWidth(2);
    D_mass_mixed->Draw("Csame");
	D_mass_like->SetMarkerStyle(7);
	D_mass_like->SetMarkerColor(4);
	D_mass_like->SetLineColor(4);
	D_mass_like->SetLineWidth(2);
	D_mass_like->Draw("hist C same");
	D_mass_rot->SetMarkerStyle(7);
	D_mass_rot->SetMarkerColor(8);
	D_mass_rot->SetLineColor(8);
	D_mass_rot->SetLineWidth(2);
	D_mass_rot->Draw("hist C same");


    TCanvas *c9 = new TCanvas("c9");
    c9->cd();
    D_rapidity_unlike->SetMarkerStyle(7);
    D_rapidity_unlike->SetMarkerColor(1);
    D_rapidity_unlike->SetLineColor(1);
	D_rapidity_unlike->SetLineWidth(2);
	D_rapidity_unlike->GetYaxis()->SetTitleFont(42);
	D_rapidity_unlike->GetYaxis()->SetLabelFont(42);
    D_rapidity_unlike->GetXaxis()->SetTitleFont(42);
	D_rapidity_unlike->GetXaxis()->SetLabelFont(42);
	D_rapidity_unlike->Draw("");
    D_rapidity_mixed->SetMarkerStyle(7);
    D_rapidity_mixed->SetMarkerColor(2);
	D_rapidity_mixed->SetLineColor(2);
	D_rapidity_mixed->SetLineWidth(2);
    D_rapidity_mixed->Draw("Csame");
	D_rapidity_like->SetMarkerStyle(7);
	D_rapidity_like->SetMarkerColor(4);
	D_rapidity_like->SetLineColor(4);
	D_rapidity_like->SetLineWidth(2);
	D_rapidity_like->Draw("hist C same");
	D_rapidity_rot->SetMarkerStyle(7);
	D_rapidity_rot->SetMarkerColor(8);
	D_rapidity_rot->SetLineColor(8);
	D_rapidity_rot->SetLineWidth(2);
	D_rapidity_rot->Draw("hist C same");

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
