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



void Danalysis()
{
    TChain *ntp_signal = new TChain("ntp_UnlikeSign");
    ntp_signal->Add("output_all.root");
    TChain *ntp_Rotated = new TChain("ntp_Rotated");
    ntp_Rotated->Add("output_all.root");
    TChain *ntp_background = new TChain("ntp_LikeSign");
    ntp_background->Add("output_all.root");


    Double_t pTforCos  = 0.0;
    Double_t cosThetaStarCutD0 = 1.0;
    Double_t cosThetaStarCutDstar = 0.77;
    Int_t nSoftPionsTof, nSoftPionsBemc;

    Int_t nD0Candidates = 0;
    Int_t nSoftPionCandidates = 0;

    Int_t nSoftPions = 0;


    Float_t D0_mass, D0_rapidity, D0_decayL, D0_theta,  D0_pt, pi1_pt, k_pt,  pi1_dca, k_dca, k_nSigma, pi1_nSigma, pi1_TOFinvbeta, k_TOFinvbeta, pi1_betaBase, k_betaBase, D_cosThetaStar, dcaD0ToPv, dcaDaughters, primVz, primVzVpd, k_nHitFit, dcaMax, pi1_eventId, k_eventId,  pi1_nHitFit, pi1_p, k_p, pi1_charge, k_charge, VzVPD_VzTPC;
    ntp_signal -> SetBranchAddress("D_mass", &D0_mass);
    ntp_signal -> SetBranchAddress("D_rapidity", &D0_rapidity);
    ntp_signal -> SetBranchAddress("D_decayL", &D0_decayL);
    ntp_signal -> SetBranchAddress("D_theta", &D0_theta);
    ntp_signal -> SetBranchAddress("D_pt", &D0_pt);
    ntp_signal -> SetBranchAddress("pi1_pt", &pi1_pt);
    ntp_signal -> SetBranchAddress("k_pt", &k_pt);
    ntp_signal -> SetBranchAddress("pi1_dca", &pi1_dca);
    ntp_signal -> SetBranchAddress("k_dca", &k_dca);
    ntp_signal -> SetBranchAddress("k_nSigma", &k_nSigma);
    ntp_signal -> SetBranchAddress("pi1_nSigma", &pi1_nSigma);
    ntp_signal -> SetBranchAddress("pi1_TOFinvbeta", &pi1_TOFinvbeta);
    ntp_signal -> SetBranchAddress("k_TOFinvbeta", &k_TOFinvbeta);
    ntp_signal -> SetBranchAddress("pi1_betaBase", &pi1_betaBase);
    ntp_signal -> SetBranchAddress("k_betaBase", &k_betaBase);
    ntp_signal -> SetBranchAddress("D_cosThetaStar", &D_cosThetaStar);
    ntp_signal -> SetBranchAddress("dcaD0ToPv", &dcaD0ToPv);
    ntp_signal -> SetBranchAddress("dcaDaughters", &dcaDaughters);
    ntp_signal -> SetBranchAddress("primVz", &primVz);
    ntp_signal -> SetBranchAddress("primVzVpd", &primVzVpd);
    ntp_signal -> SetBranchAddress("k_nHitFit", &k_nHitFit);
    ntp_signal -> SetBranchAddress("pi1_nHitFit", &pi1_nHitFit);
    ntp_signal -> SetBranchAddress("pi1_p", &pi1_p);
    ntp_signal -> SetBranchAddress("k_p", &k_p);
    ntp_signal -> SetBranchAddress("pi1_charge", &pi1_charge);
    ntp_signal -> SetBranchAddress("k_charge", &k_charge);
    ntp_signal -> SetBranchAddress("VzVPD_VzTPC", &VzVPD_VzTPC);
    
    
    
    Float_t bD0_mass, bD0_rapidity, bD0_decayL, bD0_theta,  bD0_pt, bpi1_pt, bk_pt,  bpi1_dca, bk_dca, bk_nSigma, bpi1_nSigma, bpi1_TOFinvbeta, bk_TOFinvbeta, bpi1_betaBase, bk_betaBase, bD_cosThetaStar, bdcaD0ToPv, bdcaDaughters, bprimVz, bprimVzVpd, bk_nHitFit, bdcaMax, bpi1_eventId, bk_eventId,  bpi1_nHitFit, bpi1_p, bk_p, bpi1_charge, bk_charge, bVzVPD_VzTPC;
    ntp_Rotated -> SetBranchAddress("D_mass", &bD0_mass);
    ntp_Rotated -> SetBranchAddress("D_rapidity", &bD0_rapidity);
    ntp_Rotated -> SetBranchAddress("D_decayL", &bD0_decayL);
    ntp_Rotated -> SetBranchAddress("D_theta", &bD0_theta);
    ntp_Rotated -> SetBranchAddress("D_pt", &bD0_pt);
    ntp_Rotated -> SetBranchAddress("pi1_pt", &bpi1_pt);
    ntp_Rotated -> SetBranchAddress("k_pt", &bk_pt);
    ntp_Rotated -> SetBranchAddress("pi1_dca", &bpi1_dca);
    ntp_Rotated -> SetBranchAddress("k_dca", &bk_dca);
    ntp_Rotated -> SetBranchAddress("k_nSigma", &bk_nSigma);
    ntp_Rotated -> SetBranchAddress("pi1_nSigma", &bpi1_nSigma);
    ntp_Rotated -> SetBranchAddress("pi1_TOFinvbeta", &bpi1_TOFinvbeta);
    ntp_Rotated -> SetBranchAddress("k_TOFinvbeta", &bk_TOFinvbeta);
    ntp_Rotated -> SetBranchAddress("pi1_betaBase", &bpi1_betaBase);
    ntp_Rotated -> SetBranchAddress("k_betaBase", &bk_betaBase);
    ntp_Rotated -> SetBranchAddress("D_cosThetaStar", &bD_cosThetaStar);
    ntp_Rotated -> SetBranchAddress("dcaD0ToPv", &bdcaD0ToPv);
    ntp_Rotated -> SetBranchAddress("dcaDaughters", &bdcaDaughters);
    ntp_Rotated -> SetBranchAddress("primVz", &bprimVz);
    ntp_Rotated -> SetBranchAddress("primVzVpd", &bprimVzVpd);
    ntp_Rotated -> SetBranchAddress("k_nHitFit", &bk_nHitFit);
    ntp_Rotated -> SetBranchAddress("pi1_nHitFit", &bpi1_nHitFit);
    ntp_Rotated -> SetBranchAddress("pi1_p", &bpi1_p);
    ntp_Rotated -> SetBranchAddress("k_p", &bk_p);
    ntp_Rotated -> SetBranchAddress("pi1_charge", &bpi1_charge);
    ntp_Rotated -> SetBranchAddress("k_charge", &bk_charge);
    ntp_Rotated -> SetBranchAddress("VzVPD_VzTPC", &bVzVPD_VzTPC);


    Float_t lD0_mass, lD0_rapidity, lD0_decayL, lD0_theta,  lD0_pt, lpi1_pt, lk_pt,  lpi1_dca, lk_dca, lk_nSigma, lpi1_nSigma, lpi1_TOFinvbeta, lk_TOFinvbeta, lpi1_betaBase, lk_betaBase, lD_cosThetaStar, ldcaD0ToPv, ldcaDaughters, lprimVz, lprimVzVpd, lk_nHitFit, ldcaMax, lpi1_eventId, lk_eventId,  lpi1_nHitFit, lpi1_p, lk_p, lpi1_charge, lk_charge, lVzVPD_VzTPC;
    ntp_background -> SetBranchAddress("D_mass", &lD0_mass);
    ntp_background -> SetBranchAddress("D_rapidity", &lD0_rapidity);
    ntp_background -> SetBranchAddress("D_decayL", &lD0_decayL);
    ntp_background -> SetBranchAddress("D_theta", &lD0_theta);
    ntp_background -> SetBranchAddress("D_pt", &lD0_pt);
    ntp_background -> SetBranchAddress("pi1_pt", &lpi1_pt);
    ntp_background -> SetBranchAddress("k_pt", &lk_pt);
    ntp_background -> SetBranchAddress("pi1_dca", &lpi1_dca);
    ntp_background -> SetBranchAddress("k_dca", &lk_dca);
    ntp_background -> SetBranchAddress("k_nSigma", &lk_nSigma);
    ntp_background -> SetBranchAddress("pi1_nSigma", &lpi1_nSigma);
    ntp_background -> SetBranchAddress("pi1_TOFinvbeta", &lpi1_TOFinvbeta);
    ntp_background -> SetBranchAddress("k_TOFinvbeta", &lk_TOFinvbeta);
    ntp_background -> SetBranchAddress("pi1_betaBase", &lpi1_betaBase);
    ntp_background -> SetBranchAddress("k_betaBase", &lk_betaBase);
    ntp_background -> SetBranchAddress("D_cosThetaStar", &lD_cosThetaStar);
    ntp_background -> SetBranchAddress("dcaD0ToPv", &ldcaD0ToPv);
    ntp_background -> SetBranchAddress("dcaDaughters", &ldcaDaughters);
    ntp_background -> SetBranchAddress("primVz", &lprimVz);
    ntp_background -> SetBranchAddress("primVzVpd", &lprimVzVpd);
    ntp_background -> SetBranchAddress("k_nHitFit", &lk_nHitFit);
    ntp_background -> SetBranchAddress("pi1_nHitFit", &lpi1_nHitFit);
    ntp_background -> SetBranchAddress("pi1_p", &lpi1_p);
    ntp_background -> SetBranchAddress("k_p", &lk_p);
    ntp_background -> SetBranchAddress("pi1_charge", &lpi1_charge);
    ntp_background -> SetBranchAddress("k_charge", &lk_charge);
    ntp_background -> SetBranchAddress("VzVPD_VzTPC", &lVzVPD_VzTPC);
    

    TH2F *D = new TH2F("D","Invariant mass of unlike pairs ",100,0,10.,400,0.5,2.5);
    TH2F *Dbar = new TH2F("Dbar","Invariant mass of unlike pairs",100,0,10.,400,0.5,2.5); 
    TH2F *DRotate = new TH2F("DRotate","D Rotational background",100,0,10.,400,0.5,2.5);
    TH2F *DbarRotate = new TH2F("DbarRotate","D bar Rotational background",100,0,10.,400,0.5,2.5);
    TH2F *MixedD = new TH2F("MixedD","Mixed event background",100,0,10.,400,0.5,2.5);
    TH2F *MixedDbar = new TH2F("MixedDbar","Mixed event background",100,0,10.,400,0.5,2.5);
    TH2F *LikeBgD = new TH2F("LikeBgD","Combinatorial background",100,0,10.,400,0.5,2.5); 
    TH2F *LikeBgDbar = new TH2F("LikeBgDbar","Combinatorial background",100,0,10.,400,0.5,2.5);
    TH2F *hRapidityD = new TH2F("RapidityD","#pi K pair rapidity",100,0.,5.,100,-5.,5.);
    TH2F *hDptPiPt = new TH2F("hDptPiPt","X:D^{0} p_{T}, Y:#pi p_{T}",200,0.,10.,200,0.,10.);
    TH2F *hDptKPt = new TH2F("hDptKPt","X:D^{0} p_{T}, Y:K p_{T}",200,0.,10.,200,0.,10.);
    TH2F *hDptPiP = new TH2F("hDptPiP","X:D^{0} p_{T}, Y:#pi momentum",200,0.,10.,200,0.,10.);
    TH2F *hDptKP = new TH2F("hDptKP","X:D^{0} p_{T}, Y:K momentum",200,0.,10.,200,0.,10.);
    TH2F *hPtPionVsKaon = new TH2F("hPtPionVsKaon","hPtPionVsKaon",500,0.,5.,500,0.,5.);
    TH2F *hNHitsKvsD0pt = new TH2F("NHitsKvsD0pt","Number of hits of Kaon Track vs D0 p_{#perp}",200,0.,10.,50,0,50);
    TH2F *hNHitsPivsD0pt = new TH2F("NHitsPivsD0pt","Number of hits of Pion Track vs D0 p_{#perp}",200,0.,10.,50,0,50);
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
        if(D0_pt<=0) continue;
        if(TMath::Abs(D0_rapidity)>1.0) continue;
        if(pi1_charge == 1 && k_charge == -1 ) {
            if(D0_mass>1.84 && D0_mass<1.89 && D0_pt>pTforCos && D_cosThetaStar<cosThetaStarCutDstar){
               nD0Candidates++;
               hNHitsKvsD0pt->Fill(D0_pt,k_nHitFit);
			   hNHitsPivsD0pt->Fill(D0_pt,pi1_nHitFit);
               hDptPiPt->Fill(D0_pt,pi1_pt);
               hDptKPt->Fill(D0_pt,k_pt);
               hDptPiP->Fill(D0_pt,pi1_p);
			   hDptKP->Fill(D0_pt,k_p);
               hPtPionVsKaon->Fill(pi1_pt,k_pt);
            }
            if(D0_pt>pTforCos) {
                if(D_cosThetaStar<cosThetaStarCutD0) {
                    D->Fill(D0_pt,D0_mass);
                    hRapidityD->Fill(D0_pt,D0_rapidity);
                }               
            }
            else{
					D->Fill(D0_pt,D0_mass);	
                    hRapidityD->Fill(D0_pt,D0_rapidity);
				}
        }

        if(pi1_charge == -1 && k_charge == 1 ) {
             if(D0_mass>1.84 && D0_mass<1.89 && D0_pt>pTforCos && D_cosThetaStar<cosThetaStarCutDstar){
               hNHitsKvsD0pt->Fill(D0_pt,k_nHitFit);
			   hNHitsPivsD0pt->Fill(D0_pt,pi1_nHitFit);
               hDptPiPt->Fill(D0_pt,pi1_pt);
               hDptKPt->Fill(D0_pt,k_pt);
               hDptPiP->Fill(D0_pt,pi1_p);
			   hDptKP->Fill(D0_pt,k_p);
               hPtPionVsKaon->Fill(pi1_pt,k_pt);
            }
            if(D0_pt>pTforCos) {
                if(D_cosThetaStar<cosThetaStarCutD0) {
                    Dbar->Fill(D0_pt,D0_mass);
                    hRapidityD->Fill(D0_pt,D0_rapidity);
                }               
            }
            else{
					Dbar->Fill(D0_pt,D0_mass);	
                    hRapidityD->Fill(D0_pt,D0_rapidity);
				}
        }

    }


    numberEntr = ntp_Rotated -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<< "Rotated D0 "<<i<<endl;}
        ntp_Rotated -> GetEntry(i);
        if(bD0_pt<=0) continue;
        if(TMath::Abs(bD0_rapidity)>1.0) continue;
        if(bpi1_charge == 1 && bk_charge == -1) {
            DRotate->Fill(bD0_pt,bD0_mass);
        }
        if(bpi1_charge == -1 && bk_charge == 1) {
            DbarRotate->Fill(bD0_pt,bD0_mass);
        }
   
    }

    numberEntr = ntp_background -> GetEntries();
    cout<<"Number of entries in Ntuple: "<<numberEntr<<endl;
    for (Long64_t i = 0; i < numberEntr; i++) {
        if (i%10000000==0) {cout<< "Like D0 "<<i<<endl;}
        ntp_background -> GetEntry(i);
        if(lD0_pt<=0) continue;
        if(TMath::Abs(lD0_rapidity)>1.0) continue;
        if(lpi1_charge == 1 && lk_charge == 1 ) {
            if(lD0_pt>pTforCos) {
                if(lD_cosThetaStar<cosThetaStarCutD0) {
                    LikeBgD->Fill(lD0_pt,lD0_mass);
                }               
            }
            else{
					LikeBgD->Fill(lD0_pt,lD0_mass);	
				}
        }
        if(lpi1_charge == -1 && lk_charge == -1 ) {
            if(lD0_pt>pTforCos) {
                if(lD_cosThetaStar<cosThetaStarCutD0) {
                    LikeBgDbar->Fill(lD0_pt,lD0_mass);
                }               
            }
            else{
					LikeBgDbar->Fill(lD0_pt,lD0_mass);	
				}
        }
   
    }
   
    TFile* dataRes = new TFile("pp500_D0_hybridTOF_prereqhigheta_all.root","RECREATE");


    D->Write();
    Dbar->Write();
    DRotate->Write();
    DbarRotate->Write();
    MixedD->Write();
    MixedDbar->Write();
    LikeBgD->Write();
    LikeBgDbar->Write();
    hRapidityD->Write();
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
    hD0SoftPionCandidatesCorr->Write();


    dataRes->Close();

    
    cout<<"Finished"<<endl;


}

