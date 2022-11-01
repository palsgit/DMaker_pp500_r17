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


void OneOverbeta()
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


 //   ntp->Add("2022-10-10_19-12_D0_7.picoD0AnaMaker.root");


    Float_t D0_theta, D0_mass, D0_pt, D0_decayL, k_pt, pi1_pt, pi1_dca, k_dca, k_nSigma, pi1_nSigma, pi1_TOFinvbeta, k_TOFinvbeta, dcaMax, pi1_eventId, k_eventId, dcaDaughters, D_cosThetaStar, dcaD0ToPv, primVz, primVzVpd, k_nHitFit, pi1_nHitFit, pi1_p, k_p;
    ntp -> SetBranchAddress("D_mass", &D0_mass);
    ntp -> SetBranchAddress("D_decayL", &D0_decayL);
    ntp -> SetBranchAddress("D_theta", &D0_theta);
    ntp -> SetBranchAddress("D_pt", &D0_pt);
    ntp -> SetBranchAddress("pi1_pt", &pi1_pt);
    ntp -> SetBranchAddress("k_pt", &k_pt);
    ntp -> SetBranchAddress("pi1_p", &pi1_p);
    ntp -> SetBranchAddress("k_p", &k_p);
    ntp -> SetBranchAddress("pi1_TOFinvbeta", &pi1_TOFinvbeta);
    ntp -> SetBranchAddress("k_TOFinvbeta", &k_TOFinvbeta);


    //vzor TH1* h1 = new TH1I("h1", "h1 title", 100-počet binů, 0.0, 4.0 -rozsah);
    TH2F *kaonnsigma  = new TH2F("kaonnsigma","kaonnsigma",500,0,3.5,100,-10,10);
    TH2F *pionnsigma  = new TH2F("pionnsigma","pionnsigma",500,0,3.5,100,-10,10);

    TH2F *kaonnsigma50  = new TH2F("kaonnsigma50","kaonnsigma50",50,0,2,100,-10,10);
    TH2F *pionnsigma50  = new TH2F("pionnsigma50","pionnsigma50",50,0,2,100,-10,10);



    float Pion_nsigma = 3.0;
    float Kaon_nsigma = 2.0;
    float Pion_invbeta = 0.03;
    float Kaon_invbeta = 0.03;
    float Pion_nHits = 20;
    float Kaon_nHits = 20;
    float VzTPCVzVPD = 6;
    float cosThetaStarCut = 0.8;  //minimum



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
               // float f_res = pow(0.929095 + 0.0779541 / (k_p - 0.113628), 1.62916);  //sigma
               // float f_pos = pow(-0.0538389 + 0.0439373 / (k_p - 0.0651247), 2.27704);  //mean

              //  float kaon_higher = 3 * f_res + f_pos;
              //  float kaon_lower = -2 * f_res + f_pos;

                float VzVPD_VzTPC = primVzVpd - primVz;
                kaonnsigma->Fill(k_p, nk_TOFinvbeta);
                pionnsigma->Fill(pi1_p, npi1_TOFinvbeta);
                kaonnsigma50->Fill(k_p,nk_TOFinvbeta);
                pionnsigma50->Fill(pi1_p,npi1_TOFinvbeta);

            }
        }
    }


   // TH1D* kaonprojection = kaonnsigma->ProjectionY();Pion_invbeta
   // TH1D* pionprojection = pionnsigma->ProjectionY();



    kaonnsigma50->FitSlicesY();
    TH1D *sigmakaon = (TH1D*)gDirectory->Get("kaonnsigma50_2");
    TH1D *meankaon = (TH1D*)gDirectory->Get("kaonnsigma50_1");
/*
    pionnsigma50->FitSlicesY();
    TH1D *sigmapion = (TH1D*)gDirectory->Get("pionnsigma50_2");
    TH1D *meanpion = (TH1D*)gDirectory->Get("pionnsigma50_1");*/


    TFile* dataRes = new TFile("Projections.root","RECREATE");

    sigmakaon->Write();
    meankaon->Write();
    kaonnsigma50->Write();
    pionnsigma50->Write();


    dataRes->Close();


    TCanvas *c1 = new TCanvas("c1", "Kaon_sigma", 1400, 1100);

    c1->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c1->Update();

    sigmakaon->SetStats(0);
    sigmakaon->SetMarkerStyle(4);
    sigmakaon->SetMarkerColor(2);

    TF1 *fitsigma = new TF1("fitsigma","[0]+[1]/(x+[2])^[3]",0.2,1.);
    fitsigma->SetParameter(0,0.9);
    fitsigma->SetParameter(1,0.02);
    fitsigma->SetParameter(2,0.08);
    fitsigma->SetParameter(3,4.23);

    sigmakaon->Fit(fitsigma,"L","LN",0.2,1.);

    double p0 = fitsigma->GetParameter(0);
    double p1 = fitsigma->GetParameter(1);
    double p2 = fitsigma->GetParameter(2);
    double p3 = fitsigma->GetParameter(3);

    fitsigma->SetLineWidth(0.1);
    fitsigma->SetLineColor(3);
    sigmakaon->Draw("P");
    fitsigma->Draw("same");
    sigmakaon->GetYaxis()->SetTitle("Gaussian resolution of the 1/#beta_{meas}-1/#beta)/0.012");
    sigmakaon->GetXaxis()->SetTitle("p_{K} [GeV/c]");

    sigmakaon->GetXaxis()->SetTitleOffset(1.);
    sigmakaon->GetXaxis()->SetLabelSize(0.03);
    sigmakaon->GetXaxis()->SetTitleSize(0.04);
    sigmakaon->GetXaxis()->SetLabelFont(42);
    sigmakaon->GetXaxis()->SetTitleFont(42);
    sigmakaon->GetXaxis()->CenterTitle(kTRUE);

    sigmakaon->GetYaxis()->SetTitleOffset(1.);
    sigmakaon->GetYaxis()->SetLabelSize(0.03);
    sigmakaon->GetYaxis()->SetTitleSize(0.04);
    sigmakaon->GetYaxis()->SetLabelFont(42);
    sigmakaon->GetYaxis()->SetTitleFont(42);
    sigmakaon->GetYaxis()->CenterTitle(kTRUE);

    c1-> Write(Form("c"));
    c1->SaveAs(Form("Kaon_projection_sigma.pdf"));
    c1->Close();

    cout << "Parametr [0] (sigma kaon) je " << p0 <<endl;
    cout << "Parametr [1] (sigma kaon) je " << p1 <<endl;
    cout << "Parametr [2] (sigma kaon) je " << p2 <<endl;
    cout << "Parametr [3] (sigma kaon) je " << p3 <<endl;



    TCanvas *c2 = new TCanvas("c2", "Kaon_mean", 1400, 1100);

    c1->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c1->Update();

    meankaon->SetStats(0);
    meankaon->SetMarkerStyle(4);
    meankaon->SetMarkerColor(2);

    TF1 *fitmean = new TF1("fitmean","[0]+[1]/(x+[2])^[3]",0.2,1.);
    fitmean->SetParameter(0,0.03);
    fitmean->SetParameter(1,0.0014);
    fitmean->SetParameter(2,0.1);
    fitmean->SetParameter(3,6.9);

    meankaon->Fit(fitmean,"L","LN",0.2,1.);

    double r0 = fitmean->GetParameter(0);
    double r1 = fitmean->GetParameter(1);
    double r2 = fitmean->GetParameter(2);
    double r3 = fitmean->GetParameter(3);

    fitmean->SetLineWidth(0.1);
    fitmean->SetLineColor(3);
    meankaon->Draw("P");
    fitmean->Draw("same");
    meankaon->GetYaxis()->SetTitle("Gaussian resolution of the 1/#beta_{meas}-1/#beta)/0.012");
    meankaon->GetXaxis()->SetTitle("p_{K} [GeV/c]");

    meankaon->GetXaxis()->SetTitleOffset(1.);
    meankaon->GetXaxis()->SetLabelSize(0.03);
    meankaon->GetXaxis()->SetTitleSize(0.04);
    meankaon->GetXaxis()->SetLabelFont(42);
    meankaon->GetXaxis()->SetTitleFont(42);
    meankaon->GetXaxis()->CenterTitle(kTRUE);

    meankaon->GetYaxis()->SetTitleOffset(1.);
    meankaon->GetYaxis()->SetLabelSize(0.03);
    meankaon->GetYaxis()->SetTitleSize(0.04);
    meankaon->GetYaxis()->SetLabelFont(42);
    meankaon->GetYaxis()->SetTitleFont(42);
    meankaon->GetYaxis()->CenterTitle(kTRUE);

    c2-> Write(Form("c"));
    c2->SaveAs(Form("Kaon_projection_mean.pdf"));
    c2->Close();

    cout << "Parametr [0] (mean kaon) je " << r0 <<endl;
    cout << "Parametr [1] (mean kaon) je " << r1 <<endl;
    cout << "Parametr [2] (mean kaon) je " << r2 <<endl;
    cout << "Parametr [3] (mean kaon) je " << r3 <<endl;


    cout<<"Michale jsi dobrej"<<endl;
    cout<<"Hotovo, Jarvis"<<endl;

}

