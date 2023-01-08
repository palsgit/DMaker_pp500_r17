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
#include "TFile.h"
#include "TNtuple.h"

void KaonPionBeta(){

    TFile* data = new TFile("pp500_D0.root");
    TH2F *h1  = (TH2F*) data -> Get("kaonnsigma;1");
    h1->SetTitle(" ; ");


    TCanvas *c1 = new TCanvas("c1", "Kaon 1/#beta", 1400, 1100);

    c1->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c1->Update();
    h1->SetStats(0);
    h1->Draw("colz");

    TF1 *f_res = new TF1("f_res","[0]+[1]/(x+[2])^[3]",0.2,3.5); // sigma
    f_res->SetParameters(0.929095,0.0779541,-0.113628,1.62916);
    TF1 *f_pos= new TF1("f_pos","[0]+[1]/(x+[2])^[3]",0.2,3.5); // mean
    f_pos->SetParameters(-0.0538389,0.0439373,-0.0651247,2.27704);


    TF1 *f1 = new TF1("Kupper","3*f_res+f_pos",0.2,3.5);
    f1->SetLineColor(kBlue);
    f1->Draw("same");

    TF1 *f2 = new TF1("Klower","-2*f_res+f_pos",0.2,3.5);
    f2->SetLineColor(kBlue);
    f2->Draw("same");




    h1->GetYaxis()->SetTitle("#frac{(1/#beta_{TOF}-1/#beta)}{0.012}");
    h1->GetXaxis()->SetTitle("#it{p}_{K} [GeV/#it{c}]");

    h1->GetXaxis()->SetTitleOffset(1.);
    h1->GetXaxis()->SetLabelSize(0.03);
    h1->GetXaxis()->SetTitleSize(0.04);
    h1->GetXaxis()->SetLabelFont(42);
    h1->GetXaxis()->SetTitleFont(42);
    h1->GetXaxis()->CenterTitle(kTRUE);


    h1->GetYaxis()->SetRangeUser(-7.5,7.5);

    h1->GetYaxis()->SetTitleOffset(1.);
    h1->GetYaxis()->SetLabelSize(0.03);
    h1->GetYaxis()->SetTitleSize(0.04);
    h1->GetYaxis()->SetLabelFont(42);
    h1->GetYaxis()->SetTitleFont(42);

    c1->SetLogz();

    TLine *l1 = new TLine(0,2.5,3.5,2.5);
    l1->Draw("same");
    TLine *l2 = new TLine(0,-2.5,3.5,-2.5);
    l2->Draw("same");

    h1->GetYaxis()->CenterTitle(kTRUE);


    c1->SaveAs(Form("pp500_Beta_kaons.pdf"));
    c1->Close();



    TH2F *h2  = (TH2F*) data -> Get("pionnsigma;1");
    h2->SetTitle(" ; ");



    TCanvas *c2 = new TCanvas("c2", "Pion 1/#beta", 1400, 1100);

    c2->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c2->Update();
    h2->SetStats(0);
    h2->Draw("colz");


    TF1 *f5 = new TF1("Pupper","6 - 2.666666666666666*x",0.2,1.5);
    f5->SetLineColor(kBlue);
    f5->Draw("same");

    TF1 *f6 = new TF1("Plower","-6+2.666666666666666*x",0.2,1.5);
    f6->SetLineColor(kBlue);
    f6->Draw("same");

    TF1 *f7 = new TF1("PupperL","2",1.5,3.5);
    f7->SetLineColor(kBlue);
    f7->Draw("same");

    TF1 *f8 = new TF1("PlowerL","-2",1.5,3.5);
    f8->SetLineColor(kBlue);
    f8->Draw("same");



    h2->GetYaxis()->SetTitle("#frac{(1/#beta_{TOF}-1/#beta)}{0.012}");
    h2->GetXaxis()->SetTitle("#it{p}_{#pi} [GeV/#it{c}]");

    h2->GetXaxis()->SetTitleOffset(1.);
    h2->GetXaxis()->SetLabelSize(0.03);
    h2->GetXaxis()->SetTitleSize(0.04);
    h2->GetXaxis()->SetLabelFont(42);
    h2->GetXaxis()->SetTitleFont(42);
    h2->GetXaxis()->CenterTitle(kTRUE);


    h2->GetYaxis()->SetRangeUser(-7.5,7.5);

    h2->GetYaxis()->SetTitleOffset(1.);
    h2->GetYaxis()->SetLabelSize(0.03);
    h2->GetYaxis()->SetTitleSize(0.04);
    h2->GetYaxis()->SetLabelFont(42);
    h2->GetYaxis()->SetTitleFont(42);
    // h1->GetYaxis()->SetRangeUser(0,0.06);

    c2->SetLogz();


    h2->GetYaxis()->CenterTitle(kTRUE);


    TLine *l3 = new TLine(0.,2.5,3.5,2.5);
    l3->Draw("same");
    TLine *l4 = new TLine(0.,-2.5,3.5,-2.5);
    l4->Draw("same");


    c2->SaveAs(Form("pp500_Beta_pions.pdf"));
    c2->Close();
}
