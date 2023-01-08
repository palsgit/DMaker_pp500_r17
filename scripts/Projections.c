#include <iostream>
#include "TFile.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"
#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TNtuple.h"

void Projections(){

    TFile* data = new TFile("pp500_D0.root");
    TH1F *h1  = (TH1F*) data -> Get("kaonnsigma50_1;1");


    TCanvas *c1 = new TCanvas("c1", "Kaon_mean", 1400, 1000);

    c1->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c1->Update();
    h1->SetStats(0);
    h1->SetMarkerStyle(4);
    h1->SetMarkerColor(kBlue);


    TF1 *fitmean = new TF1("fitmean","[0]+[1]/(x+[2])^[3]",0.2,1.);
    fitmean->SetParameter(0,0.03);
    fitmean->SetParameter(1,0.0014);
    fitmean->SetParameter(2,0.1);
    fitmean->SetParameter(3,6.9);

    h1->Fit(fitmean,"L","LN",0.2,1.);

    double p0 = fitmean->GetParameter(0);
    double p1 = fitmean->GetParameter(1);
    double p2 = fitmean->GetParameter(2);
    double p3 = fitmean->GetParameter(3);

    fitmean->SetLineWidth(0.5);
    fitmean->SetLineColor(3);
    h1->Draw("PE");
    fitmean->Draw("same");
    h1->GetYaxis()->SetTitle("Gaussian means of the #frac{(1/#beta_{TOF}-1/#beta)}{0.012}");
    h1->GetXaxis()->SetTitle("#it{p}_{K} [GeV/#it{c}]");

    TLegend *l1 = new TLegend(0.55,0.69, 0.73, 0.89,"","brNDC");
    l1->AddEntry((TObject*)0, "", "");
    l1->AddEntry((TObject*)0, "THIS THESIS", "");
    l1->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");

    l1->SetFillStyle(0);
    l1->SetLineColor(0);
    l1->SetTextSize(0.03);
    l1->Draw("same");

    h1->GetXaxis()->SetTitleOffset(1.);
    h1->GetXaxis()->SetLabelSize(0.03);
    h1->GetXaxis()->SetTitleSize(0.04);
    h1->GetXaxis()->SetLabelFont(42);
    h1->GetXaxis()->SetTitleFont(42);
    h1->GetXaxis()->CenterTitle(kTRUE);
    h1->GetXaxis()->SetRangeUser(0.2,1.2);



    h1->GetYaxis()->SetTitleOffset(1.);
    h1->GetYaxis()->SetLabelSize(0.03);
    h1->GetYaxis()->SetTitleSize(0.04);
    h1->GetYaxis()->SetLabelFont(42);
    h1->GetYaxis()->SetTitleFont(42);
    // histo->GetYaxis()->SetRangeUser(0,0.06);


    h1->GetYaxis()->CenterTitle(kTRUE);



    c1-> Write(Form("c"));
    c1->SaveAs(Form("pp500_Kaon_projection_mean.pdf"));
    c1->Close();

    cout << "Parametr [0] je " << p0 <<endl;
    cout << "Parametr [1] je " << p1 <<endl;
    cout << "Parametr [2] je " << p2 <<endl;
    cout << "Parametr [3] je " << p3 <<endl;

    TH2F *h2  = (TH2F*) data -> Get("kaonnsigma50_2;1");


    TCanvas *c2 = new TCanvas("c2", "Kaon_sigma", 1400, 1000);

    c2->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c2->Update();
    h2->SetStats(0);
    h2->SetMarkerStyle(4);
    h2->SetMarkerColor(kBlue);


    TF1 *fitsigma = new TF1("fitsigma","[0]+[1]/(x+[2])^[3]",0.2,1.);
    fitsigma->SetParameter(0,0.9);
    fitsigma->SetParameter(1,0.02);
    fitsigma->SetParameter(2,0.08);
    fitsigma->SetParameter(3,4.23);

    h2->Fit(fitsigma,"L","LN",0.2,1.);

    double r0 = fitsigma->GetParameter(0);
    double r1 = fitsigma->GetParameter(1);
    double r2 = fitsigma->GetParameter(2);
    double r3 = fitsigma->GetParameter(3);

    fitsigma->SetLineWidth(0.5);
    fitsigma->SetLineColor(3);
    h2->Draw("PE");
    fitsigma->Draw("same");
    h2->GetYaxis()->SetTitle("Gaussian resolution of the #frac{(1/#beta_{TOF}-1/#beta)}{0.012}");
    h2->GetXaxis()->SetTitle("#it{p}_{K} [GeV/#it{c}]");

    TLegend *l2 = new TLegend(0.55,0.69, 0.73, 0.89,"","brNDC");
    l2->AddEntry((TObject*)0, "", "");
    l2->AddEntry((TObject*)0, "THIS THESIS", "");
    l2->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");

    l2->SetFillStyle(0);
    l2->SetLineColor(0);
    l2->SetTextSize(0.03);
    l2->Draw("same");

    h2->GetXaxis()->SetTitleOffset(1.);
    h2->GetXaxis()->SetLabelSize(0.03);
    h2->GetXaxis()->SetTitleSize(0.04);
    h2->GetXaxis()->SetLabelFont(42);
    h2->GetXaxis()->SetTitleFont(42);
    h2->GetXaxis()->CenterTitle(kTRUE);
    h2->GetXaxis()->SetRangeUser(0.2,1.2);



    h2->GetYaxis()->SetTitleOffset(1.);
    h2->GetYaxis()->SetLabelSize(0.03);
    h2->GetYaxis()->SetTitleSize(0.04);
    h2->GetYaxis()->SetLabelFont(42);
    h2->GetYaxis()->SetTitleFont(42);
    // histo->GetYaxis()->SetRangeUser(0,0.06);


    h2->GetYaxis()->CenterTitle(kTRUE);



    c2-> Write(Form("c"));
    c2->SaveAs(Form("pp500_Kaon_projection_sigma.pdf"));
    c2->Close();

    cout << "Parametr [0] je " << r0 <<endl;
    cout << "Parametr [1] je " << r1 <<endl;
    cout << "Parametr [2] je " << r2 <<endl;
    cout << "Parametr [3] je " << r3 <<endl;
    cout << "Pokud tohle funguje, jsi dobrej Michale" <<endl;


}
