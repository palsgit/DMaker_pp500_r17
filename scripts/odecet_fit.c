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

Bool_t reject;
Double_t fline(Double_t *x, Double_t *par)
{
    if (reject && x[0] > 1.833 && x[0] < 1.899) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0];
}

Double_t fquadratic(Double_t *x, Double_t *par)
{
    if (reject && x[0] > 1.833 && x[0] < 1.899) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t fcubic(Double_t *x, Double_t *par)
{
    if (reject && x[0] > 1.833 && x[0] < 1.899) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}



void odecet_fit() {

    TFile* d1 = new TFile("pp500_D0.root");


// Vizualizace fitu - kubický -----------------------------------------------------------------------------------------------------------------------------
    TH1F* h1 = (TH1F*) d1 -> Get("D0 signal;1");
    TCanvas *c1 = new TCanvas("c1", "Invariant mass", 1400, 1000);

    c1->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    
    h1->Draw("PE");
    c1->Update();

    TF1 *fc = new TF1("fc",fcubic,1.7,2.05,4);
    reject = kTRUE;
    h1->Fit(fc,"0");

    Double_t AAAA=0,BBBB=0,CCCC=0,DDDD=0;

    AAAA = fc->GetParameter(0);
    BBBB = fc->GetParameter(1); //parametr u x
    CCCC = fc->GetParameter(2); //parametr u x^2
    DDDD = fc->GetParameter(3); //parametr u x^3

    fc->Draw("same");
    fc->SetLineWidth(0.5);
    fc->SetLineColor(kRed);

    h1->SetMarkerStyle(8);
    h1->SetMarkerColor(kBlue);
    h1->SetLineWidth(0);
    h1->SetTitle(" ; ");
    h1->SetStats(0);
    h1->SetLineColor(kBlue);
    TLegend *l1 = new TLegend(0.55,0.65, 0.8, 0.85,"","brNDC");
    l1->AddEntry(h1, "Unlike-sign K^{#pm}#pi^{#mp} pairs", "pl");
    l1->AddEntry(fc, "Third order polynomial", "pl");
    l1->AddEntry((TObject*)0, "", "");
    l1->AddEntry((TObject*)0, "THIS THESIS", "");
    l1->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l1->SetFillStyle(0);
    l1->SetLineColor(0);
    l1->SetTextSize(0.04);
    l1->Draw("same");
    h1->GetYaxis()->SetTitle("Counts [-]");
    h1->GetYaxis()->SetTitleOffset(1.);
    h1->GetYaxis()->SetLabelSize(0.03);
    h1->GetYaxis()->SetTitleSize(0.04);
    h1->GetYaxis()->SetLabelFont(42);
    h1->GetYaxis()->SetTitleFont(42);
    h1->GetYaxis()->CenterTitle(kTRUE);
    h1->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    h1->GetXaxis()->SetLabelFont(42);
    h1->GetXaxis()->SetTitleFont(42);
    h1->GetXaxis()->SetTitleSize(0.04);
    h1->GetXaxis()->SetLabelSize(0.03);
    h1->GetXaxis()->CenterTitle(kTRUE);
    c1->SaveAs(Form("pp500_fit_cubic.pdf"));
    c1->Close();



// Fit Gauss -----------------------------------------------------------------------------------------------------------------------------
/*    TCanvas *c6 = new TCanvas("c6", "Invariant mass", 1400, 1000);



    TH1F *h6 = new TH1F("h6","h6",32, 1.7, 2.05);


    TF1 *gauss = new TF1("gauss","[0]*exp(-(x-[1])^2/(2*[2]^2))",1.7,2.05);
    gauss->SetParameters(1000,1.85,0.05);

    for(int m=0; m<=32;m++){
        RR=h1->GetBinContent(m);
        UU=h1->GetBinError(m);

        SS=h5->GetBinCenter(m);
        TT=bckg_funct_cubic->TF1::Eval(SS);

        h6->SetBinContent(m,RR-TT);
        h6->SetBinError(m,UU);

    }



    h6->Draw("PE");
    h6->Fit("gauss");
    TLine *line1 = new TLine(1.8,7000,1.8,-4000);
    line1->Draw("same");
    TLine *line2 = new TLine(1.898,7000,1.898,-4000);
    line2->Draw("same");
    c6->Update();

    gauss->Draw("same");


    h6->SetMarkerStyle(8);
    h6->SetMarkerColor(kBlue);
    h6->SetLineWidth(0);
    h6->SetTitle(" ; ");
    h6->SetStats(0);
    h6->SetLineColor(kBlue);
    TLegend *l6 = new TLegend(0.55,0.69, 0.8, 0.907,"","brNDC");
    l6->AddEntry(h6, "Unlike-sign pairs", "pl");
    l6->SetFillStyle(0);
    l6->SetLineColor(0);
    l6->SetTextSize(0.04);
    l6->Draw("same");
    h6->GetYaxis()->SetTitle("Counts [-]");
    h6->GetYaxis()->SetTitleOffset(1.);
    h6->GetYaxis()->SetLabelSize(0.03);
    h6->GetYaxis()->SetTitleSize(0.04);
    h6->GetYaxis()->SetLabelFont(42);
    h6->GetYaxis()->SetTitleFont(42);
    h6->GetYaxis()->CenterTitle(kTRUE);
    h6->GetXaxis()->SetTitle("m_{inv} [GeV/#it{c}^{2}]");
    h6->GetXaxis()->SetLabelFont(42);
    h6->GetXaxis()->SetTitleFont(42);
    h6->GetXaxis()->SetTitleSize(0.04);
    h6->GetXaxis()->SetLabelSize(0.03);
    h6->GetXaxis()->CenterTitle(kTRUE);
    c6->SaveAs(Form("Gauss_fit.pdf"));
    c6->Close();*/











// Fit cubic-----------------------------------------------------------------------------------------------------------------------------
    TH1F* h2 = (TH1F*) d1 -> Get("odecet_cubic;1");
    TCanvas *c2 = new TCanvas("c2", "Odecet kub. funkce", 1400, 1000);
    h2->Rebin();
    h2->Draw("PE");

    c2->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);

    c2->Update();
    h2->SetMarkerStyle(8);
    h2->SetMarkerColor(kBlue);
    h2->SetLineWidth(0);
    h2->SetTitle(" ; ");
    h2->SetStats(0);
    h2->SetLineColor(kBlue);
    TLegend *l2 = new TLegend(0.45,0.69, 0.7, 0.907,"","brNDC");
    l2->AddEntry(h2, "(Unlike-Sign K^{#pm}#pi^{#mp})-Fit", "pl");
    l2->AddEntry((TObject*)0, "", "");
    l2->AddEntry((TObject*)0, "THIS THESIS", "");
    l2->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l2->SetFillStyle(0);
    l2->SetLineColor(0);
    l2->SetTextSize(0.04);
    l2->Draw("same");
    h2->GetYaxis()->SetTitle("Counts [-]");
    h2->GetYaxis()->SetTitleOffset(1.);
    h2->GetYaxis()->SetLabelSize(0.03);
    h2->GetYaxis()->SetTitleSize(0.04);
    h2->GetYaxis()->SetLabelFont(42);
    h2->GetYaxis()->SetTitleFont(42);
    h2->GetYaxis()->CenterTitle(kTRUE);
    h2->GetYaxis()->SetRangeUser(-3000,5000);
    h2->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    h2->GetXaxis()->SetLabelFont(42);
    h2->GetXaxis()->SetTitleFont(42);
    h2->GetXaxis()->SetTitleSize(0.04);
    h2->GetXaxis()->SetLabelSize(0.03);
    h2->GetXaxis()->CenterTitle(kTRUE);
    c2->SaveAs(Form("pp500_Cubic_SideBand.pdf"));
    c2->Close();



}
