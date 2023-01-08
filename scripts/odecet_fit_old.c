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
    if (reject && x[0] > 1.761 && x[0] < 1.935) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0];
}

Double_t fquadratic(Double_t *x, Double_t *par)
{
    if (reject && x[0] > 1.761 && x[0] < 1.935) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t fcubic(Double_t *x, Double_t *par)
{
    if (reject && x[0] > 1.761 && x[0] < 1.935) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
}



void odecet_fit() {

    TFile* d1 = new TFile("D0.root");

// Vizualizace fitu - lineární -----------------------------------------------------------------------------------------------------------------------------
    TH1F* h1 = (TH1F*) d1 -> Get("signal;1");
    TCanvas *c1 = new TCanvas("c1", "Invariant mass", 1400, 1100);
    h1->Draw("PE");
    c1->Update();

    TF1 *fl = new TF1("fl",fline,1.7,2.05,2);
    fl->SetParameters(2,-1);
    reject = kTRUE;
    h1->Fit(fl,"0");

    Double_t AA=0,BB=0;

    AA = fl->GetParameter(1);
    BB = fl->GetParameter(0);
    fl->Draw("same");

    h1->SetMarkerStyle(8);
    h1->SetMarkerColor(kBlue);
    h1->SetLineWidth(0);
    h1->SetTitle(" ; ");
    h1->SetStats(0);
    h1->SetLineColor(kBlue);
    TLegend *l1 = new TLegend(0.55,0.80, 0.8, 0.907,"","brNDC");
    l1->AddEntry(h1, "Unlike-sign pairs", "pl");
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
    h1->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h1->GetXaxis()->SetLabelFont(42);
    h1->GetXaxis()->SetTitleFont(42);
    h1->GetXaxis()->SetTitleSize(0.04);
    h1->GetXaxis()->SetLabelSize(0.03);
    h1->GetXaxis()->CenterTitle(kTRUE);
    c1->SaveAs(Form("Fit_linearni.pdf"));
    c1->Close();

// Vizualizace fitu - kvadratický -----------------------------------------------------------------------------------------------------------------------------
    TH1F* h2 = (TH1F*) d1 -> Get("signal;1");
    TCanvas *c2 = new TCanvas("c2", "Invariant mass", 1400, 1100);

    h2->Draw("PE");
    c2->Update();

    TF1 *fq = new TF1("fq",fquadratic,1.7,2.05,3);
    reject = kTRUE;
    h2->Fit(fq,"0");

    Double_t AAA=0,BBB=0,CCC=0;

    AAA = fq->GetParameter(0);
    BBB = fq->GetParameter(1); //parametr u x
    CCC = fq->GetParameter(2); //parametr u x^2

    fq->Draw("same");

    h2->SetMarkerStyle(8);
    h2->SetMarkerColor(kBlue);
    h2->SetLineWidth(0);
    h2->SetTitle(" ; ");
    h2->SetStats(0);
    h2->SetLineColor(kBlue);
    TLegend *l2 = new TLegend(0.55,0.80, 0.8, 0.907,"","brNDC");
    l2->AddEntry(h2, "Unlike-sign pairs", "pl");
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
    h2->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h2->GetXaxis()->SetLabelFont(42);
    h2->GetXaxis()->SetTitleFont(42);
    h2->GetXaxis()->SetTitleSize(0.04);
    h2->GetXaxis()->SetLabelSize(0.03);
    h2->GetXaxis()->CenterTitle(kTRUE);
    c2->SaveAs(Form("Fit_kvadraticky.pdf"));
    c2->Close();

// Vizualizace fitu - kubický -----------------------------------------------------------------------------------------------------------------------------
    TH1F* h4 = (TH1F*) d1 -> Get("signal;1");
    TCanvas *c4 = new TCanvas("c4", "Invariant mass", 1400, 1100);

    h4->Draw("PE");
    c4->Update();

    TF1 *fc = new TF1("fc",fcubic,1.7,2.05,4);
    reject = kTRUE;
    h4->Fit(fc,"0");

    Double_t AAAA=0,BBBB=0,CCCC=0,DDDD=0;

    AAAA = fc->GetParameter(0);
    BBBB = fc->GetParameter(1); //parametr u x
    CCCC = fc->GetParameter(2); //parametr u x^2
    DDDD = fc->GetParameter(3); //parametr u x^3

    fc->Draw("same");

    h4->SetMarkerStyle(8);
    h4->SetMarkerColor(kBlue);
    h4->SetLineWidth(0);
    h4->SetTitle(" ; ");
    h4->SetStats(0);
    h4->SetLineColor(kBlue);
    TLegend *l4 = new TLegend(0.55,0.80, 0.8, 0.907,"","brNDC");
    l4->AddEntry(h4, "Unlike-sign pairs", "pl");
    l4->SetFillStyle(0);
    l4->SetLineColor(0);
    l4->SetTextSize(0.04);
    l4->Draw("same");
    h4->GetYaxis()->SetTitle("Counts [-]");
    h4->GetYaxis()->SetTitleOffset(1.);
    h4->GetYaxis()->SetLabelSize(0.03);
    h4->GetYaxis()->SetTitleSize(0.04);
    h4->GetYaxis()->SetLabelFont(42);
    h4->GetYaxis()->SetTitleFont(42);
    h4->GetYaxis()->CenterTitle(kTRUE);
    h4->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h4->GetXaxis()->SetLabelFont(42);
    h4->GetXaxis()->SetTitleFont(42);
    h4->GetXaxis()->SetTitleSize(0.04);
    h4->GetXaxis()->SetLabelSize(0.03);
    h4->GetXaxis()->CenterTitle(kTRUE);
    c4->SaveAs(Form("Fit_kubicky.pdf"));
    c4->Close();


// Odecet kvadraticky -----------------------------------------------------------------------------------------------------------------------------
    TCanvas *c3 = new TCanvas("c3", "Invariant mass", 1400, 1100);

    TF1 *bckg_funct_quadratic = new TF1("bckg_funct_quadratic","[0]+[1]*x+[2]*x*x",1.7,2.05);
    bckg_funct_quadratic->SetParameters(AAA,BBB,CCC);

    TH1F *h3 = new TH1F("h3","h3",32, 1.7, 2.05);

    Double_t NN=0, OO=0, PP=0, QQ=0;

    for(int l=0; l<=32;l++){
        NN=h1->GetBinContent(l);
        QQ=h1->GetBinError(l);

        OO=h3->GetBinCenter(l);
        PP=bckg_funct_quadratic->TF1::Eval(OO);

        h3->SetBinContent(l,NN-PP);
        h3->SetBinError(l,QQ);

    }

    h3->Draw("PE");
    c3->Update();

    h3->SetMarkerStyle(8);
    h3->SetMarkerColor(kBlue);
    h3->SetLineWidth(0);
    h3->SetTitle(" ; ");
    h3->SetStats(0);
    h3->SetLineColor(kBlue);
    TLegend *l3 = new TLegend(0.55,0.80, 0.8, 0.907,"","brNDC");
    l3->AddEntry(h3, "Unlike-sign pairs", "pl");
    l3->SetFillStyle(0);
    l3->SetLineColor(0);
    l3->SetTextSize(0.04);
    l3->Draw("same");
    h3->GetYaxis()->SetTitle("Counts [-]");
    h3->GetYaxis()->SetTitleOffset(1.);
    h3->GetYaxis()->SetLabelSize(0.03);
    h3->GetYaxis()->SetTitleSize(0.04);
    h3->GetYaxis()->SetLabelFont(42);
    h3->GetYaxis()->SetTitleFont(42);
    h3->GetYaxis()->CenterTitle(kTRUE);
    h3->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h3->GetXaxis()->SetLabelFont(42);
    h3->GetXaxis()->SetTitleFont(42);
    h3->GetXaxis()->SetTitleSize(0.04);
    h3->GetXaxis()->SetLabelSize(0.03);
    h3->GetXaxis()->CenterTitle(kTRUE);
    c3->SaveAs(Form("Odecet_funkce_kvadraticky.pdf"));
    c3->Close();



// Odecet_kubicky -----------------------------------------------------------------------------------------------------------------------------
    TCanvas *c5 = new TCanvas("c5", "Invariant mass", 1400, 1100);


    TF1 *bckg_funct_cubic = new TF1("bckg_funct_cubic","[0]+[1]*x+[2]*x*x+[3]*x*x*x",1.7,2.05);
    bckg_funct_cubic->SetParameters(AAAA,BBBB,CCCC,DDDD);

    TH1F *h5 = new TH1F("h5","h5",32, 1.7, 2.05);

    Double_t RR=0, SS=0, TT=0, UU=0;

    for(int m=0; m<=32;m++){
        RR=h1->GetBinContent(m);
        UU=h1->GetBinError(m);

        SS=h5->GetBinCenter(m);
        TT=bckg_funct_cubic->TF1::Eval(SS);

        h5->SetBinContent(m,RR-TT);
        h5->SetBinError(m,UU);

    }

    h5->Draw("PE");
    c5->Update();

    h5->SetMarkerStyle(8);
    h5->SetMarkerColor(kBlue);
    h5->SetLineWidth(0);
    h5->SetTitle(" ; ");
    h5->SetStats(0);
    h5->SetLineColor(kBlue);
    TLegend *l5 = new TLegend(0.55,0.80, 0.8, 0.907,"","brNDC");
    l5->AddEntry(h5, "Unlike-sign pairs", "pl");
    l5->SetFillStyle(0);
    l5->SetLineColor(0);
    l5->SetTextSize(0.04);
    l5->Draw("same");
    h5->GetYaxis()->SetTitle("Counts [-]");
    h5->GetYaxis()->SetTitleOffset(1.);
    h5->GetYaxis()->SetLabelSize(0.03);
    h5->GetYaxis()->SetTitleSize(0.04);
    h5->GetYaxis()->SetLabelFont(42);
    h5->GetYaxis()->SetTitleFont(42);
    h5->GetYaxis()->CenterTitle(kTRUE);
    h5->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h5->GetXaxis()->SetLabelFont(42);
    h5->GetXaxis()->SetTitleFont(42);
    h5->GetXaxis()->SetTitleSize(0.04);
    h5->GetXaxis()->SetLabelSize(0.03);
    h5->GetXaxis()->CenterTitle(kTRUE);
    c5->SaveAs(Form("Odecet_funkce_kubicky.pdf"));
    c5->Close();



// Fit Gauss -----------------------------------------------------------------------------------------------------------------------------
    TCanvas *c6 = new TCanvas("c6", "Invariant mass", 1400, 1100);



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
    TLegend *l6 = new TLegend(0.55,0.80, 0.8, 0.907,"","brNDC");
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
    h6->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h6->GetXaxis()->SetLabelFont(42);
    h6->GetXaxis()->SetTitleFont(42);
    h6->GetXaxis()->SetTitleSize(0.04);
    h6->GetXaxis()->SetLabelSize(0.03);
    h6->GetXaxis()->CenterTitle(kTRUE);
    c6->SaveAs(Form("Gauss_fit.pdf"));
    c6->Close();









// Fit lineární-----------------------------------------------------------------------------------------------------------------------------
    TH1F* h10 = (TH1F*) d1 -> Get("odecet_linear;1");
    TCanvas *c10 = new TCanvas("c10", "Odecet lin. funkce", 1400, 1100);
    h10->Draw("PE");
    c10->Update();
    h10->SetMarkerStyle(8);
    h10->SetMarkerColor(kBlue);
    h10->SetLineWidth(0);
    h10->SetTitle(" ; ");
    h10->SetStats(0);
    h10->SetLineColor(kBlue);
    TLegend *l10 = new TLegend(0.65,0.80, 0.8, 0.907,"","brNDC");
    l10->AddEntry(h10, "US-L. Fit", "pl");
    l10->SetFillStyle(0);
    l10->SetLineColor(0);
    l10->SetTextSize(0.04);
    l10->Draw("same");
    h10->GetYaxis()->SetTitle("Counts [-]");
    h10->GetYaxis()->SetTitleOffset(1.);
    h10->GetYaxis()->SetLabelSize(0.03);
    h10->GetYaxis()->SetTitleSize(0.04);
    h10->GetYaxis()->SetLabelFont(42);
    h10->GetYaxis()->SetTitleFont(42);
    h10->GetYaxis()->CenterTitle(kTRUE);
    h10->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h10->GetXaxis()->SetLabelFont(42);
    h10->GetXaxis()->SetTitleFont(42);
    h10->GetXaxis()->SetTitleSize(0.04);
    h10->GetXaxis()->SetLabelSize(0.03);
    h10->GetXaxis()->CenterTitle(kTRUE);
    c10->SaveAs(Form("Odecet_funkce_linearni.pdf"));
    c10->Close();

// Fit kvadratický-----------------------------------------------------------------------------------------------------------------------------
    TH1F* h11 = (TH1F*) d1 -> Get("odecet_quadratic;1");
    TCanvas *c11 = new TCanvas("c11", "Odecet kv. funkce", 1400, 1100);
    h11->Draw("PE");
    c11->Update();
    h11->SetMarkerStyle(8);
    h11->SetMarkerColor(kBlue);
    h11->SetLineWidth(0);
    h11->SetTitle(" ; ");
    h11->SetStats(0);
    h11->SetLineColor(kBlue);
    TLegend *l11 = new TLegend(0.65,0.80, 0.8, 0.907,"","brNDC");
    l11->AddEntry(h11, "US-Q. Fit", "pl");
    l11->SetFillStyle(0);
    l11->SetLineColor(0);
    l11->SetTextSize(0.04);
    l11->Draw("same");
    h11->GetYaxis()->SetTitle("Counts [-]");
    h11->GetYaxis()->SetTitleOffset(1.);
    h11->GetYaxis()->SetLabelSize(0.03);
    h11->GetYaxis()->SetTitleSize(0.04);
    h11->GetYaxis()->SetLabelFont(42);
    h11->GetYaxis()->SetTitleFont(42);
    h11->GetYaxis()->CenterTitle(kTRUE);
    h11->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h11->GetXaxis()->SetLabelFont(42);
    h11->GetXaxis()->SetTitleFont(42);
    h11->GetXaxis()->SetTitleSize(0.04);
    h11->GetXaxis()->SetLabelSize(0.03);
    h11->GetXaxis()->CenterTitle(kTRUE);
    c11->SaveAs(Form("Odecet_funkce_kvadraticka.pdf"));
    c11->Close();



}
