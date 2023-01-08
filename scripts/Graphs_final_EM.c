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
#include "TStyle.h"
#include "TLegend.h"

using namespace std;

void Graphs_final_EM() {

    TFile* d1 = new TFile("pp500_D0.root");
    TFile* d2 = new TFile("pp500_Dstar.root");
    TFile* d3 = new TFile("pp500_D0_EM.root");


// Signal narrow window all p_T -----------------------------------------------------------------------------------------------------------------------------
    TH1F* s1 = (TH1F*) d1 -> Get("D0 signal;1");
    TH1F* back1 = (TH1F*) d3 -> Get("D0 signal;1");
//    TH1F* diff1 = (TH1F*) d3 -> Get("D0 rozdil;1");


    TCanvas *c1 = new TCanvas("c1", "Invariant mass", 1400, 1000);
    //  gPad->SetMargin(2.5,0.1,0.03,0.03);
    // pionpt->SetLogy();

    double si1= s1->Integral();
    double bi1= back1->Integral();

    double norm1= si1/bi1;

    TH1F *b1 = new TH1F();
    *b1 = norm1*(*back1);

    TH1F *roz1 = new TH1F("Rall","Rall",32, 1.7, 2.05);
    roz1->Add(s1,b1,1,-1);

    TH1F *r1 = new TH1F();
    *r1 = 10*(*roz1);

    s1->Draw("PE");
    b1->Draw("PE same");
    r1->Draw("PE same");


    c1->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c1->Update();

    s1->SetMarkerStyle(22);
    b1->SetMarkerStyle(23);
    r1->SetMarkerStyle(33);

    s1->SetMarkerSize(1.2);
    b1->SetMarkerSize(1.2);
    r1->SetMarkerSize(1.8);

    s1->SetMarkerColor(kBlue);
    b1->SetMarkerColor(kRed);
    r1->SetMarkerColor(kBlack);
    
    s1->SetTitle(" ; ");
    s1->SetStats(0);

    /* TText *text = new TText(1.291488,746304.2,"THIS THESIS");
     text->SetTextAlign(23);
     text->SetTextSize(.04);
     text->Draw();*/
    TLegend *l1 = new TLegend(0.55,0.69, 0.73, 0.89,"","brNDC");
    l1->AddEntry(s1, "Unlike-sign K^{#pm}#pi^{#mp} pairs", "pl");
    l1->AddEntry(b1, "Mixing K^{#pm}#pi^{#mp} pairs (background)", "pl");
    l1->AddEntry(r1, "Difference US-Mixing scaled by 10", "pl");
    l1->AddEntry((TObject*)0, "", "");
    l1->AddEntry((TObject*)0, "THIS THESIS", "");
    l1->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");

    l1->SetFillStyle(0);
    l1->SetLineColor(0);
    l1->SetTextSize(0.03);
    l1->Draw("same");

    s1->GetYaxis()->SetTitle("Counts [-]");
    s1->GetYaxis()->SetTitleOffset(1.);
    s1->GetYaxis()->SetLabelSize(0.03);
    s1->GetYaxis()->SetTitleSize(0.04);
    s1->GetYaxis()->SetLabelFont(42);
    s1->GetYaxis()->SetTitleFont(42);
    s1->GetYaxis()->SetRangeUser(-250000,2100000);
    s1->GetYaxis()->CenterTitle(kTRUE);
    s1->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s1->GetXaxis()->SetLabelFont(42);
    s1->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s1->GetXaxis()->SetTitleSize(0.04);
    s1->GetXaxis()->SetLabelSize(0.03);
    s1->GetXaxis()->CenterTitle(kTRUE);
    c1->SaveAs(Form("pp500_EM_Mass_narrow_all_pT.pdf"));
    c1->Close();

// Signal narrow window p_T 1-2 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s2 = (TH1F*) d1 -> Get("D0 S12;1");
    TH1F* back2 = (TH1F*) d3 -> Get("D0 S12;1");
//    TH1F* diff2 = (TH1F*) d3 -> Get("D0 rozdil12;1");

    TCanvas *c2 = new TCanvas("c2", "Invariant mass", 1400, 1000);

    double si2= s2->Integral();
    double bi2= back2->Integral();

    double norm2= si2/bi2;

    TH1F *b2 = new TH1F();
    *b2 = norm2*(*back2);

    TH1F *roz2 = new TH1F("R12","R12",32, 1.7, 2.05);
    roz2->Add(s2,b2,1,-1);

    TH1F *r2 = new TH1F();
    *r2 = 10*(*roz2);

    s2->Draw("PE");
    b2->Draw("PE same");
    r2->Draw("PE same");


    c2->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c2->Update();

    s2->SetMarkerStyle(22);
    b2->SetMarkerStyle(23);
    r2->SetMarkerStyle(33);

    s2->SetMarkerSize(1.2);
    b2->SetMarkerSize(1.2);
    r2->SetMarkerSize(1.8);

    s2->SetMarkerColor(kBlue);
    b2->SetMarkerColor(kRed);
    r2->SetMarkerColor(kBlack);

    s2->SetTitle(" ; ");
    s2->SetStats(0);

    TLegend *l2 = new TLegend(0.55,0.65, 0.73, 0.89,"","brNDC");
    l2->AddEntry(s2, "Unlike-sign K^{#pm}#pi^{#mp} pairs", "pl");
    l2->AddEntry(b2, "Mixing K^{#pm}#pi^{#mp} pairs (background)", "pl");
    l2->AddEntry(r2, "Difference US-Mixing scaled by 10", "pl");
    l2->AddEntry((TObject*)0, "", "");
    l2->AddEntry((TObject*)0, "THIS THESIS", "");
    l2->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l2->AddEntry((TObject*)0, "1 < #it{p}_{T}(K#pi) < 2 GeV/#it{c}^{2}", "");

    l2->SetFillStyle(0);
    l2->SetLineColor(0);
    l2->SetTextSize(0.03);
    l2->Draw("same");

    s2->GetYaxis()->SetTitle("Counts [-]");
    s2->GetYaxis()->SetTitleOffset(1.);
    s2->GetYaxis()->SetLabelSize(0.03);
    s2->GetYaxis()->SetTitleSize(0.04);
    s2->GetYaxis()->SetLabelFont(42);
    s2->GetYaxis()->SetTitleFont(42);
    s2->GetYaxis()->SetRangeUser(-100000,850000);
    s2->GetYaxis()->CenterTitle(kTRUE);
    s2->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s2->GetXaxis()->SetLabelFont(42);
    s2->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s2->GetXaxis()->SetTitleSize(0.04);
    s2->GetXaxis()->SetLabelSize(0.03);
    s2->GetXaxis()->CenterTitle(kTRUE);
    c2->SaveAs(Form("pp500_EM_Mass_narrow_12_pT.pdf"));
    c2->Close();

// Signal narrow window p_T 2-3 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s3 = (TH1F*) d1 -> Get("D0 S23;1");
    TH1F* back3 = (TH1F*) d3 -> Get("D0 S23;1");
//    TH1F* diff3 = (TH1F*) d3 -> Get("D0 rozdil23;1");

    TCanvas *c3 = new TCanvas("c3", "Invariant mass", 1400, 1000);

    double si3= s3->Integral();
    double bi3= back3->Integral();

    double norm3= si3/bi3;

    TH1F *b3 = new TH1F();
    *b3 = norm3*(*back3);

    TH1F *roz3 = new TH1F("R23","R23",32, 1.7, 2.05);
    roz3->Add(s3,b3,1,-1);

    TH1F *r3 = new TH1F();
    *r3 = 10*(*roz3);

    s3->Draw("PE");
    b3->Draw("PE same");
    r3->Draw("PE same");
    c3->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c3->Update();

    s3->SetMarkerStyle(22);
    b3->SetMarkerStyle(23);
    r3->SetMarkerStyle(33);

    s3->SetMarkerSize(1.2);
    b3->SetMarkerSize(1.2);
    r3->SetMarkerSize(1.8);

    s3->SetMarkerColor(kBlue);
    b3->SetMarkerColor(kRed);
    r3->SetMarkerColor(kBlack);

    s3->SetTitle(" ; ");
    s3->SetStats(0);

    TLegend *l3 = new TLegend(0.55,0.65, 0.73, 0.89,"","brNDC");
    l3->AddEntry(s3, "Unlike-sign K^{#pm}#pi^{#mp} pairs", "pl");
    l3->AddEntry(b3, "Mixing K^{#pm}#pi^{#mp} pairs (background)", "pl");
    l3->AddEntry(r3, "Difference US-Mixing scaled by 10", "pl");
    l3->AddEntry((TObject*)0, "", "");
    l3->AddEntry((TObject*)0, "THIS THESIS", "");
    l3->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l3->AddEntry((TObject*)0, "2 < #it{p}_{T}(K#pi) < 3 GeV/#it{c}^{2}", "");

    l3->SetFillStyle(0);
    l3->SetLineColor(0);
    l3->SetTextSize(0.03);
    l3->Draw("same");

    s3->GetYaxis()->SetTitle("Counts [-]");
    s3->GetYaxis()->SetTitleOffset(1.);
    s3->GetYaxis()->SetLabelSize(0.03);
    s3->GetYaxis()->SetTitleSize(0.04);
    s3->GetYaxis()->SetLabelFont(42);
    s3->GetYaxis()->SetTitleFont(42);
    s3->GetYaxis()->SetRangeUser(-50000,300000);
    s3->GetYaxis()->CenterTitle(kTRUE);
    s3->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s3->GetXaxis()->SetLabelFont(42);
    s3->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s3->GetXaxis()->SetTitleSize(0.04);
    s3->GetXaxis()->SetLabelSize(0.03);
    s3->GetXaxis()->CenterTitle(kTRUE);
    c3->SaveAs(Form("pp500_EM_Mass_narrow_23_pT.pdf"));
    c3->Close();

// Signal narrow window p_T 3-4 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s4 = (TH1F*) d1 -> Get("D0 S34;1");
    TH1F* back4 = (TH1F*) d3 -> Get("D0 S34;1");
//    TH1F* diff4 = (TH1F*) d3 -> Get("D0 rozdil34;1");

    TCanvas *c4 = new TCanvas("c4", "Invariant mass", 1400, 1000);

    double si4= s4->Integral();
    double bi4= back4->Integral();

    double norm4= si4/bi4;

    TH1F *b4 = new TH1F();
    *b4 = norm4*(*back4);

    TH1F *roz4 = new TH1F("R34","R34",32, 1.7, 2.05);
    roz4->Add(s4,b4,1,-1);

    TH1F *r4 = new TH1F();
    *r4 = 10*(*roz4);

    s4->Draw("PE");
    b4->Draw("PE same");
    r4->Draw("PE same");
    c4->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c4->Update();

    s4->SetMarkerStyle(22);
    b4->SetMarkerStyle(23);
    r4->SetMarkerStyle(33);

    s4->SetMarkerSize(1.2);
    b4->SetMarkerSize(1.2);
    r4->SetMarkerSize(1.8);

    s4->SetMarkerColor(kBlue);
    b4->SetMarkerColor(kRed);
    r4->SetMarkerColor(kBlack);

    s4->SetTitle(" ; ");
    s4->SetStats(0);

    TLegend *l4 = new TLegend(0.55,0.65, 0.73, 0.89,"","brNDC");
    l4->AddEntry(s4, "Unlike-sign K^{#pm}#pi^{#mp} pairs", "pl");
    l4->AddEntry(b4, "Mixing K^{#pm}#pi^{#mp} pairs (background)", "pl");
    l4->AddEntry(r4, "Difference US-Mixing scaled by 10", "pl");
    l4->AddEntry((TObject*)0, "", "");
    l4->AddEntry((TObject*)0, "THIS THESIS", "");
    l4->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l4->AddEntry((TObject*)0, "3 < #it{p}_{T}(K#pi) < 4 GeV/#it{c}^{2}", "");

    l4->SetFillStyle(0);
    l4->SetLineColor(0);
    l4->SetTextSize(0.03);
    l4->Draw("same");

    s4->GetYaxis()->SetTitle("Counts [-]");
    s4->GetYaxis()->SetTitleOffset(1.2);
    s4->GetYaxis()->SetLabelSize(0.03);
    s4->GetYaxis()->SetTitleSize(0.04);
    s4->GetYaxis()->SetLabelFont(42);
    s4->GetYaxis()->SetTitleFont(42);
    s4->GetYaxis()->SetRangeUser(-10000,60000);
    s4->GetYaxis()->CenterTitle(kTRUE);
    s4->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s4->GetXaxis()->SetLabelFont(42);
    s4->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s4->GetXaxis()->SetTitleSize(0.04);
    s4->GetXaxis()->SetLabelSize(0.03);
    s4->GetXaxis()->CenterTitle(kTRUE);
    c4->SaveAs(Form("pp500_EM_Mass_narrow_34_pT.pdf"));
    c4->Close();

// Signal narrow window p_T >4 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s5 = (TH1F*) d1 -> Get("D0 S45;1");
    TH1F* back5 = (TH1F*) d3 -> Get("D0 S45;1");
//    TH1F* diff5 = (TH1F*) d3 -> Get("D0 rozdil45;1");

    TCanvas *c5 = new TCanvas("c5", "Invariant mass", 1400, 1000);

    double si5= s5->Integral();
    double bi5= back5->Integral();

    double norm5= si5/bi5;

    TH1F *b5 = new TH1F();
    *b5 = norm5*(*back5);


    TH1F *roz5 = new TH1F("R45","R45",32, 1.7, 2.05);
    roz5->Add(s5,b5,1,-1);

    TH1F *r5 = new TH1F();
    *r5 = 10*(*roz5);

    s5->Draw("PE");
    b5->Draw("PE same");
    r5->Draw("PE same");
    c5->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c5->Update();

    s5->SetMarkerStyle(22);
    b5->SetMarkerStyle(23);
    r5->SetMarkerStyle(33);

    s5->SetMarkerSize(1.2);
    b5->SetMarkerSize(1.2);
    r5->SetMarkerSize(1.8);

    s5->SetMarkerColor(kBlue);
    b5->SetMarkerColor(kRed);
    r5->SetMarkerColor(kBlack);

    s5->SetTitle(" ; ");
    s5->SetStats(0);

    TLegend *l5 = new TLegend(0.55,0.65, 0.73, 0.89,"","brNDC");
    l5->AddEntry(s5, "Unlike-sign K^{#pm}#pi^{#mp} pairs", "pl");
    l5->AddEntry(b5, "Mixing K^{#pm}#pi^{#mp} pairs (background)", "pl");
    l5->AddEntry(r5, "Difference US-Mixing scaled by 10", "pl");
    l5->AddEntry((TObject*)0, "", "");
    l5->AddEntry((TObject*)0, "THIS THESIS", "");
    l5->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l5->AddEntry((TObject*)0, "#it{p}_{T}(K#pi) > 4 GeV/#it{c}^{2}", "");

    l5->SetFillStyle(0);
    l5->SetLineColor(0);
    l5->SetTextSize(0.03);
    l5->Draw("same");

    s5->GetYaxis()->SetTitle("Counts [-]");
    s5->GetYaxis()->SetTitleOffset(1.);
    s5->GetYaxis()->SetLabelSize(0.03);
    s5->GetYaxis()->SetTitleSize(0.04);
    s5->GetYaxis()->SetLabelFont(42);
    s5->GetYaxis()->SetTitleFont(42);
    s5->GetYaxis()->SetRangeUser(-5000,20000);
    s5->GetYaxis()->CenterTitle(kTRUE);
    s5->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s5->GetXaxis()->SetLabelFont(42);
    s5->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s5->GetXaxis()->SetTitleSize(0.04);
    s5->GetXaxis()->SetLabelSize(0.03);
    s5->GetXaxis()->CenterTitle(kTRUE);
    c5->SaveAs(Form("pp500_EM_Mass_narrow_45_pT.pdf"));
    c5->Close();

// Signal wide window all p_T -----------------------------------------------------------------------------------------------------------------------------
    TH1F* s6 = (TH1F*) d1 -> Get("D0 signalw;1");
    TH1F* b6 = (TH1F*) d1 -> Get("D0 backgroundw;1");
    TH1F* diff6 = (TH1F*) d1 -> Get("D0 rozdilw;1");
    TH1F* e6 = (TH1F*) d3 -> Get("D0 signalw;1");

    TCanvas *c6 = new TCanvas("c6", "Invariant mass", 1400, 1000);

    TH1F *r6 = new TH1F();
    *r6 = 1*(*diff6);

    s6->Scale(1/s6->GetEntries());
    b6->Scale(1/b6->GetEntries());
    e6->Scale(1/e6->GetEntries());


    s6->Draw("PE");
    b6->Draw("PE same");
    r6->Draw("PE same");
    e6->Draw("HIST same C");

    e6->SetLineColor(kGreen);
    e6->SetLineWidth(1.5);

    c6->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c6->Update();

    s6->SetMarkerStyle(22);
    b6->SetMarkerStyle(23);
    r6->SetMarkerStyle(33);

    s6->SetMarkerSize(0.7);
    b6->SetMarkerSize(0.7);
    r6->SetMarkerSize(0.9);

    s6->SetMarkerColor(kBlue);
    b6->SetMarkerColor(kRed);
    r6->SetMarkerColor(kBlack);

    s6->SetTitle(" ; ");
    s6->SetStats(0);

    TLegend *l6 = new TLegend(0.55,0.65, 0.73, 0.89,"","brNDC");
    l6->AddEntry(s6, "Unlike-sign K^{#pm}#pi^{#mp} pairs", "pl");
    l6->AddEntry(b6, "Like-sign K^{#pm}#pi^{#pm} pairs (background)", "pl");
//    l6->AddEntry(r6, "Difference US-LS", "pl");
    l6->AddEntry(e6, "Event-Mixing (background)", "pl");
    l6->AddEntry((TObject*)0, "", "");
    l6->AddEntry((TObject*)0, "THIS THESIS", "");
    l6->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
//    l6->AddEntry((TObject*)0, "#it{p}_{T}(K#pi) > 4 GeV/#it{c}^{2}", "");

    l6->SetFillStyle(0);
    l6->SetLineColor(0);
    l6->SetTextSize(0.03);
    l6->Draw("same");

    s6->GetYaxis()->SetTitle("Counts [-]");
    s6->GetYaxis()->SetTitleOffset(1.);
    s6->GetYaxis()->SetLabelSize(0.03);
    s6->GetYaxis()->SetTitleSize(0.04);
    s6->GetYaxis()->SetLabelFont(42);
    s6->GetYaxis()->SetTitleFont(42);
 //   s6->GetYaxis()->SetRangeUser(0,370000);
    s6->GetYaxis()->CenterTitle(kTRUE);
    s6->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s6->GetXaxis()->SetLabelFont(42);
    s6->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s6->GetXaxis()->SetTitleSize(0.04);
    s6->GetXaxis()->SetLabelSize(0.03);
    s6->GetXaxis()->CenterTitle(kTRUE);
    c6->SaveAs(Form("pp500_EM_Mass_wide_all_pT.pdf"));
    c6->Close();

// Signal wide window p_T 1-2 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s7 = (TH1F*) d1 -> Get("D0 S12w;1");
    TH1F* b7 = (TH1F*) d1 -> Get("D0 B12w;1");
    TH1F* diff7 = (TH1F*) d1 -> Get("D0 rozdil12w;1");
    TH1F* e7 = (TH1F*) d3 -> Get("D0 S12w;1");

    TCanvas *c7 = new TCanvas("c7", "Invariant mass", 1400, 1000);

    TH1F *r7 = new TH1F();
    *r7 = 1*(*diff7);

    s7->Scale(1/s6->GetEntries());
    b7->Scale(1/b6->GetEntries());
    e7->Scale(1/e6->GetEntries());
    
    s7->Draw("PE");
    b7->Draw("PE same");
    r7->Draw("PE same");
    e7->Draw("HIST same C");

    e7->SetLineColor(kGreen);
    e7->SetLineWidth(1.5);

    c7->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c7->Update();

    s7->SetMarkerStyle(22);
    b7->SetMarkerStyle(23);
    r7->SetMarkerStyle(33);

    s7->SetMarkerSize(0.7);
    b7->SetMarkerSize(0.7);
    r7->SetMarkerSize(0.9);

    s7->SetMarkerColor(kBlue);
    b7->SetMarkerColor(kRed);
    r7->SetMarkerColor(kBlack);

    s7->SetTitle(" ; ");
    s7->SetStats(0);

    TLegend *l7 = new TLegend(0.55,0.65, 0.73, 0.89,"","brNDC");
    l7->AddEntry(s7, "Unlike-sign K^{#pm}#pi^{#mp} pairs", "pl");
    l7->AddEntry(b7, "Like-sign K^{#pm}#pi^{#pm} pairs (background)", "pl");
//    l7->AddEntry(r7, "Difference US-LS", "pl");
    l7->AddEntry(e7, "Event-Mixing (background)", "pl");
    l7->AddEntry((TObject*)0, "", "");
    l7->AddEntry((TObject*)0, "THIS THESIS", "");
    l7->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l7->AddEntry((TObject*)0, "1 < #it{p}_{T}(K#pi) < 2 GeV/#it{c}^{2}", "");

    l7->SetFillStyle(0);
    l7->SetLineColor(0);
    l7->SetTextSize(0.03);
    l7->Draw("same");

    s7->GetYaxis()->SetTitle("Counts [-]");
    s7->GetYaxis()->SetTitleOffset(1.);
    s7->GetYaxis()->SetLabelSize(0.03);
    s7->GetYaxis()->SetTitleSize(0.04);
    s7->GetYaxis()->SetLabelFont(42);
    s7->GetYaxis()->SetTitleFont(42);
//    s7->GetYaxis()->SetRangeUser(0,250000);
    s7->GetYaxis()->CenterTitle(kTRUE);
    s7->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s7->GetXaxis()->SetLabelFont(42);
    s7->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s7->GetXaxis()->SetTitleSize(0.04);
    s7->GetXaxis()->SetLabelSize(0.03);
    s7->GetXaxis()->CenterTitle(kTRUE);
    c7->SaveAs(Form("pp500_EM_Mass_wide_12_pT.pdf"));
    c7->Close();

// Signal wide window p_T 2-3 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s8 = (TH1F*) d1 -> Get("D0 S23w;1");
    TH1F* b8 = (TH1F*) d1 -> Get("D0 B23w;1");
    TH1F* diff8 = (TH1F*) d1 -> Get("D0 rozdil23w;1");
    TH1F* e8 = (TH1F*) d3 -> Get("D0 S23w;1");

    TCanvas *c8 = new TCanvas("c8", "Invariant mass", 1400, 1000);

    TH1F *r8 = new TH1F();
    *r8 = 1*(*diff8);

    s8->Scale(1/s6->GetEntries());
    b8->Scale(1/b6->GetEntries());
    e8->Scale(1/e6->GetEntries());
    
    s8->Draw("PE");
    b8->Draw("PE same");
    r8->Draw("PE same");    
    e8->Draw("HIST same C");

    e8->SetLineColor(kGreen);
    e8->SetLineWidth(1.5);
    
    c8->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c8->Update();

    s8->SetMarkerStyle(22);
    b8->SetMarkerStyle(23);
    r8->SetMarkerStyle(33);

    s8->SetMarkerSize(0.7);
    b8->SetMarkerSize(0.7);
    r8->SetMarkerSize(0.9);

    s8->SetMarkerColor(kBlue);
    b8->SetMarkerColor(kRed);
    r8->SetMarkerColor(kBlack);

    s8->SetTitle(" ; ");
    s8->SetStats(0);

    TLegend *l8 = new TLegend(0.55,0.65, 0.73, 0.89,"","brNDC");
    l8->AddEntry(s8, "Unlike-sign K^{#pm}#pi^{#mp} pairs", "pl");
    l8->AddEntry(b8, "Like-sign K^{#pm}#pi^{#pm} pairs (background)", "pl");
    l8->AddEntry(r8, "Difference US-LS", "pl");
    l8->AddEntry((TObject*)0, "", "");
    l8->AddEntry((TObject*)0, "THIS THESIS", "");
    l8->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l8->AddEntry((TObject*)0, "2 < #it{p}_{T}(K#pi) < 3 GeV/#it{c}^{2}", "");

    l8->SetFillStyle(0);
    l8->SetLineColor(0);
    l8->SetTextSize(0.03);
    l8->Draw("same");

    s8->GetYaxis()->SetTitle("Counts [-]");
    s8->GetYaxis()->SetTitleOffset(1.);
    s8->GetYaxis()->SetLabelSize(0.03);
    s8->GetYaxis()->SetTitleSize(0.04);
    s8->GetYaxis()->SetLabelFont(42);
    s8->GetYaxis()->SetTitleFont(42);
//    s8->GetYaxis()->SetRangeUser(0,100000);
    s8->GetYaxis()->CenterTitle(kTRUE);
    s8->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s8->GetXaxis()->SetLabelFont(42);
    s8->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s8->GetXaxis()->SetTitleSize(0.04);
    s8->GetXaxis()->SetLabelSize(0.03);
    s8->GetXaxis()->CenterTitle(kTRUE);
    c8->SaveAs(Form("pp500_EM_Mass_wide_23_pT.pdf"));
    c8->Close();

// Signal wide window p_T 3-4 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s9 = (TH1F*) d1 -> Get("D0 S34w;1");
    TH1F* b9 = (TH1F*) d1 -> Get("D0 B34w;1");
    TH1F* diff9 = (TH1F*) d1 -> Get("D0 rozdil34w;1");
    TH1F* e9 = (TH1F*) d3 -> Get("D0 S34w;1");

    TCanvas *c9 = new TCanvas("c9", "Invariant mass", 1400, 1000);

    TH1F *r9 = new TH1F();
    *r9 = 1*(*diff9);

    s9->Scale(1/s6->GetEntries());
    b9->Scale(1/b6->GetEntries());
    e9->Scale(1/e6->GetEntries());

    s9->Draw("PE");
    b9->Draw("PE same");
    r9->Draw("PE same");
    e9->Draw("HIST same C");

    e9->SetLineColor(kGreen);
    e9->SetLineWidth(1.5);
    
    c9->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c9->Update();

    s9->SetMarkerStyle(22);
    b9->SetMarkerStyle(23);
    r9->SetMarkerStyle(33);

    s9->SetMarkerSize(0.7);
    b9->SetMarkerSize(0.7);
    r9->SetMarkerSize(0.9);

    s9->SetMarkerColor(kBlue);
    b9->SetMarkerColor(kRed);
    r9->SetMarkerColor(kBlack);

    s9->SetTitle(" ; ");
    s9->SetStats(0);

    TLegend *l9 = new TLegend(0.55,0.65, 0.73, 0.89,"","brNDC");
    l9->AddEntry(s9, "Unlike-sign K^{#pm}#pi^{#mp} pairs", "pl");
    l9->AddEntry(b9, "Like-sign K^{#pm}#pi^{#pm} pairs (background)", "pl");
    l9->AddEntry(r9, "Difference US-LS", "pl");
    l9->AddEntry((TObject*)0, "", "");
    l9->AddEntry((TObject*)0, "THIS THESIS", "");
    l9->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l9->AddEntry((TObject*)0, "3 < #it{p}_{T}(K#pi) < 4 GeV/#it{c}^{2}", "");

    l9->SetFillStyle(0);
    l9->SetLineColor(0);
    l9->SetTextSize(0.03);
    l9->Draw("same");

    s9->GetYaxis()->SetTitle("Counts [-]");
    s9->GetYaxis()->SetTitleOffset(1.);
    s9->GetYaxis()->SetLabelSize(0.03);
    s9->GetYaxis()->SetTitleSize(0.04);
    s9->GetYaxis()->SetLabelFont(42);
    s9->GetYaxis()->SetTitleFont(42);
//    s9->GetYaxis()->SetRangeUser(0,50000);
    s9->GetYaxis()->CenterTitle(kTRUE);
    s9->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s9->GetXaxis()->SetLabelFont(42);
    s9->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s9->GetXaxis()->SetTitleSize(0.04);
    s9->GetXaxis()->SetLabelSize(0.03);
    s9->GetXaxis()->CenterTitle(kTRUE);
    c9->SaveAs(Form("pp500_EM_Mass_wide_34_pT.pdf"));
    c9->Close();

// Signal wide window p_T >4 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s10 = (TH1F*) d1 -> Get("D0 S45w;1");
    TH1F* b10 = (TH1F*) d1 -> Get("D0 B45w;1");
    TH1F* diff10 = (TH1F*) d1 -> Get("D0 rozdil45w;1");
    TH1F* e10 = (TH1F*) d3 -> Get("D0 S45w;1");

    TCanvas *c10 = new TCanvas("c10", "Invariant mass", 1400, 1000);

    TH1F *r10 = new TH1F();
    *r10 = 1*(*diff10);

    s10->Scale(1/s6->GetEntries());
    b10->Scale(1/b6->GetEntries());
    e10->Scale(1/e6->GetEntries());

    s10->Draw("PE");
    b10->Draw("PE same");
    r10->Draw("PE same");
    e10->Draw("HIST same C");

    e10->SetLineColor(kGreen);
    e10->SetLineWidth(1.5);
    
    c10->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c10->Update();

    s10->SetMarkerStyle(22);
    b10->SetMarkerStyle(23);
    r10->SetMarkerStyle(33);

    s10->SetMarkerSize(0.7);
    b10->SetMarkerSize(0.7);
    r10->SetMarkerSize(0.9);

    s10->SetMarkerColor(kBlue);
    b10->SetMarkerColor(kRed);
    r10->SetMarkerColor(kBlack);

    s10->SetTitle(" ; ");
    s10->SetStats(0);

    TLegend *l10 = new TLegend(0.55,0.65, 0.73, 0.89,"","brNDC");
    l10->AddEntry(s10, "Unlike-sign K^{#pm}#pi^{#mp} pairs", "pl");
    l10->AddEntry(b10, "Like-sign K^{#pm}#pi^{#pm} pairs (background)", "pl");
    l10->AddEntry(r10, "Difference US-LS", "pl");
    l10->AddEntry((TObject*)0, "", "");
    l10->AddEntry((TObject*)0, "THIS THESIS", "");
    l10->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l10->AddEntry((TObject*)0, "#it{p}_{T}(K#pi) > 4 GeV/#it{c}^{2}", "");

    l10->SetFillStyle(0);
    l10->SetLineColor(0);
    l10->SetTextSize(0.03);
    l10->Draw("same");

    s10->GetYaxis()->SetTitle("Counts [-]");
    s10->GetYaxis()->SetTitleOffset(1.);
    s10->GetYaxis()->SetLabelSize(0.03);
    s10->GetYaxis()->SetTitleSize(0.04);
    s10->GetYaxis()->SetLabelFont(42);
    s10->GetYaxis()->SetTitleFont(42);
//    s10->GetYaxis()->SetRangeUser(0,50000);
    s10->GetYaxis()->CenterTitle(kTRUE);
    s10->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s10->GetXaxis()->SetLabelFont(42);
    s10->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s10->GetXaxis()->SetTitleSize(0.04);
    s10->GetXaxis()->SetLabelSize(0.03);
    s10->GetXaxis()->CenterTitle(kTRUE);
    c10->SaveAs(Form("pp500_EM_Mass_wide_45_pT.pdf"));
    c10->Close();

















    cout<<"Michale jsi dobrej"<<endl;
}
