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

void Graphs_final() {

    TFile* d1 = new TFile("pp500_D0.root");
    TFile* d2 = new TFile("pp500_Dstar.root");

// Signal narrow window all p_T -----------------------------------------------------------------------------------------------------------------------------
    TH1F* s1 = (TH1F*) d1 -> Get("D0 signal;1");
    TH1F* b1 = (TH1F*) d1 -> Get("D0 background;1");
    TH1F* diff1 = (TH1F*) d1 -> Get("D0 rozdil;1");

    TCanvas *c1 = new TCanvas("c1", "Invariant mass", 1400, 1000);
    //  gPad->SetMargin(2.5,0.1,0.03,0.03);
    // pionpt->SetLogy();

    // s1->Scale(1/s1->GetEntries());
    // b1->Scale(1/b1->GetEntries());

    TH1F *r1 = new TH1F();
    *r1 = 10*(*diff1);

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
    l1->AddEntry(b1, "Like-sign K^{#pm}#pi^{#pm} pairs (background)", "pl");
    l1->AddEntry(r1, "Difference US-LS scaled by 10", "pl");
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
    s1->GetYaxis()->SetRangeUser(0,2100000);
    s1->GetYaxis()->CenterTitle(kTRUE);
    s1->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s1->GetXaxis()->SetLabelFont(42);
    s1->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s1->GetXaxis()->SetTitleSize(0.04);
    s1->GetXaxis()->SetLabelSize(0.03);
    s1->GetXaxis()->CenterTitle(kTRUE);
    c1->SaveAs(Form("pp500_PionMatched_Mass_narrow_all_pT.pdf"));
    c1->Close();

// Signal narrow window p_T 1-2 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s2 = (TH1F*) d1 -> Get("D0 S12;1");
    TH1F* b2 = (TH1F*) d1 -> Get("D0 B12;1");
    TH1F* diff2 = (TH1F*) d1 -> Get("D0 rozdil12;1");

    TCanvas *c2 = new TCanvas("c2", "Invariant mass", 1400, 1000);

    TH1F *r2 = new TH1F();
    *r2 = 10*(*diff2);

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
    l2->AddEntry(b2, "Like-sign K^{#pm}#pi^{#pm} pairs (background)", "pl");
    l2->AddEntry(r2, "Difference US-LS scaled by 10", "pl");
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
    s2->GetYaxis()->SetRangeUser(0,800000);
    s2->GetYaxis()->CenterTitle(kTRUE);
    s2->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s2->GetXaxis()->SetLabelFont(42);
    s2->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s2->GetXaxis()->SetTitleSize(0.04);
    s2->GetXaxis()->SetLabelSize(0.03);
    s2->GetXaxis()->CenterTitle(kTRUE);
    c2->SaveAs(Form("pp500_PionMatched_Mass_narrow_12_pT.pdf"));
    c2->Close();

// Signal narrow window p_T 2-3 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s3 = (TH1F*) d1 -> Get("D0 S23;1");
    TH1F* b3 = (TH1F*) d1 -> Get("D0 B23;1");
    TH1F* diff3 = (TH1F*) d1 -> Get("D0 rozdil23;1");

    TCanvas *c3 = new TCanvas("c3", "Invariant mass", 1400, 1000);

    TH1F *r3 = new TH1F();
    *r3 = 10*(*diff3);

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
    l3->AddEntry(b3, "Like-sign K^{#pm}#pi^{#pm} pairs (background)", "pl");
    l3->AddEntry(r3, "Difference US-LS scaled by 10", "pl");
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
    s3->GetYaxis()->SetRangeUser(0,250000);
    s3->GetYaxis()->CenterTitle(kTRUE);
    s3->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s3->GetXaxis()->SetLabelFont(42);
    s3->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s3->GetXaxis()->SetTitleSize(0.04);
    s3->GetXaxis()->SetLabelSize(0.03);
    s3->GetXaxis()->CenterTitle(kTRUE);
    c3->SaveAs(Form("pp500_PionMatched_Mass_narrow_23_pT.pdf"));
    c3->Close();

// Signal narrow window p_T 3-4 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s4 = (TH1F*) d1 -> Get("D0 S34;1");
    TH1F* b4 = (TH1F*) d1 -> Get("D0 B34;1");
    TH1F* diff4 = (TH1F*) d1 -> Get("D0 rozdil34;1");

    TCanvas *c4 = new TCanvas("c4", "Invariant mass", 1400, 1000);

    TH1F *r4 = new TH1F();
    *r4 = 5*(*diff4);

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
    l4->AddEntry(b4, "Like-sign K^{#pm}#pi^{#pm} pairs (background)", "pl");
    l4->AddEntry(r4, "Difference US-LS scaled by 5", "pl");
    l4->AddEntry((TObject*)0, "", "");
    l4->AddEntry((TObject*)0, "THIS THESIS", "");
    l4->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l4->AddEntry((TObject*)0, "3 < #it{p}_{T}(K#pi) < 4 GeV/#it{c}^{2}", "");

    l4->SetFillStyle(0);
    l4->SetLineColor(0);
    l4->SetTextSize(0.03);
    l4->Draw("same");

    s4->GetYaxis()->SetTitle("Counts [-]");
    s4->GetYaxis()->SetTitleOffset(1.);
    s4->GetYaxis()->SetLabelSize(0.03);
    s4->GetYaxis()->SetTitleSize(0.04);
    s4->GetYaxis()->SetLabelFont(42);
    s4->GetYaxis()->SetTitleFont(42);
    s4->GetYaxis()->SetRangeUser(0,50000);
    s4->GetYaxis()->CenterTitle(kTRUE);
    s4->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s4->GetXaxis()->SetLabelFont(42);
    s4->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s4->GetXaxis()->SetTitleSize(0.04);
    s4->GetXaxis()->SetLabelSize(0.03);
    s4->GetXaxis()->CenterTitle(kTRUE);
    c4->SaveAs(Form("pp500_PionMatched_Mass_narrow_34_pT.pdf"));
    c4->Close();

// Signal narrow window p_T >4 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s5 = (TH1F*) d1 -> Get("D0 S45;1");
    TH1F* b5 = (TH1F*) d1 -> Get("D0 B45;1");
    TH1F* diff5 = (TH1F*) d1 -> Get("D0 rozdil45;1");

    TCanvas *c5 = new TCanvas("c5", "Invariant mass", 1400, 1000);

    TH1F *r5 = new TH1F();
    *r5 = 5*(*diff5);

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
    l5->AddEntry(b5, "Like-sign K^{#pm}#pi^{#pm} pairs (background)", "pl");
    l5->AddEntry(r5, "Difference US-LS scaled by 5", "pl");
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
    s5->GetYaxis()->SetRangeUser(0,15000);
    s5->GetYaxis()->CenterTitle(kTRUE);
    s5->GetXaxis()->SetTitle("#it{M}_{inv} [GeV/#it{c}^{2}]");
    s5->GetXaxis()->SetLabelFont(42);
    s5->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s5->GetXaxis()->SetTitleSize(0.04);
    s5->GetXaxis()->SetLabelSize(0.03);
    s5->GetXaxis()->CenterTitle(kTRUE);
    c5->SaveAs(Form("pp500_PionMatched_Mass_narrow_45_pT.pdf"));
    c5->Close();

// Signal wide window all p_T -----------------------------------------------------------------------------------------------------------------------------
    TH1F* s6 = (TH1F*) d1 -> Get("D0 signalw;1");
    TH1F* b6 = (TH1F*) d1 -> Get("D0 backgroundw;1");
    TH1F* diff6 = (TH1F*) d1 -> Get("D0 rozdilw;1");

    TCanvas *c6 = new TCanvas("c6", "Invariant mass", 1400, 1000);

    TH1F *r6 = new TH1F();
    *r6 = 1*(*diff6);

    s6->Draw("PE");
    b6->Draw("PE same");
    r6->Draw("PE same");
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
    l6->AddEntry(r6, "Difference US-LS", "pl");
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
    c6->SaveAs(Form("pp500_PionMatched_Mass_wide_all_pT.pdf"));
    c6->Close();

// Signal wide window p_T 1-2 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s7 = (TH1F*) d1 -> Get("D0 S12w;1");
    TH1F* b7 = (TH1F*) d1 -> Get("D0 B12w;1");
    TH1F* diff7 = (TH1F*) d1 -> Get("D0 rozdil12w;1");

    TCanvas *c7 = new TCanvas("c7", "Invariant mass", 1400, 1000);

    TH1F *r7 = new TH1F();
    *r7 = 1*(*diff7);

    s7->Draw("PE");
    b7->Draw("PE same");
    r7->Draw("PE same");
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
    l7->AddEntry(r7, "Difference US-LS", "pl");
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
    c7->SaveAs(Form("pp500_PionMatched_Mass_wide_12_pT.pdf"));
    c7->Close();

// Signal wide window p_T 2-3 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s8 = (TH1F*) d1 -> Get("D0 S23w;1");
    TH1F* b8 = (TH1F*) d1 -> Get("D0 B23w;1");
    TH1F* diff8 = (TH1F*) d1 -> Get("D0 rozdil23w;1");

    TCanvas *c8 = new TCanvas("c8", "Invariant mass", 1400, 1000);

    TH1F *r8 = new TH1F();
    *r8 = 1*(*diff8);

    s8->Draw("PE");
    b8->Draw("PE same");
    r8->Draw("PE same");
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
    c8->SaveAs(Form("pp500_PionMatched_Mass_wide_23_pT.pdf"));
    c8->Close();

// Signal wide window p_T 3-4 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s9 = (TH1F*) d1 -> Get("D0 S34w;1");
    TH1F* b9 = (TH1F*) d1 -> Get("D0 B34w;1");
    TH1F* diff9 = (TH1F*) d1 -> Get("D0 rozdil34w;1");

    TCanvas *c9 = new TCanvas("c9", "Invariant mass", 1400, 1000);

    TH1F *r9 = new TH1F();
    *r9 = 1*(*diff9);

    s9->Draw("PE");
    b9->Draw("PE same");
    r9->Draw("PE same");
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
    c9->SaveAs(Form("pp500_PionMatched_Mass_wide_34_pT.pdf"));
    c9->Close();

// Signal wide window p_T >4 GeV-----------------------------------------------------------------------------------------------------------------------------
    TH1F* s10 = (TH1F*) d1 -> Get("D0 S45w;1");
    TH1F* b10 = (TH1F*) d1 -> Get("D0 B45w;1");
    TH1F* diff10 = (TH1F*) d1 -> Get("D0 rozdil45w;1");

    TCanvas *c10 = new TCanvas("c10", "Invariant mass", 1400, 1000);

    TH1F *r10 = new TH1F();
    *r10 = 1*(*diff10);

    s10->Draw("PE");
    b10->Draw("PE same");
    r10->Draw("PE same");
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
    c10->SaveAs(Form("pp500_PionMatched_Mass_wide_45_pT.pdf"));
    c10->Close();

// Dstar all p_T -----------------------------------------------------------------------------------------------------------------------------
    TH1F* s11 = (TH1F*) d2 -> Get("Dstar-D0;1");
    TH1F* back11 = (TH1F*) d2 -> Get("bDstar-D0;1");
    TH1F* diff11 = (TH1F*) d2 -> Get("Dstar-D0 rozdil;1");

    TCanvas *c11 = new TCanvas("c11", "Invariant mass", 1400, 1000);

    TH1F *r11 = new TH1F();
    *r11 = 1*(*diff11);

    TH1F *b11 = new TH1F();
    *b11 = 0.33333333333333333333333333333333*(*back11);

    s11->Draw("PE");
    b11->Draw("PE same");
    r11->Draw("PE same");
    c11->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c11->Update();
    

    s11->SetMarkerStyle(22);
    b11->SetMarkerStyle(23);
    r11->SetMarkerStyle(33);

    s11->SetMarkerSize(1.2);
    b11->SetMarkerSize(1.2);
    r11->SetMarkerSize(1.8);

    s11->SetMarkerColor(kBlue);
    b11->SetMarkerColor(kRed);
    r11->SetMarkerColor(kBlack);

    s11->SetTitle(" ; ");
    s11->SetStats(0);

    TLegend *l11 = new TLegend(0.15,0.69, 0.35, 0.89,"","brNDC");
    l11->AddEntry(s11, "Correct-sign K^{#pm}#pi^{#mp}#pi^{#mp} triplets", "pl");
    l11->AddEntry(b11, "Wrong-sign K^{#pm}#pi^{#pm}#pi^{#pm} triplets (background)", "pl");
    l11->AddEntry(r11, "Difference US-LS", "pl");
    l11->AddEntry((TObject*)0, "", "");
    l11->AddEntry((TObject*)0, "THIS THESIS", "");
    l11->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");

    l11->SetFillStyle(0);
    l11->SetLineColor(0);
    l11->SetTextSize(0.03);
    l11->Draw("same");

    s11->GetYaxis()->SetTitle("Counts [-]");
    s11->GetYaxis()->SetTitleOffset(1.);
    s11->GetYaxis()->SetLabelSize(0.03);
    s11->GetYaxis()->SetTitleSize(0.04);
    s11->GetYaxis()->SetLabelFont(42);
    s11->GetYaxis()->SetTitleFont(42);
    s11->GetYaxis()->SetRangeUser(-100,1800);
    s11->GetYaxis()->CenterTitle(kTRUE);
    s11->GetXaxis()->SetTitle("#it{M}_{K#pi#pi}-#it{M}_{K#pi} [GeV/#it{c}^{2}]");
    s11->GetXaxis()->SetLabelFont(42);
    s11->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s11->GetXaxis()->SetTitleSize(0.04);
    s11->GetXaxis()->SetLabelSize(0.03);
    s11->GetXaxis()->CenterTitle(kTRUE);
    c11->SaveAs(Form("pp500_PionMatched_Mass_Dstar-D0_all_pT.pdf"));
    c11->Close();

// Dstar all p_T 2-3 GeV ----------------------------------------------------------------------------------------------------------------------------
    TH1F* s12 = (TH1F*) d2 -> Get("Dstar-D0 S23;1");
    TH1F* back12 = (TH1F*) d2 -> Get("bDstar-D0 B23;1");
    TH1F* diff12 = (TH1F*) d2 -> Get("Dstar-D0 rozdil23;1");

    TCanvas *c12 = new TCanvas("c12", "Invariant mass", 1400, 1000);

    TH1F *r12 = new TH1F();
    *r12 = 1*(*diff12);

    TH1F *b12 = new TH1F();
    *b12 = 0.33333333333333333333333333333333*(*back12);

    s12->Draw("PE");
    b12->Draw("PE same");
    r12->Draw("PE same");
    c12->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c12->Update();

    s12->SetMarkerStyle(22);
    b12->SetMarkerStyle(23);
    r12->SetMarkerStyle(33);

    s12->SetMarkerSize(1.2);
    b12->SetMarkerSize(1.2);
    r12->SetMarkerSize(1.8);

    s12->SetMarkerColor(kBlue);
    b12->SetMarkerColor(kRed);
    r12->SetMarkerColor(kBlack);

    s12->SetTitle(" ; ");
    s12->SetStats(0);

    TLegend *l12 = new TLegend(0.15,0.69, 0.35, 0.89,"","brNDC");
    l12->AddEntry(s12, "Correct-sign K^{#pm}#pi^{#mp}#pi^{#mp} triplets", "pl");
    l12->AddEntry(b12, "Wrong-sign K^{#pm}#pi^{#pm}#pi^{#pm} triplets (background)", "pl");
    l12->AddEntry(r12, "Difference US-LS", "pl");
    l12->AddEntry((TObject*)0, "", "");
    l12->AddEntry((TObject*)0, "THIS THESIS", "");
    l12->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l12->AddEntry((TObject*)0, "2 < #it{p}_{T}(K#pi#pi) < 3 GeV/#it{c}^{2}", "");

    l12->SetFillStyle(0);
    l12->SetLineColor(0);
    l12->SetTextSize(0.03);
    l12->Draw("same");

    s12->GetYaxis()->SetTitle("Counts [-]");
    s12->GetYaxis()->SetTitleOffset(1.);
    s12->GetYaxis()->SetLabelSize(0.03);
    s12->GetYaxis()->SetTitleSize(0.04);
    s12->GetYaxis()->SetLabelFont(42);
    s12->GetYaxis()->SetTitleFont(42);
    s12->GetYaxis()->SetRangeUser(-100,1400);
    s12->GetYaxis()->CenterTitle(kTRUE);
    s12->GetXaxis()->SetTitle("#it{M}_{K#pi#pi}-#it{M}_{K#pi} [GeV/#it{c}^{2}]");
    s12->GetXaxis()->SetLabelFont(42);
    s12->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s12->GetXaxis()->SetTitleSize(0.04);
    s12->GetXaxis()->SetLabelSize(0.03);
    s12->GetXaxis()->CenterTitle(kTRUE);
    c12->SaveAs(Form("pp500_PionMatched_Mass_Dstar-D0_23_pT.pdf"));
    c12->Close();

// Dstar all p_T 3-4 GeV ----------------------------------------------------------------------------------------------------------------------------
    TH1F* s13 = (TH1F*) d2 -> Get("Dstar-D0 S34;1");
    TH1F* back13 = (TH1F*) d2 -> Get("bDstar-D0 B34;1");
    TH1F* diff13 = (TH1F*) d2 -> Get("Dstar-D0 rozdil34;1");

    TCanvas *c13 = new TCanvas("c13", "Invariant mass", 1400, 1000);

    TH1F *r13 = new TH1F();
    *r13 = 1*(*diff13);

    TH1F *b13 = new TH1F();
    *b13 = 0.33333333333333333333333333333333*(*back13);

    s13->Draw("PE");
    b13->Draw("PE same");
    r13->Draw("PE same");
    c13->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c13->Update();

    s13->SetMarkerStyle(22);
    b13->SetMarkerStyle(23);
    r13->SetMarkerStyle(33);

    s13->SetMarkerSize(1.2);
    b13->SetMarkerSize(1.2);
    r13->SetMarkerSize(1.8);

    s13->SetMarkerColor(kBlue);
    b13->SetMarkerColor(kRed);
    r13->SetMarkerColor(kBlack);

    s13->SetTitle(" ; ");
    s13->SetStats(0);

    TLegend *l13 = new TLegend(0.15,0.69, 0.35, 0.89,"","brNDC");
    l13->AddEntry(s13, "Correct-sign K^{#pm}#pi^{#mp}#pi^{#mp} triplets", "pl");
    l13->AddEntry(b13, "Wrong-sign K^{#pm}#pi^{#pm}#pi^{#pm} triplets (background)", "pl");
    l13->AddEntry(r13, "Difference US-LS", "pl");
    l13->AddEntry((TObject*)0, "", "");
    l13->AddEntry((TObject*)0, "THIS THESIS", "");
    l13->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l13->AddEntry((TObject*)0, "3 < #it{p}_{T}(K#pi#pi) < 4 GeV/#it{c}^{2}", "");

    l13->SetFillStyle(0);
    l13->SetLineColor(0);
    l13->SetTextSize(0.03);
    l13->Draw("same");

    s13->GetYaxis()->SetTitle("Counts [-]");
    s13->GetYaxis()->SetTitleOffset(1.);
    s13->GetYaxis()->SetLabelSize(0.03);
    s13->GetYaxis()->SetTitleSize(0.04);
    s13->GetYaxis()->SetLabelFont(42);
    s13->GetYaxis()->SetTitleFont(42);
    s13->GetYaxis()->SetRangeUser(-50,400);
    s13->GetYaxis()->CenterTitle(kTRUE);
    s13->GetXaxis()->SetTitle("#it{M}_{K#pi#pi}-#it{M}_{K#pi} [GeV/#it{c}^{2}]");
    s13->GetXaxis()->SetLabelFont(42);
    s13->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s13->GetXaxis()->SetTitleSize(0.04);
    s13->GetXaxis()->SetLabelSize(0.03);
    s13->GetXaxis()->CenterTitle(kTRUE);
    c13->SaveAs(Form("pp500_PionMatched_Mass_Dstar-D0_34_pT.pdf"));
    c13->Close();

// Dstar all p_T >4 GeV ----------------------------------------------------------------------------------------------------------------------------
    TH1F* s14 = (TH1F*) d2 -> Get("Dstar-D0 S45;1");
    TH1F* back14 = (TH1F*) d2 -> Get("bDstar-D0 B45;1");
    TH1F* diff14 = (TH1F*) d2 -> Get("Dstar-D0 rozdil45;1");

    TCanvas *c14 = new TCanvas("c14", "Invariant mass", 1400, 1000);

    TH1F *r14 = new TH1F();
    *r14 = 1*(*diff14);

    TH1F *b14 = new TH1F();
    *b14 = 0.33333333333333333333333333333333*(*back14);

    s14->Draw("PE");
    b14->Draw("PE same");
    r14->Draw("PE same");
    c14->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c14->Update();

    s14->SetMarkerStyle(22);
    b14->SetMarkerStyle(23);
    r14->SetMarkerStyle(33);

    s14->SetMarkerSize(1.2);
    b14->SetMarkerSize(1.2);
    r14->SetMarkerSize(1.8);

    s14->SetMarkerColor(kBlue);
    b14->SetMarkerColor(kRed);
    r14->SetMarkerColor(kBlack);

    s14->SetTitle(" ; ");
    s14->SetStats(0);

    TLegend *l14 = new TLegend(0.15,0.69, 0.35, 0.89,"","brNDC");
    l14->AddEntry(s14, "Correct-sign K^{#pm}#pi^{#mp}#pi^{#mp} triplets", "pl");
    l14->AddEntry(b14, "Wrong-sign K^{#pm}#pi^{#pm}#pi^{#pm} triplets (background)", "pl");
    l14->AddEntry(r14, "Difference US-LS", "pl");
    l14->AddEntry((TObject*)0, "", "");
    l14->AddEntry((TObject*)0, "THIS THESIS", "");
    l14->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");
    l14->AddEntry((TObject*)0, "#it{p}_{T}(K#pi#pi) > 4 GeV/#it{c}^{2}", "");

    l14->SetFillStyle(0);
    l14->SetLineColor(0);
    l14->SetTextSize(0.03);
    l14->Draw("same");

    s14->GetYaxis()->SetTitle("Counts [-]");
    s14->GetYaxis()->SetTitleOffset(1.);
    s14->GetYaxis()->SetLabelSize(0.03);
    s14->GetYaxis()->SetTitleSize(0.04);
    s14->GetYaxis()->SetLabelFont(42);
    s14->GetYaxis()->SetTitleFont(42);
    s14->GetYaxis()->SetRangeUser(-30,130);
    s14->GetYaxis()->CenterTitle(kTRUE);
    s14->GetXaxis()->SetTitle("#it{M}_{K#pi#pi}-#it{M}_{K#pi} [GeV/#it{c}^{2}]");
    s14->GetXaxis()->SetLabelFont(42);
    s14->GetXaxis()->SetTitleFont(42);
    //  s1->GetXaxis()->SetRangeUser(0,6);
    s14->GetXaxis()->SetTitleSize(0.04);
    s14->GetXaxis()->SetLabelSize(0.03);
    s14->GetXaxis()->CenterTitle(kTRUE);
    c14->SaveAs(Form("pp500_PionMatched_Mass_Dstar-D0_45_pT.pdf"));
    c14->Close();























    cout<<"Michale jsi dobrej"<<endl;
}
