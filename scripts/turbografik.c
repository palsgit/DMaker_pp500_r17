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

void turbografik() {

    TFile* d1 = new TFile("D0.root");

// Signal narrow window -----------------------------------------------------------------------------------------------------------------------------
    TH1F* h1 = (TH1F*) d1 -> Get("signal;1");
    TH1F* b1 = (TH1F*) d1 -> Get("background;1");
    TCanvas *c1 = new TCanvas("c1", "Invariant mass", 1400, 1100);
    //  gPad->SetMargin(2.5,0.1,0.03,0.03);
    // pionpt->SetLogy();

    // h1->Scale(1/h1->GetEntries());
    // b1->Scale(1/b1->GetEntries());
    h1->Draw("PE");
    b1->Draw("PE same");
    c1->Update();
    h1->SetMarkerStyle(8);
    b1->SetMarkerStyle(8);
    h1->SetMarkerColor(kBlue);
    b1->SetMarkerColor(kRed);
    h1->SetLineWidth(0);
    b1->SetLineWidth(0);
    h1->SetTitle(" ; ");
    h1->SetStats(0);
    h1->SetLineColor(kBlue);
    b1->SetLineColor(kRed);
    /* TText *text = new TText(1.291488,746304.2,"THIS THESIS");
     text->SetTextAlign(23);
     text->SetTextSize(.04);
     text->Draw();*/
    TLegend *l1 = new TLegend(0.55,0.80, 0.8, 0.907,"","brNDC");
    l1->AddEntry(h1, "Unlike-sign pairs", "pl");
    l1->AddEntry(b1, "Like-sign pairs - background", "pl");
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
    //  h1->GetYaxis()->SetRangeUser(0,30000);
    h1->GetYaxis()->CenterTitle(kTRUE);
    h1->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h1->GetXaxis()->SetLabelFont(42);
    h1->GetXaxis()->SetTitleFont(42);
    //  h1->GetXaxis()->SetRangeUser(0,6);
    h1->GetXaxis()->SetTitleSize(0.04);
    h1->GetXaxis()->SetLabelSize(0.03);
    h1->GetXaxis()->CenterTitle(kTRUE);
    c1->SaveAs(Form("Mass_narrow.pdf"));
    c1->Close();


// Signal full window -----------------------------------------------------------------------------------------------------------------------------
    TH1F* h2 = (TH1F*) d1 -> Get("signalw;1");
    TH1F* b2 = (TH1F*) d1 -> Get("backgroundw;1");
    TCanvas *c2 = new TCanvas("c2", "Invariant Mass", 1400, 1100);
    //  gPad->SetMargin(2.5,0.1,0.03,0.03);
    // h2->Scale(1/h2->GetEntries());
    // b2->Scale(1/b2->GetEntries());
    h2->Draw("PE");
    b2->Draw("PE same");
    c2->Update();
    h2->SetMarkerStyle(8);
    b2->SetMarkerStyle(8);
    h2->SetMarkerColor(kBlue);
    b2->SetMarkerColor(kRed);
    h2->SetLineWidth(0);
    b2->SetLineWidth(0);
    h2->SetTitle(" ; ");
    h2->SetStats(0);
    h2->SetLineColor(kBlue);
    b2->SetLineColor(kRed);
    /* TText *text = new TText(1.291488,746304.2,"THIS THESIS");
     text->SetTextAlign(23);
     text->SetTextSize(.04);
     text->Draw();*/
    TLegend *l2 = new TLegend(0.55,0.80, 0.8, 0.907,"","brNDC");
    l2->AddEntry(h2, "Unlike-sign pairs", "pl");
    l2->AddEntry(b2, "Like-sign pairs - background", "pl");
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
    //  h2->GetYaxis()->SetRangeUser(0,30000);
    h2->GetYaxis()->CenterTitle(kTRUE);
    h2->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h2->GetXaxis()->SetLabelFont(42);
    h2->GetXaxis()->SetTitleFont(42);
    //  h2->GetXaxis()->SetRangeUser(0,6);
    h2->GetXaxis()->SetTitleSize(0.04);
    h2->GetXaxis()->SetLabelSize(0.03);
    h2->GetXaxis()->CenterTitle(kTRUE);
    c2->SaveAs(Form("Mass_full.pdf"));
    c2->Close();


// Signal narrow window(1-2 GeV) -----------------------------------------------------------------------------------------------------------------------------
    TH1F* h3 = (TH1F*) d1 -> Get("S12;1");
    TH1F* b3 = (TH1F*) d1 -> Get("B12;1");
    TCanvas *c3 = new TCanvas("c3", "Invariant Mass (1-2 GeV)", 1400, 1100);
    //  gPad->SetMargin(2.5,0.1,0.03,0.03);
    // h3->Scale(1/h3->GetEntries());
    // b3->Scale(1/b3->GetEntries());
    h3->Draw("PE");
    b3->Draw("PE same");
    c3->Update();
    h3->SetMarkerStyle(8);
    b3->SetMarkerStyle(8);
    h3->SetMarkerColor(kBlue);
    b3->SetMarkerColor(kRed);
    h3->SetLineWidth(0);
    b3->SetLineWidth(0);
    h3->SetTitle(" ; ");
    h3->SetStats(0);
    h3->SetLineColor(kBlue);
    b3->SetLineColor(kRed);
    /* TText *text = new TText(1.291488,746304.2,"THIS THESIS");
     text->SetTextAlign(23);
     text->SetTextSize(.04);
     text->Draw();*/
    TLegend *l3 = new TLegend(0.55,0.80, 0.8, 0.907,"","brNDC");
    l3->AddEntry(h3, "Unlike-sign pairs", "pl");
    l3->AddEntry(b3, "Like-sign pairs - background", "pl");
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
    //  h3->GetYaxis()->SetRangeUser(0,30000);
    h3->GetYaxis()->CenterTitle(kTRUE);
    h3->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h3->GetXaxis()->SetLabelFont(42);
    h3->GetXaxis()->SetTitleFont(42);
    //  h3->GetXaxis()->SetRangeUser(0,6);
    h3->GetXaxis()->SetTitleSize(0.04);
    h3->GetXaxis()->SetLabelSize(0.03);
    h3->GetXaxis()->CenterTitle(kTRUE);
    c3->SaveAs(Form("Mass_narrow12.pdf"));
    c3->Close();


// Signal full window(1-2 GeV) -----------------------------------------------------------------------------------------------------------------------------
    TH1F* h4 = (TH1F*) d1 -> Get("S12w;1");
    TH1F* b4 = (TH1F*) d1 -> Get("B12w;1");
    TCanvas *c4 = new TCanvas("c4", "Invariant Mass (1-2 GeV)", 1400, 1100);
    //  gPad->SetMargin(2.5,0.1,0.03,0.03);
    // h4->Scale(1/h4->GetEntries());
    // b4->Scale(1/b4->GetEntries());
    h4->Draw("PE");
    b4->Draw("PE same");
    c4->Update();
    h4->SetMarkerStyle(8);
    b4->SetMarkerStyle(8);
    h4->SetMarkerColor(kBlue);
    b4->SetMarkerColor(kRed);
    h4->SetLineWidth(0);
    b4->SetLineWidth(0);
    h4->SetTitle(" ; ");
    h4->SetStats(0);
    h4->SetLineColor(kBlue);
    b4->SetLineColor(kRed);
    /* TText *text = new TText(1.291488,746304.2,"THIS THESIS");
     text->SetTextAlign(23);
     text->SetTextSize(.04);
     text->Draw();*/
    TLegend *l4 = new TLegend(0.55,0.80, 0.8, 0.907,"","brNDC");
    l4->AddEntry(h4, "Unlike-sign pairs", "pl");
    l4->AddEntry(b4, "Like-sign pairs - background", "pl");
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
    //  h4->GetYaxis()->SetRangeUser(0,30000);
    h4->GetYaxis()->CenterTitle(kTRUE);
    h4->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h4->GetXaxis()->SetLabelFont(42);
    h4->GetXaxis()->SetTitleFont(42);
    //  h4->GetXaxis()->SetRangeUser(0,6);
    h4->GetXaxis()->SetTitleSize(0.04);
    h4->GetXaxis()->SetLabelSize(0.03);
    h4->GetXaxis()->CenterTitle(kTRUE);
    c4->SaveAs(Form("Mass_full12.pdf"));
    c4->Close();


// Subtracted narrow window -----------------------------------------------------------------------------------------------------------------------------
    TH1F* h5 = (TH1F*) d1 -> Get("rozdil;1");
    TCanvas *c5 = new TCanvas("c5", "Background-Signal", 1400, 1100);
    h5->Draw("PE");
    c5->Update();
    h5->SetMarkerStyle(8);
    h5->SetMarkerColor(kBlue);
    h5->SetLineWidth(0);
    h5->SetTitle(" ; ");
    h5->SetStats(0);
    h5->SetLineColor(kBlue);
    TLegend *l5 = new TLegend(0.65,0.80, 0.8, 0.907,"","brNDC");
    l5->AddEntry(h5, "US-LS", "pl");
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
    //  h5->GetYaxis()->SetRangeUser(0,30000);
    h5->GetYaxis()->CenterTitle(kTRUE);
    h5->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h5->GetXaxis()->SetLabelFont(42);
    h5->GetXaxis()->SetTitleFont(42);
    //  h5->GetXaxis()->SetRangeUser(0,6);
    h5->GetXaxis()->SetTitleSize(0.04);
    h5->GetXaxis()->SetLabelSize(0.03);
    h5->GetXaxis()->CenterTitle(kTRUE);
    c5->SaveAs(Form("Mass_narrowD.pdf"));
    c5->Close();


// Subtracted full window -----------------------------------------------------------------------------------------------------------------------------
    TH1F* h6 = (TH1F*) d1 -> Get("rozdilw;1");
    TCanvas *c6 = new TCanvas("c6", "Background-Signal", 1400, 1100);
    h6->Draw("PE");
    c6->Update();
    h6->SetMarkerStyle(8);
    h6->SetMarkerColor(kBlue);
    h6->SetLineWidth(0);
    h6->SetTitle(" ; ");
    h6->SetStats(0);
    h6->SetLineColor(kBlue);
    TLegend *l6 = new TLegend(0.65,0.80, 0.8, 0.907,"","brNDC");
    l6->AddEntry(h6, "US-LS", "pl");
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
    //  h6->GetYaxis()->SetRangeUser(0,30000);
    h6->GetYaxis()->CenterTitle(kTRUE);
    h6->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h6->GetXaxis()->SetLabelFont(42);
    h6->GetXaxis()->SetTitleFont(42);
    //  h6->GetXaxis()->SetRangeUser(0,6);
    h6->GetXaxis()->SetTitleSize(0.04);
    h6->GetXaxis()->SetLabelSize(0.03);
    h6->GetXaxis()->CenterTitle(kTRUE);
    c6->SaveAs(Form("Mass_fullD.pdf"));
    c6->Close();


// Subtracted narrow window(1-2 GeV) -----------------------------------------------------------------------------------------------------------------------------
    TH1F* h7 = (TH1F*) d1 -> Get("rozdil12;1");
    TCanvas *c7 = new TCanvas("c7", "Background-Signal(1-2 GeV)", 1400, 1100);
    h7->Draw("PE");
    c7->Update();
    h7->SetMarkerStyle(8);
    h7->SetMarkerColor(kBlue);
    h7->SetLineWidth(0);
    h7->SetTitle(" ; ");
    h7->SetStats(0);
    h7->SetLineColor(kBlue);
    TLegend *l7 = new TLegend(0.65,0.80, 0.8, 0.907,"","brNDC");
    l7->AddEntry(h7, "US-LS", "pl");
    l7->SetFillStyle(0);
    l7->SetLineColor(0);
    l7->SetTextSize(0.04);
    l7->Draw("same");
    h7->GetYaxis()->SetTitle("Counts [-]");
    h7->GetYaxis()->SetTitleOffset(1.);
    h7->GetYaxis()->SetLabelSize(0.03);
    h7->GetYaxis()->SetTitleSize(0.04);
    h7->GetYaxis()->SetLabelFont(42);
    h7->GetYaxis()->SetTitleFont(42);
    //  h7->GetYaxis()->SetRangeUser(0,30000);
    h7->GetYaxis()->CenterTitle(kTRUE);
    h7->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h7->GetXaxis()->SetLabelFont(42);
    h7->GetXaxis()->SetTitleFont(42);
    //  h7->GetXaxis()->SetRangeUser(0,6);
    h7->GetXaxis()->SetTitleSize(0.04);
    h7->GetXaxis()->SetLabelSize(0.03);
    h7->GetXaxis()->CenterTitle(kTRUE);
    c7->SaveAs(Form("Mass_narrow12D.pdf"));
    c7->Close();


// Subtracted full window(1-2 GeV) -----------------------------------------------------------------------------------------------------------------------------
    TH1F* h8 = (TH1F*) d1 -> Get("rozdil12w;1");
    TCanvas *c8 = new TCanvas("c8", "Background-Signal(1-2 GeV)", 1400, 1100);
    h8->Draw("PE");
    c8->Update();
    h8->SetMarkerStyle(8);
    h8->SetMarkerColor(kBlue);
    h8->SetLineWidth(0);
    h8->SetTitle(" ; ");
    h8->SetStats(0);
    h8->SetLineColor(kBlue);
    TLegend *l8 = new TLegend(0.65,0.80, 0.8, 0.907,"","brNDC");
    l8->AddEntry(h8, "US-LS", "pl");
    l8->SetFillStyle(0);
    l8->SetLineColor(0);
    l8->SetTextSize(0.04);
    l8->Draw("same");
    h8->GetYaxis()->SetTitle("Counts [-]");
    h8->GetYaxis()->SetTitleOffset(1.);
    h8->GetYaxis()->SetLabelSize(0.03);
    h8->GetYaxis()->SetTitleSize(0.04);
    h8->GetYaxis()->SetLabelFont(42);
    h8->GetYaxis()->SetTitleFont(42);
    //  h8->GetYaxis()->SetRangeUser(0,30000);
    h8->GetYaxis()->CenterTitle(kTRUE);
    h8->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
    h8->GetXaxis()->SetLabelFont(42);
    h8->GetXaxis()->SetTitleFont(42);
    //  h8->GetXaxis()->SetRangeUser(0,6);
    h8->GetXaxis()->SetTitleSize(0.04);
    h8->GetXaxis()->SetLabelSize(0.03);
    h8->GetXaxis()->CenterTitle(kTRUE);
    c8->SaveAs(Form("Mass_full12D.pdf"));
    c8->Close();



}
