#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TString.h"
#include "TROOT.h"
#include "TChain.h"
#include "TNtuple.h"

using namespace std;


void EventStatistics()
{

    TFile *f1 = TFile::Open("1merged_output_3_000.root");
    TList *list1 = (TList*)f1->Get("picoD0AnaMaker;1");
    TH1F *hEventStat1 = (TH1F*)list1->FindObject("hEventStat1");
    TH2F *BetaPion1 = (TH2F*)list1->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon1 = (TH2F*)list1->FindObject("h_QA_OneOverBetaDiffKaon");


    TFile *f2 = TFile::Open("1merged_output_3_001.root");
    TList *list2 = (TList*)f2->Get("picoD0AnaMaker;1");
    TH1F *hEventStat2 = (TH1F*)list2->FindObject("hEventStat1");
    TH2F *BetaPion2 = (TH2F*)list2->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon2 = (TH2F*)list2->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f3 = TFile::Open("1merged_output_3_002.root");
    TList *list3 = (TList*)f3->Get("picoD0AnaMaker;1");
    TH1F *hEventStat3 = (TH1F*)list3->FindObject("hEventStat1");
    TH2F *BetaPion3 = (TH2F*)list3->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon3 = (TH2F*)list3->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f4 = TFile::Open("1merged_output_3_003.root");
    TList *list4 = (TList*)f4->Get("picoD0AnaMaker;1");
    TH1F *hEventStat4 = (TH1F*)list4->FindObject("hEventStat1");
    TH2F *BetaPion4 = (TH2F*)list4->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon4 = (TH2F*)list4->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f5 = TFile::Open("1merged_output_3_004.root");
    TList *list5 = (TList*)f5->Get("picoD0AnaMaker;1");
    TH1F *hEventStat5 = (TH1F*)list5->FindObject("hEventStat1");
    TH2F *BetaPion5 = (TH2F*)list5->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon5 = (TH2F*)list5->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f6 = TFile::Open("1merged_output_3_005.root");
    TList *list6 = (TList*)f6->Get("picoD0AnaMaker;1");
    TH1F *hEventStat6 = (TH1F*)list6->FindObject("hEventStat1");
    TH2F *BetaPion6 = (TH2F*)list6->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon6 = (TH2F*)list6->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f7 = TFile::Open("1merged_output_3_006.root");
    TList *list7 = (TList*)f7->Get("picoD0AnaMaker;1");
    TH1F *hEventStat7 = (TH1F*)list7->FindObject("hEventStat1");
    TH2F *BetaPion7 = (TH2F*)list7->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon7 = (TH2F*)list7->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f8 = TFile::Open("2merged_output_3_000.root");
    TList *list8 = (TList*)f8->Get("picoD0AnaMaker;1");
    TH1F *hEventStat8 = (TH1F*)list8->FindObject("hEventStat1");
    TH2F *BetaPion8 = (TH2F*)list8->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon8 = (TH2F*)list8->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f9 = TFile::Open("2merged_output_3_001.root");
    TList *list9 = (TList*)f9->Get("picoD0AnaMaker;1");
    TH1F *hEventStat9 = (TH1F*)list9->FindObject("hEventStat1");
    TH2F *BetaPion9 = (TH2F*)list9->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon9 = (TH2F*)list9->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f10 = TFile::Open("2merged_output_3_002.root");
    TList *list10 = (TList*)f10->Get("picoD0AnaMaker;1");
    TH1F *hEventStat10 = (TH1F*)list10->FindObject("hEventStat1");
    TH2F *BetaPion10 = (TH2F*)list10->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon10 = (TH2F*)list10->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f11 = TFile::Open("2merged_output_3_003.root");
    TList *list11 = (TList*)f11->Get("picoD0AnaMaker;1");
    TH1F *hEventStat11 = (TH1F*)list11->FindObject("hEventStat1");
    TH2F *BetaPion11 = (TH2F*)list11->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon11 = (TH2F*)list11->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f12 = TFile::Open("2merged_output_3_004.root");
    TList *list12 = (TList*)f12->Get("picoD0AnaMaker;1");
    TH1F *hEventStat12 = (TH1F*)list12->FindObject("hEventStat1");
    TH2F *BetaPion12 = (TH2F*)list12->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon12 = (TH2F*)list12->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f13 = TFile::Open("2merged_output_3_005.root");
    TList *list13 = (TList*)f13->Get("picoD0AnaMaker;1");
    TH1F *hEventStat13 = (TH1F*)list13->FindObject("hEventStat1");
    TH2F *BetaPion13 = (TH2F*)list13->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon13 = (TH2F*)list13->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f14 = TFile::Open("2merged_output_3_006.root");
    TList *list14 = (TList*)f14->Get("picoD0AnaMaker;1");
    TH1F *hEventStat14 = (TH1F*)list14->FindObject("hEventStat1");
    TH2F *BetaPion14 = (TH2F*)list14->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon14 = (TH2F*)list14->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f15 = TFile::Open("3merged_output_3_000.root");
    TList *list15 = (TList*)f15->Get("picoD0AnaMaker;1");
    TH1F *hEventStat15 = (TH1F*)list15->FindObject("hEventStat1");
    TH2F *BetaPion15 = (TH2F*)list15->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon15 = (TH2F*)list15->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f16 = TFile::Open("3merged_output_3_001.root");
    TList *list16 = (TList*)f16->Get("picoD0AnaMaker;1");
    TH1F *hEventStat16 = (TH1F*)list16->FindObject("hEventStat1");
    TH2F *BetaPion16 = (TH2F*)list16->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon16 = (TH2F*)list16->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f17 = TFile::Open("3merged_output_3_002.root");
    TList *list17 = (TList*)f17->Get("picoD0AnaMaker;1");
    TH1F *hEventStat17 = (TH1F*)list17->FindObject("hEventStat1");
    TH2F *BetaPion17 = (TH2F*)list17->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon17 = (TH2F*)list17->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f18 = TFile::Open("3merged_output_3_003.root");
    TList *list18 = (TList*)f18->Get("picoD0AnaMaker;1");
    TH1F *hEventStat18 = (TH1F*)list18->FindObject("hEventStat1");
    TH2F *BetaPion18 = (TH2F*)list18->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon18 = (TH2F*)list18->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f19 = TFile::Open("3merged_output_3_004.root");
    TList *list19 = (TList*)f19->Get("picoD0AnaMaker;1");
    TH1F *hEventStat19 = (TH1F*)list19->FindObject("hEventStat1");
    TH2F *BetaPion19 = (TH2F*)list19->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon19 = (TH2F*)list19->FindObject("h_QA_OneOverBetaDiffKaon");


    TFile *f20 = TFile::Open("1merged_output_3_007.root");
    TList *list20 = (TList*)f20->Get("picoD0AnaMaker;1");
    TH1F *hEventStat20 = (TH1F*)list20->FindObject("hEventStat1");
    TH2F *BetaPion20 = (TH2F*)list20->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon20 = (TH2F*)list20->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f21 = TFile::Open("2merged_output_3_007.root");
    TList *list21 = (TList*)f21->Get("picoD0AnaMaker;1");
    TH1F *hEventStat21 = (TH1F*)list21->FindObject("hEventStat1");
    TH2F *BetaPion21 = (TH2F*)list21->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon21 = (TH2F*)list21->FindObject("h_QA_OneOverBetaDiffKaon");

    TFile *f22 = TFile::Open("3merged_output_3_005.root");
    TList *list22 = (TList*)f22->Get("picoD0AnaMaker;1");
    TH1F *hEventStat22 = (TH1F*)list22->FindObject("hEventStat1");
    TH2F *BetaPion22 = (TH2F*)list22->FindObject("h_QA_OneOverBetaDiffPion");
    TH2F *BetaKaon22 = (TH2F*)list22->FindObject("h_QA_OneOverBetaDiffKaon");
    
    TList *list = new TList;
    list->Add(hEventStat1);
    list->Add(hEventStat2);
    list->Add(hEventStat3);
    list->Add(hEventStat4);
    list->Add(hEventStat5);
    list->Add(hEventStat6);
    list->Add(hEventStat7);
    list->Add(hEventStat8);
    list->Add(hEventStat9);
    list->Add(hEventStat10);
    list->Add(hEventStat11);
    list->Add(hEventStat12);
    list->Add(hEventStat13);
    list->Add(hEventStat14);
    list->Add(hEventStat15);
    list->Add(hEventStat16);
    list->Add(hEventStat17);
    list->Add(hEventStat18);
    list->Add(hEventStat19);
    list->Add(hEventStat20);
    list->Add(hEventStat21);
    list->Add(hEventStat22);

    TList *BetaPion = new TList;
    BetaPion->Add(BetaPion1);
    BetaPion->Add(BetaPion2);
    BetaPion->Add(BetaPion3);
    BetaPion->Add(BetaPion4);
    BetaPion->Add(BetaPion5);
    BetaPion->Add(BetaPion6);
    BetaPion->Add(BetaPion7);
    BetaPion->Add(BetaPion8);
    BetaPion->Add(BetaPion9);
    BetaPion->Add(BetaPion10);
    BetaPion->Add(BetaPion11);
    BetaPion->Add(BetaPion12);
    BetaPion->Add(BetaPion13);
    BetaPion->Add(BetaPion14);
    BetaPion->Add(BetaPion15);
    BetaPion->Add(BetaPion16);
    BetaPion->Add(BetaPion17);
    BetaPion->Add(BetaPion18);
    BetaPion->Add(BetaPion19);
    BetaPion->Add(BetaPion20);
    BetaPion->Add(BetaPion21);
    BetaPion->Add(BetaPion22);

    TList *BetaKaon = new TList;
    BetaKaon->Add(BetaKaon1);
    BetaKaon->Add(BetaKaon2);
    BetaKaon->Add(BetaKaon3);
    BetaKaon->Add(BetaKaon4);
    BetaKaon->Add(BetaKaon5);
    BetaKaon->Add(BetaKaon6);
    BetaKaon->Add(BetaKaon7);
    BetaKaon->Add(BetaKaon8);
    BetaKaon->Add(BetaKaon9);
    BetaKaon->Add(BetaKaon10);
    BetaKaon->Add(BetaKaon11);
    BetaKaon->Add(BetaKaon12);
    BetaKaon->Add(BetaKaon13);
    BetaKaon->Add(BetaKaon14);
    BetaKaon->Add(BetaKaon15);
    BetaKaon->Add(BetaKaon16);
    BetaKaon->Add(BetaKaon17);
    BetaKaon->Add(BetaKaon18);
    BetaKaon->Add(BetaKaon19);
    BetaKaon->Add(BetaKaon20);
    BetaKaon->Add(BetaKaon21);
    BetaKaon->Add(BetaKaon22);
    
    TCanvas *c1 = new TCanvas("c1", "Event Statistics", 1500, 1000);

//    TPaveLabel *t = new TPaveLabel(0.0, 0.9, 0.3, 1.0, "", "brNDC");
//    t->SetBorderSize(0);

    c1->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c1->Update();

    TH1F *hEventStatFull = (TH1F*)hEventStat1->Clone("hEventStat1");
    hEventStatFull->Reset();
    hEventStatFull->Merge(list);
    hEventStatFull->Draw();
    hEventStatFull->SetTitle(" ; ");
    hEventStatFull->SetStats(0);

    TLegend *l1 = new TLegend(0.55,0.65, 0.73, 0.89,"","brNDC");

    l1->AddEntry((TObject*)0, "", "");
    l1->AddEntry((TObject*)0, "THIS THESIS", "");
    l1->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 510 GeV", "");

    l1->SetFillStyle(0);
    l1->SetLineColor(0);
    l1->SetTextSize(0.03);
    l1->Draw("same");

    c1->SaveAs(Form("pp500_EventStatistics.pdf"));
    c1->Close();



    TCanvas *c2 = new TCanvas("c2", "Beta Pion", 1500, 1000);

//    TPaveLabel *t = new TPaveLabel(0.0, 0.9, 0.3, 1.0, "", "brNDC");
//    t->SetBorderSize(0);

    c2->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c2->Update();

    TH2F *OneOverBetaPion = (TH2F*)BetaPion1->Clone("h_QA_OneOverBetaDiffPion");
    OneOverBetaPion->Reset();
    OneOverBetaPion->Merge(BetaPion);
    OneOverBetaPion->Draw("colz");
    OneOverBetaPion->SetTitle(" ; ");
    OneOverBetaPion->SetStats(0);

    TLegend *l2 = new TLegend(0.55,0.65, 0.73, 0.89,"","brNDC");

    l2->AddEntry((TObject*)0, "", "");
    l2->AddEntry((TObject*)0, "THIS THESIS", "");
    l2->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 500 GeV", "");

    l2->SetFillStyle(0);
    l2->SetLineColor(0);
    l2->SetTextSize(0.03);
    l2->Draw("same");

    c2->SaveAs(Form("pp500_OneOverBetaPion.pdf"));
    c2->Close();


    TCanvas *c3 = new TCanvas("c3", "Beta Kaon", 1500, 1000);

//    TPaveLabel *t = new TPaveLabel(0.0, 0.9, 0.3, 1.0, "", "brNDC");
//    t->SetBorderSize(0);

    c3->SetGrid(0,0);
    gStyle->SetOptTitle(kFALSE);
    gStyle->SetOptDate(0);
    c3->Update();

    TH2F *OneOverBetaKaon = (TH2F*)BetaKaon1->Clone("h_QA_OneOverBetaDiffKaon");
    OneOverBetaKaon->Reset();
    OneOverBetaKaon->Merge(BetaKaon);
    OneOverBetaKaon->Draw("colz");
    OneOverBetaKaon->SetTitle(" ; ");
    OneOverBetaKaon->SetStats(0);

    TLegend *l3 = new TLegend(0.55,0.65, 0.73, 0.89,"","brNDC");

    l3->AddEntry((TObject*)0, "", "");
    l3->AddEntry((TObject*)0, "THIS THESIS", "");
    l3->AddEntry((TObject*)0, "p+p #sqrt{#it{s}} = 500 GeV", "");

    l3->SetFillStyle(0);
    l3->SetLineColor(0);
    l3->SetTextSize(0.03);
    l3->Draw("same");

    c3->SaveAs(Form("pp500_OneOverBetaKaon.pdf"));
    c3->Close();










    cout<<"Michale jsi dobrej"<<endl;
    cout<<"Hotovo, Jarvis"<<endl;

}

