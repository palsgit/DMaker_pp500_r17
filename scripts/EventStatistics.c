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

    TFile *f2 = TFile::Open("1merged_output_3_001.root");
    TList *list2 = (TList*)f2->Get("picoD0AnaMaker;1");
    TH1F *hEventStat2 = (TH1F*)list2->FindObject("hEventStat1");

    TFile *f3 = TFile::Open("1merged_output_3_002.root");
    TList *list3 = (TList*)f3->Get("picoD0AnaMaker;1");
    TH1F *hEventStat3 = (TH1F*)list3->FindObject("hEventStat1");

    TFile *f4 = TFile::Open("1merged_output_3_003.root");
    TList *list4 = (TList*)f4->Get("picoD0AnaMaker;1");
    TH1F *hEventStat4 = (TH1F*)list4->FindObject("hEventStat1");

    TFile *f5 = TFile::Open("1merged_output_3_004.root");
    TList *list5 = (TList*)f5->Get("picoD0AnaMaker;1");
    TH1F *hEventStat5 = (TH1F*)list5->FindObject("hEventStat1");

    TFile *f6 = TFile::Open("1merged_output_3_005.root");
    TList *list6 = (TList*)f6->Get("picoD0AnaMaker;1");
    TH1F *hEventStat6 = (TH1F*)list6->FindObject("hEventStat1");

    TFile *f7 = TFile::Open("1merged_output_3_006.root");
    TList *list7 = (TList*)f7->Get("picoD0AnaMaker;1");
    TH1F *hEventStat7 = (TH1F*)list7->FindObject("hEventStat1");

    TFile *f8 = TFile::Open("2merged_output_3_000.root");
    TList *list8 = (TList*)f8->Get("picoD0AnaMaker;1");
    TH1F *hEventStat8 = (TH1F*)list8->FindObject("hEventStat1");

    TFile *f9 = TFile::Open("2merged_output_3_001.root");
    TList *list9 = (TList*)f9->Get("picoD0AnaMaker;1");
    TH1F *hEventStat9 = (TH1F*)list9->FindObject("hEventStat1");

    TFile *f10 = TFile::Open("2merged_output_3_002.root");
    TList *list10 = (TList*)f10->Get("picoD0AnaMaker;1");
    TH1F *hEventStat10 = (TH1F*)list10->FindObject("hEventStat1");

    TFile *f11 = TFile::Open("2merged_output_3_003.root");
    TList *list11 = (TList*)f11->Get("picoD0AnaMaker;1");
    TH1F *hEventStat11 = (TH1F*)list11->FindObject("hEventStat1");

    TFile *f12 = TFile::Open("2merged_output_3_004.root");
    TList *list12 = (TList*)f12->Get("picoD0AnaMaker;1");
    TH1F *hEventStat12 = (TH1F*)list12->FindObject("hEventStat1");

    TFile *f13 = TFile::Open("2merged_output_3_005.root");
    TList *list13 = (TList*)f13->Get("picoD0AnaMaker;1");
    TH1F *hEventStat13 = (TH1F*)list13->FindObject("hEventStat1");

    TFile *f14 = TFile::Open("2merged_output_3_006.root");
    TList *list14 = (TList*)f14->Get("picoD0AnaMaker;1");
    TH1F *hEventStat14 = (TH1F*)list14->FindObject("hEventStat1");

    TFile *f15 = TFile::Open("3merged_output_3_000.root");
    TList *list15 = (TList*)f15->Get("picoD0AnaMaker;1");
    TH1F *hEventStat15 = (TH1F*)list15->FindObject("hEventStat1");

    TFile *f16 = TFile::Open("3merged_output_3_001.root");
    TList *list16 = (TList*)f16->Get("picoD0AnaMaker;1");
    TH1F *hEventStat16 = (TH1F*)list16->FindObject("hEventStat1");

    TFile *f17 = TFile::Open("3merged_output_3_002.root");
    TList *list17 = (TList*)f17->Get("picoD0AnaMaker;1");
    TH1F *hEventStat17 = (TH1F*)list17->FindObject("hEventStat1");

    TFile *f18 = TFile::Open("3merged_output_3_003.root");
    TList *list18 = (TList*)f18->Get("picoD0AnaMaker;1");
    TH1F *hEventStat18 = (TH1F*)list18->FindObject("hEventStat1");

    TFile *f19 = TFile::Open("3merged_output_3_004.root");
    TList *list19 = (TList*)f19->Get("picoD0AnaMaker;1");
    TH1F *hEventStat19 = (TH1F*)list19->FindObject("hEventStat1");


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



    c1->SaveAs(Form("EventStatistics.pdf"));
    c1->Close();


    cout<<"Michale jsi dobrej"<<endl;
    cout<<"Hotovo, Jarvis"<<endl;

}

