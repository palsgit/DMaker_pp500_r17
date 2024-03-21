#include <exception>
#include <assert.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <string>

#include "TString.h"
#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TFile.h"
#include "TArrayI.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TBranch.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TKey.h"
#include "TProfile.h"
#include "TGaxis.h"
#include "TList.h"

//Int_t padw = 800, padh = 600;

#include "style.C"
#include "DrawHistFileReco.h"
#include "HistDrawOptReco.h"

using namespace std;




/////////////////////////////////////////////////////////////
Int_t drawRecoHist() {

	style();
	//gStyle->SetOptStat(0);
	gROOT->ForceStyle();
	//gStyle->SetLineColor(kBlue);
	//gStyle->SetHistLineColor(kBlue);
	gStyle->SetHistLineColor(kBlack);
	gStyle->SetHistLineWidth(2);
	gStyle->SetOptTitle(1);
	gStyle->SetOptStat(1);

	TGaxis::SetExponentOffset(0.02, -0.07, "x");


	vector<TString> *list = new vector<TString>;


	list->push_back("h_nSigmaOneOverBetaK_tr");
	list->push_back("h_nSigmaOneOverBetaPi_tr");
	list->push_back("h_nSigmadEdxK_tr");
	list->push_back("h_nSigmadEdxPi_tr");

	list->push_back("hDedx_tr");
	list->push_back("hBetavsP_tr");


	TFile *file = new TFile("Projections.root", "read");
	//TFile *file = new TFile("data/output_bhcal_pion_test.root", "read");

	//TList *list = (TList*)file->Get("picoD0V2AnaMaker");
	//TH1F *hEventStat0 = (TH1F*)list->FindObject("hEventStat0");


	//TCanvas *cnv = new TCanvas();
	//cnv->cd();

	TString outputdir = "nsigmaPID";
	//outputdir = "outputReco_proton";
	//outputdir = "outputReco_pi-";

	gSystem->mkdir(outputdir);
	gSystem->cd(outputdir);
	//cnv->cd()->SaveAs("hEventStat0.png");
	gSystem->cd("../");

	//delete cnv;
	file->Close();

	drawAny("nsigmaPID/", "Projections.root", list);


	return 1.0;

}
