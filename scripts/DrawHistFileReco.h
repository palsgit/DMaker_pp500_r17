/*
 * DrawHistFileReco.h
 *
 *  Created on: 2 mar 2022
 *      Author: Khaless
 */

#ifndef DRAWHISTFILERECO_H_
#define DRAWHISTFILERECO_H_

#include "TString.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMath.h"
#include "TKey.h"
#include "TProfile.h"
#include "TGaxis.h"
#include "TList.h"

#include "HistDrawOptReco.h"

bool drawFromList(TString name, TString *list, int max = 0);
bool drawFromVector(TString name, vector<TString> *vec);
int drawAny(TString dir, TString fname, TString *list, int max = 0);
int drawAny(TString dir, TString fname, vector<TString> *list);
int drawAnySubDir(TString dir, TString subdir, TString fname, TString *list, int max = 0);
int drawWithProfile(TString dir, TString fname, TString *list, int max = 0);
int drawProfile(TString dir, TString fname, TString *list, int max = 0);

Int_t padw = 1100, padh = 600;


Double_t nSigmaTOFPionTPC_upper(double x) { return 2.2; }
Double_t nSigmaTOFPionTPC_lower(double x) { return -2.0; }

Double_t nSigmaTOFKaonTPC_upper(double x) { 
	double f_nsigmaTOFKaonTPC_res =  1.07078 + 0.0194723/pow((x + 0.0183138), 3.41348);
    double f_nsigmaTOFKaonTPC_pos = 0.0283416 + 0.00230815/pow((x + 0.0270787), 4.96298); 
	double kaon_higher = 1.8*f_nsigmaTOFKaonTPC_res + f_nsigmaTOFKaonTPC_pos;
    return kaon_higher; 
	}
Double_t nSigmaTOFKaonTPC_lower(double x) {
	double f_nsigmaTOFKaonTPC_res =  1.07078 + 0.0194723/pow((x + 0.0183138), 3.41348);
    double f_nsigmaTOFKaonTPC_pos = 0.0283416 + 0.00230815/pow((x + 0.0270787), 4.96298); 
	 double kaon_lower = -1.2*f_nsigmaTOFKaonTPC_res + f_nsigmaTOFKaonTPC_pos;
	 return kaon_lower; 
	 }

/////////////////////////////////////////////////////////////
bool drawFromList(TString name, TString *list, int max)
{

	for (int i = 0; i < max; ++i) {
		if(name == list[i]) return true;
	}
	return false;

}
/////////////////////////////////////////////////////////////
bool drawFromVector(TString name, vector<TString> *vec)
{
	int max = vec->size();

	for (int i = 0; i < max; ++i) {
		if(name == vec->at(i)) return true;
	}
	return false;

}
/////////////////////////////////////////////////////////////
int drawAny(TString dir, TString fname, TString *list, int max)
{
	//Int_t padw = 600, padh = 400;

	TFile *file = new TFile(fname, "READ");
	TIter nextkey(file->GetListOfKeys());
	TKey *key = NULL;
	gSystem->cd(dir);

	TCanvas *cnv = new TCanvas("cnv", "cnv", padw, padh);
	TCanvas *cnv_log = new TCanvas("cnv_log", "cnv_log", padw, padh);

	while (key = (TKey*)nextkey()) {

		TObject *obj = dynamic_cast<TObject*>(key->ReadObj());

		cnv_log->cd()->SetLogy(0);
		cnv_log->cd()->SetLogz(0);

		if(obj != NULL)
		{
			TString drawOpt = "colz";
			Bool_t logScale = kTRUE;

			TString class_name = obj->ClassName();

			TString name = obj->GetName();
			TString name_log = "log_"+name;
			TString names[] = {name};
			TString names_log[] = {name_log};

			//cout<<name<<endl;

			if(!drawFromList(name, list, max))	continue;

			drawOpt = getHistDrawOpt(obj, logScale);

			cout<<drawOpt<<endl;

			cnv->cd();
			obj->Draw(drawOpt);
			//save(cnv, names, 0, "png");
			cnv->cd()->SaveAs(Form("%s.png", name.Data()));

			if(logScale)
			{
				if(!class_name.CompareTo("TH1D")||!class_name.CompareTo("TH1F"))	cnv_log->cd()->SetLogy();
				//if(!class_name.CompareTo("TH1D")||!class_name.CompareTo("TH1F"))	cnv_log->cd()->SetLogx();
				if(!class_name.CompareTo("TH2D")||!class_name.CompareTo("TH2F"))	cnv_log->cd()->SetLogz();
				//if(!class_name.CompareTo("TH2D")||!class_name.CompareTo("TH2F"))	cnv_log->cd()->SetLogy();
				obj->Draw(drawOpt);
				//save(cnv_log, names_log, 0, "png");
				cnv_log->cd()->SaveAs(Form("%s_log.png", name.Data()));
			}
		}
	}

	gSystem->cd("../");

	return 1;
}
/////////////////////////////////////////////////////////////
int drawAny(TString dir, TString fname, vector<TString> *list)
{
	//Int_t padw = 600, padh = 400;

	TFile *file = new TFile(fname, "READ");
	TIter nextkey(file->GetListOfKeys());
	TKey *key = NULL;
	gSystem->cd(dir);

	TCanvas *cnv = new TCanvas("cnv", "cnv", padw, padh);
	TCanvas *cnv_log = new TCanvas("cnv_log", "cnv_log", padw, padh);

	while (key = (TKey*)nextkey()) {

		TObject *obj = dynamic_cast<TObject*>(key->ReadObj());

		cnv_log->cd()->SetLogy(0);
		cnv_log->cd()->SetLogz(0);

		if(obj != NULL)
		{
			TString drawOpt = "colz";
			Bool_t logScale = kTRUE;

			TString class_name = obj->ClassName();

			TString name = obj->GetName();
			TString name_log = "log_"+name;
			TString names[] = {name};
			TString names_log[] = {name_log};

			//cout<<name<<endl;

			if(!drawFromVector(name, list))	continue;

			drawOpt = getHistDrawOpt(obj, logScale);

			cout<<drawOpt<<endl;

			cnv->cd();
			obj->Draw(drawOpt);
			TF1 *fitsigmakaonupper = new TF1("fitsigmakaonupper","nSigmaTOFKaonTPC_upper(x)",0.16,3.5);
            TF1 *fitsigmakaonlower = new TF1("fitsigmakaonlower","nSigmaTOFKaonTPC_lower(x)",0.16,3.5);

			TF1 *fitsigmapionupper = new TF1("fitsigmapionupper","nSigmaTOFPionTPC_upper(x)",0.16,3.5);
            TF1 *fitsigmapionlower = new TF1("fitsigmapionlower","nSigmaTOFPionTPC_lower(x)",0.16,3.5);

			TF1 *fitsigmaupper = new TF1();
			TF1 *fitsigmalower = new TF1();

			if (name == "h_OneOverBetaDiffPi_TPC") {
                    fitsigmaupper = fitsigmapionupper;
					fitsigmalower = fitsigmapionlower;
				}

			if (name == "h_OneOverBetaDiffK_TPC") {
                    fitsigmaupper = fitsigmakaonupper;
					fitsigmalower = fitsigmakaonlower;
				}
				
    /////h_QA_OneOverBetaDiffKaon->Draw("colz");
    fitsigmaupper->Draw("same");
    fitsigmalower->Draw("same");
			//save(cnv, names, 0, "png");
			cnv->cd()->SaveAs(Form("%s.png", name.Data()));

			if(logScale)
			{
				if(!class_name.CompareTo("TH1D")||!class_name.CompareTo("TH1F"))	cnv_log->cd()->SetLogy();
				//if(!class_name.CompareTo("TH1D")||!class_name.CompareTo("TH1F"))	cnv_log->cd()->SetLogx();
				if(!class_name.CompareTo("TH2D")||!class_name.CompareTo("TH2F"))	cnv_log->cd()->SetLogz();
				//if(!class_name.CompareTo("TH2D")||!class_name.CompareTo("TH2F"))	cnv_log->cd()->SetLogy();
				obj->Draw(drawOpt);
				TF1 *fitsigmakaonupper = new TF1("fitsigmakaonupper","nSigmaTOFKaonTPC_upper(x)",0.16,3.5);
                TF1 *fitsigmakaonlower = new TF1("fitsigmakaonlower","nSigmaTOFKaonTPC_lower(x)",0.16,3.5);
      
	            TF1 *fitsigmapionupper = new TF1("fitsigmapionupper","nSigmaTOFPionTPC_upper(x)",0.16,3.5);
                TF1 *fitsigmapionlower = new TF1("fitsigmapionlower","nSigmaTOFPionTPC_lower(x)",0.16,3.5);


				TF1 *fitsigmaupper = new TF1();
			    TF1 *fitsigmalower = new TF1();

				if (name == "h_OneOverBetaDiffPi_TPC") {
                    fitsigmaupper = fitsigmapionupper;
					fitsigmalower = fitsigmapionlower;
				}

				if (name == "h_OneOverBetaDiffK_TPC") {
                    fitsigmaupper = fitsigmakaonupper;
					fitsigmalower = fitsigmakaonlower;
				}

    /////h_QA_OneOverBetaDiffKaon->Draw("colz");
    fitsigmaupper->Draw("same");
    fitsigmalower->Draw("same");
				//save(cnv_log, names_log, 0, "png");
				cnv_log->cd()->SaveAs(Form("%s_log.png", name.Data()));
			}
		}
	}

	gSystem->cd("../");

	return 1;
}
/////////////////////////////////////////////////////////////
int drawAnySubDir(TString dir, TString subdir, TString fname, TString *list, int max)
{
	//Int_t padw = 600, padh = 400;

	TFile *file = new TFile(fname, "READ");
	TList *listdir = dynamic_cast<TList*>(file->Get(subdir));
	TIter nextkey(listdir->MakeIterator());
	TKey *key = NULL;
	gSystem->cd(dir);

	TCanvas *cnv = new TCanvas("cnv", "cnv", padw, padh);
	TCanvas *cnv_log = new TCanvas("cnv_log", "cnv_log", padw, padh);

	while (key = (TKey*)nextkey()) {

		TObject *obj = dynamic_cast<TObject*>(key->ReadObj());

		cnv_log->cd()->SetLogy(0);
		cnv_log->cd()->SetLogz(0);

		if(obj != NULL)
		{
			TString drawOpt = "colz";
			Bool_t logScale = kTRUE;

			TString class_name = obj->ClassName();

			TString name = obj->GetName();
			TString name_log = "log_"+name;
			TString names[] = {name};
			TString names_log[] = {name_log};

			//cout<<name<<endl;
/*
			if(name == subdir) key->cd();
			if(name != subdir) continue;*/

			//if(!drawFromList(name, list, max))	continue;

			drawOpt = getHistDrawOpt(obj, logScale);

			cout<<drawOpt<<endl;

			cnv->cd();
			obj->Draw(drawOpt);
			//save(cnv, names, 0, "png");
			cnv->cd()->SaveAs(Form("%s.png", name.Data()));

			if(logScale)
			{
				if(!class_name.CompareTo("TH1D")||!class_name.CompareTo("TH1F"))	cnv_log->cd()->SetLogy();
				if(!class_name.CompareTo("TH2D")||!class_name.CompareTo("TH2F"))	cnv_log->cd()->SetLogz();
				obj->Draw(drawOpt);
				//save(cnv_log, names_log, 0, "png");
				cnv_log->cd()->SaveAs(Form("%s_log.png", name.Data()));
			}
		}
	}

	gSystem->cd("../");

	return 1;
}
/////////////////////////////////////////////////////////////
int drawWithProfile(TString dir, TString fname, TString *list, int max)
{
	//Int_t padw = 600, padh = 400;

	TFile *file = new TFile(fname, "READ");
	TIter nextkey(file->GetListOfKeys());
	TKey *key = NULL;
	gSystem->cd(dir);

	TCanvas *cnv = new TCanvas("cnv", "cnv", padw, padh);
	TCanvas *cnv_log = new TCanvas("cnv_log", "cnv_log", padw, padh);

	while (key = (TKey*)nextkey()) {

		TObject *obj = dynamic_cast<TObject*>(key->ReadObj());

		if(obj != NULL)
		{
			TString drawOpt = "colz";
			Bool_t logScale = kTRUE;

			TString class_name = obj->ClassName();

			TString name = obj->GetName();
			TString name_log = "log_"+name;
			TString names[] = {name};
			TString names_log[] = {name_log};

			if(!drawFromList(name, list, max))	continue;

			drawOpt = getHistDrawOpt(obj, logScale);

			if(!class_name.CompareTo("TProfile")) continue;

			TProfile *profile = (TProfile *)file->Get(name+"_pfx");
			profile->SetLineColor(kBlack);

			cout<<drawOpt<<endl;

			cnv->cd();
			obj->Draw(drawOpt);
			profile->Draw("samee");
			//save(cnv, names, 0, "png");
			cnv->cd()->SaveAs(Form("%s.png", name.Data()));

			if(logScale)
			{
				if(!class_name.CompareTo("TH1D")||!class_name.CompareTo("TH1F"))	cnv_log->cd()->SetLogy();
				if(!class_name.CompareTo("TH2D")||!class_name.CompareTo("TH2F"))	cnv_log->cd()->SetLogz();
				obj->Draw(drawOpt);
				profile->Draw("samee");
				//cnv_log->cd()->SaveAs();
				//save(cnv_log, names_log, 0, "png");
				cnv_log->cd()->SaveAs(Form("%s_log.png", name.Data()));
			}
		}
	}

	gSystem->cd("../");

	return 1;
}

/////////////////////////////////////////////////////////////
int drawProfile(TString dir, TString fname, TString *list, int max)
{
	//Int_t padw = 600, padh = 400;

	TFile *file = new TFile(fname, "READ");
	TIter nextkey(file->GetListOfKeys());
	TKey *key = NULL;
	gSystem->cd(dir);

	TCanvas *cnv = new TCanvas("cnv", "cnv", padw, padh);
	TCanvas *cnv_log = new TCanvas("cnv_log", "cnv_log", padw, padh);

	while (key = (TKey*)nextkey()) {

		TObject *obj = dynamic_cast<TObject*>(key->ReadObj());

		if(obj != NULL)
		{
			TString drawOpt = "colz";
			Bool_t logScale = kTRUE;

			TString class_name = obj->ClassName();

			TString name = obj->GetName();
			TString name_log = "log_"+name;
			TString names[] = {name};
			TString names_log[] = {name_log};

			if(!drawFromList(name, list, max))	continue;

			drawOpt = getHistDrawOpt(obj, logScale);

			if(!class_name.CompareTo("TProfile")) continue;

			//TProfile *profile = (TProfile *)file->Get(name+"_pfx");
			//profile->SetLineColor(kBlack);

			cout<<drawOpt<<endl;

			cnv->cd();
			obj->Draw(drawOpt);
			//profile->Draw("samee");
			//save(cnv, names, 0, "png");
			cnv->cd()->SaveAs(Form("%s.png", name.Data()));

			if(logScale)
			{
				if(!class_name.CompareTo("TH1D")||!class_name.CompareTo("TH1F"))	cnv_log->cd()->SetLogy();
				if(!class_name.CompareTo("TH2D")||!class_name.CompareTo("TH2F"))	cnv_log->cd()->SetLogz();
				obj->Draw(drawOpt);
				//profile->Draw("samee");
				//cnv_log->cd()->SaveAs();
				//save(cnv_log, names_log, 0, "png");
				cnv_log->cd()->SaveAs(Form("%s_log.png", name.Data()));
			}
		}
	}

	gSystem->cd("../");

	return 1;
}


#endif /* DRAWHISTFILERECO_H_ */
