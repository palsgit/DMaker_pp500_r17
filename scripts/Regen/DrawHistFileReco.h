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


Double_t nSigmaTOFPionTrack_upper(Double_t x) { return 2.5; }
Double_t nSigmaTOFPionTrack_lower(Double_t x) { return -2.5; }

Double_t nSigmaTPCPion_upper(Double_t x) { return 3.0; }
Double_t nSigmaTPCPion_lower(Double_t x) { return -3.0; }

Double_t nSigmaTPCKaon_upper(Double_t x) { return 3.0; }
Double_t nSigmaTPCKaon_lower(Double_t x) { return -2.5; }

Double_t nSigmaTOFKaonTrack_upper(Double_t x) { 
	double f_nsigmaTOFKaonTrack_res =  1.21 + (0.10/(pow((x - 0.17), 1.21)));
    double f_nsigmaTOFKaonTrack_pos =  -0.01 + (0.02/(pow((x - 0.17), 1.73))); 

    /*float mSigma = 3.18;
    if (x > 2.0) mSigma = 2.41;*/
    /////float mSigma = 2.5; //weighted average
    float mSigma = 2.4;//midstep
	double kaon_higher = mSigma*f_nsigmaTOFKaonTrack_res + f_nsigmaTOFKaonTrack_pos;
    return kaon_higher; 
	}
Double_t nSigmaTOFKaonTrack_lower(Double_t x) {
	double f_nsigmaTOFKaonTrack_res =  1.21 + (0.10/(pow((x - 0.17), 1.21)));
    double f_nsigmaTOFKaonTrack_pos =  -0.01 + (0.02/(pow((x - 0.17), 1.73))); 
    /*float mSigma = -2.31;
    if (x > 1.25 && x < 1.50) mSigma = -1.56;
    if (x > 1.50 && x < 1.65) mSigma = -1.18;
    if (x > 1.65) mSigma = -0.8;*/
    float mSigma = -1.2; //midstep
    /////float mSigma = -1.0; //weighted average
	 double kaon_lower = mSigma*f_nsigmaTOFKaonTrack_res + f_nsigmaTOFKaonTrack_pos;
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
			TF1 *fitsigmakaonupper = new TF1("fitsigmakaonupper","nSigmaTOFKaonTrack_upper(x)",0.20,3.5);
            TF1 *fitsigmakaonlower = new TF1("fitsigmakaonlower","nSigmaTOFKaonTrack_lower(x)",0.20,3.5);

			TF1 *fitsigmapionupper = new TF1("fitsigmapionupper","nSigmaTOFPionTrack_upper(x)",0.20,3.5);
            TF1 *fitsigmapionlower = new TF1("fitsigmapionlower","nSigmaTOFPionTrack_lower(x)",0.20,3.5);

			TF1 *fitsigmaPionTrackupper = new TF1("fitsigmaPionTrackupper","nSigmaTPCPion_upper(x)",0.20,3.5);
            TF1 *fitsigmaPionTracklower = new TF1("fitsigmaPionTracklower","nSigmaTPCPion_lower(x)",0.20,3.5);

			TF1 *fitsigmaKaonTrackupper = new TF1("fitsigmaKaonTrackupper","nSigmaTPCKaon_upper(x)",0.20,3.5);
            TF1 *fitsigmaKaonTracklower = new TF1("fitsigmaKaonTracklower","nSigmaTPCKaon_lower(x)",0.20,3.5);


			TF1 *fitsigmaupper = new TF1();
			TF1 *fitsigmalower = new TF1();

			if (name == "h_nSigmaOneOverBetaPi_tr") {
                    fitsigmaupper = fitsigmapionupper;
					fitsigmalower = fitsigmapionlower;
					fitsigmaupper->Draw("same");
                    fitsigmalower->Draw("same");
				}

			if (name == "h_nSigmaOneOverBetaK_tr") {
                    fitsigmaupper = fitsigmakaonupper;
					fitsigmalower = fitsigmakaonlower;
					/////fitsigmaupper->Draw("same");
                    fitsigmalower->Draw("same");
				}

				if (name == "h_nSigmadEdxK_tr") {
                    fitsigmaupper = fitsigmaPionTrackupper;
					fitsigmalower = fitsigmaPionTracklower;
					fitsigmaupper->Draw("same");
                    fitsigmalower->Draw("same");
				}

			if (name == "h_nSigmadEdxPi_tr") {
                    fitsigmaupper = fitsigmaKaonTrackupper;
					fitsigmalower = fitsigmaKaonTracklower;
					fitsigmaupper->Draw("same");
                    fitsigmalower->Draw("same");
				}
				
    /////h_QA_OneOverBetaDiffKaon->Draw("colz");
    
			//save(cnv, names, 0, "png");
			cnv->cd()->SaveAs(Form("%s.png", name.Data()));

			if(logScale)
			{
				if(!class_name.CompareTo("TH1D")||!class_name.CompareTo("TH1F"))	cnv_log->cd()->SetLogy();
				//if(!class_name.CompareTo("TH1D")||!class_name.CompareTo("TH1F"))	cnv_log->cd()->SetLogx();
				if(!class_name.CompareTo("TH2D")||!class_name.CompareTo("TH2F"))	cnv_log->cd()->SetLogz();
				//if(!class_name.CompareTo("TH2D")||!class_name.CompareTo("TH2F"))	cnv_log->cd()->SetLogy();
				obj->Draw(drawOpt);
				TF1 *fitsigmakaonupper = new TF1("fitsigmakaonupper","nSigmaTOFKaonTrack_upper(x)",0.20,3.5);
                TF1 *fitsigmakaonlower = new TF1("fitsigmakaonlower","nSigmaTOFKaonTrack_lower(x)",0.20,3.5);
      
	            TF1 *fitsigmapionupper = new TF1("fitsigmapionupper","nSigmaTOFPionTrack_upper(x)",0.20,3.5);
                TF1 *fitsigmapionlower = new TF1("fitsigmapionlower","nSigmaTOFPionTrack_lower(x)",0.20,3.5);

			    TF1 *fitsigmaPionTrackupper = new TF1("fitsigmaPionTrackupper","nSigmaTPCPion_upper(x)",0.20,3.5);
                TF1 *fitsigmaPionTracklower = new TF1("fitsigmaPionTracklower","nSigmaTPCPion_lower(x)",0.20,3.5);

			    TF1 *fitsigmaKaonTrackupper = new TF1("fitsigmaKaonTrackupper","nSigmaTPCKaon_upper(x)",0.20,3.5);
                TF1 *fitsigmaKaonTracklower = new TF1("fitsigmaKaonTracklower","nSigmaTPCKaon_lower(x)",0.20,3.5);


			    TF1 *fitsigmaupper = new TF1();
			    TF1 *fitsigmalower = new TF1();

				if (name == "h_nSigmaOneOverBetaPi_tr") {
                    fitsigmaupper = fitsigmapionupper;
					fitsigmalower = fitsigmapionlower;
					fitsigmaupper->Draw("same");
                    fitsigmalower->Draw("same");
				}

			   if (name == "h_nSigmaOneOverBetaK_tr") {
                    fitsigmaupper = fitsigmakaonupper;
					fitsigmalower = fitsigmakaonlower;
					fitsigmaupper->Draw("same");
                    fitsigmalower->Draw("same");
				}

				if (name == "h_nSigmadEdxK_tr") {
                    fitsigmaupper = fitsigmaKaonTrackupper;
					fitsigmalower = fitsigmaKaonTracklower;
					fitsigmaupper->Draw("same");
                    fitsigmalower->Draw("same");
				}

			   if (name == "h_nSigmadEdxPi_tr") {
                    fitsigmaupper = fitsigmaPionTrackupper;
					fitsigmalower = fitsigmaPionTracklower;
					fitsigmaupper->Draw("same");
                    fitsigmalower->Draw("same");
				}

				

    /////h_QA_OneOverBetaDiffKaon->Draw("colz");
    
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
