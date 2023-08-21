/*
 * HistDrawOpt.h
 *
 *  Created on: 26 sie 2015
 *      Author: Khaless
 */

#include "TString.h"

#ifndef HISTDRAWOPTRECO_H_
#define HISTDRAWOPTRECO_H_

TString getHistDrawOpt(TObject *obj, Bool_t &logScale)
{
	TString drawOpt;

	TString name = obj->GetName();

	//cout<<name<<endl;

	

	if(!name.CompareTo("h_OneOverBetaDiffK_TPC")) {
		TH2 *objH2 = dynamic_cast<TH2*>(obj);
		drawOpt = "colz";
		objH2->GetXaxis()->SetRangeUser(0.0, 3.5);
		objH2->GetYaxis()->SetRangeUser(-10.0, 10.0);
		objH2->SetStats(kFALSE);
		//drawOpt = "e";
		//logScale=kTRUE;
	}

	if(!name.CompareTo("h_OneOverBetaDiffPi_TPC")) {
		TH2 *objH2 = dynamic_cast<TH2*>(obj);
		drawOpt = "colz";
		objH2->GetXaxis()->SetRangeUser(0.0, 3.5);
		objH2->GetYaxis()->SetRangeUser(-10.0, 10.0);
		objH2->SetStats(kFALSE);
		//drawOpt = "e";
		//logScale=kTRUE;
	}

	



/*
	if(!name.CompareTo("NtracksBvsF")) {
		TH2 *objH2 = dynamic_cast<TH2*>(obj);
		drawOpt = "colz";
		logScale=kFALSE;
		//objH2->GetXaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
		//objH2->GetYaxis()->SetRangeUser(-TMath::Pi(), TMath::Pi());
		objH2->SetTitle("Number of tracks backward vs. forward");
		objH2->GetXaxis()->SetTitle("Ntracks backward");
		objH2->GetYaxis()->SetTitle("Ntracks forward");
	}
*/



	return drawOpt;
}



#endif /* HISTDRAWOPTRECO_H_ */
