Double_t bgfitfunction(Double_t *x, Double_t *par)
	{
		if (x[0] > 1.83 && x[0] < 1.89) {
			TF1::RejectPoint();
			return 0;
		}
		return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
//		return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
	}
Double_t linearbgfitfunction(Double_t *x, Double_t *par)
{
	if (x[0] > 1.8 && x[0] < 1.92) {
		TF1::RejectPoint();
		return 0;
	}
	//		return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] +
	//		par[4]*x[0]*x[0]*x[0]*x[0];
	return par[0] + par[1]*x[0];
}
Double_t parabola(Double_t *x, Double_t *par)
{
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
//	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

Double_t linear(Double_t *x, Double_t *par)
{
	//	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] +
	//	par[4]*x[0]*x[0]*x[0]*x[0];
	return par[0] + par[1]*x[0];
}

// Quadratic background function
Double_t background(Double_t *x, Double_t *par) {
	return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
}

// Lorenzian Peak function
Double_t lorentzianPeak(Double_t *x, Double_t *par) {
	return (0.5*par[0]*par[1]/TMath::Pi()) / 
    TMath::Max( 1.e-10,(x[0]-par[2])*(x[0]-par[2]) 
			   + .25*par[1]*par[1]);
}

// Gaussian Peak function
Double_t gaussianPeak(Double_t *x, Double_t *par) {
	return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2*par[2]*par[2]));
	
}
// Gaussian Peak with area as first parameter
Double_t normalDistribution(Double_t *x, Double_t *par) {
	return par[0]/(TMath::Sqrt(2*TMath::Pi()*par[2]*par[2]))*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2*par[2]*par[2]));
	
}

// Sum of background and peak function
Double_t fitFunction(Double_t *x, Double_t *par) {
	return parabola(x,par) + normalDistribution(x,&par[4]);
}
void PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val,
              Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle);

void RawYieldPt_og(float pT1, float pT2, int rebin) {
	
	float x1,x2,xs1,xs2,lowEdge,highEdge,normEdge1,normEdge2,ptBin1,ptBin2;
	int nBinsInvMass;
	bool GeometricMean = true;
	double par[10];
	char *parstring = new char[50];
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetFrameFillColor(10);
	gStyle->SetCanvasColor(10);
	gStyle->SetPadColor(10);
	// Fitting the backgound directly
//	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	TGaxis::SetMaxDigits(3);
	gStyle->SetTitleX(0.3);
	//TF1 *bgfitfunction = new TF1("bgfitfunction",bgfitfunction,1.7,2.15,5);
	int rebinMixed = 1;
	int rebinRotate = 1;
	int rebinLike = 1;
	ptBin1 = D->FindBin(pT1);
	ptBin2 = D->FindBin(pT2);
	D->ProjectionY("_py",ptBin1,ptBin2);
	Dbar->ProjectionY("_py",ptBin1,ptBin2);
	MixedD->ProjectionY("_py",ptBin1,ptBin2);
	MixedDbar->ProjectionY("_py",ptBin1,ptBin2);
	LikeBgD->ProjectionY("_py",ptBin1,ptBin2);
	LikeBgDbar->ProjectionY("_py",ptBin1,ptBin2);
	DRotate->ProjectionY("_py",ptBin1,ptBin2);
	DbarRotate->ProjectionY("_py",ptBin1,ptBin2);
  TH1D *D_py = (TH1D*)gDirectory->Get("D_py");
  TH1D *Dbar_py = (TH1D*)gDirectory->Get("Dbar_py");
  TH1D *MixedD_py = (TH1D*)gDirectory->Get("MixedD_py");
  TH1D *MixedDbar_py = (TH1D*)gDirectory->Get("MixedDbar_py");
  TH1D *LikeBgD_py = (TH1D*)gDirectory->Get("LikeBgD_py");
  TH1D *LikeBgDbar_py = (TH1D*)gDirectory->Get("LikeBgDbar_py");
  TH1D *DRotate_py = (TH1D*)gDirectory->Get("DRotate_py");
  TH1D *DbarRotate_py = (TH1D*)gDirectory->Get("DbarRotate_py");
  
  
  D_py->Sumw2();
	Dbar_py->Sumw2();
	MixedD_py->Sumw2();
	MixedDbar_py->Sumw2();
	LikeBgD_py->Sumw2();
	LikeBgDbar_py->Sumw2();
	DRotate_py->Sumw2();
	DbarRotate_py->Sumw2();
	
	D_py->Rebin(rebin);
	Dbar_py->Rebin(rebin);
	MixedD_py->Rebin(rebin);
	MixedDbar_py->Rebin(rebin);
	LikeBgD_py->Rebin(rebin);
	LikeBgDbar_py->Rebin(rebin);
	DRotate_py->Rebin(rebin);
	DbarRotate_py->Rebin(rebin);

	nBinsInvMass = D_py->GetNbinsX();
	lowEdge = D_py->GetBinLowEdge(1);
	highEdge = D_py->GetBinLowEdge(nBinsInvMass+1);
	
	cout<<nBinsInvMass<<"   "<<lowEdge<<"    "<<highEdge<<endl;
	
	TLegend *leg = new TLegend(0.6,0.7,0.88,0.88,"");
	leg->SetFillColor(10);
	leg->SetTextFont(42);
	TH1D *Unlike = new TH1D("Unlike","",nBinsInvMass,lowEdge,highEdge);
	TH1D *MixedBg = new TH1D("Mixed","",nBinsInvMass,lowEdge,highEdge);
	TH1D *LikeBg= new TH1D("LikeBg","",nBinsInvMass,lowEdge,highEdge);
	TH1D *RotateBg= new TH1D("RotateBg","",nBinsInvMass,lowEdge,highEdge);
	
	Unlike->Sumw2();
	MixedBg->Sumw2();
	LikeBg->Sumw2();
	RotateBg->Sumw2();
	
	Unlike->Add(D_py,Dbar_py,1,1);
	MixedBg->Add(MixedD_py,MixedDbar_py,1,1);
	if(GeometricMean){
		float a,b,da,db;
		for(int i=1; i<=nBinsInvMass; i++){
			a = LikeBgD_py->GetBinContent(i);
			b = LikeBgDbar_py->GetBinContent(i);
			da = LikeBgD_py->GetBinError(i);
			db = LikeBgDbar_py->GetBinError(i);	
			if(a==0 || b==0) 
			{
				LikeBg->SetBinContent(i,2*TMath::Sqrt(a*b));
				LikeBg->SetBinError(i,0);
			}else{
				LikeBg->SetBinContent(i,2*TMath::Sqrt(a*b));
				LikeBg->SetBinError(i,TMath::Sqrt(da*da*b/a+db*db*a/b));
			}
		}
	}else LikeBg->Add(LikeBgD_py,LikeBgDbar_py,1,1);	
	RotateBg->Add(DRotate_py,DbarRotate_py,1,1);
	
	Double_t Scalefactor;
	normEdge1 = Unlike->FindBin(1.7);
	normEdge2 = Unlike->FindBin(1.8);
//	Scalefactor = Unlike->Integral(normEdge1,normEdge2)/MixedBg->Integral(normEdge1,normEdge2);
//	MixedBg->Scale(Scalefactor);
//	cout<<"Mixed event Scale factor: "<<Scalefactor<<endl;
	Scalefactor = Unlike->Integral(normEdge1,normEdge2)/RotateBg->Integral(normEdge1,normEdge2); 
	RotateBg->Scale(Scalefactor);
  cout<<"Rotated momentum scale factor: "<<Scalefactor<<endl;
	Scalefactor = Unlike->Integral(normEdge1,normEdge2)/LikeBg->Integral(normEdge1,normEdge2);
	LikeBg->Scale(Scalefactor);
  cout<<"Like sign scale factor: "<<Scalefactor<<endl;

	//********************************************* Draw Invariant mass **********************************************
	TCanvas *c1 = new TCanvas("c1");
	c1->cd();
	c1->SetRightMargin(0.05);
	Unlike->GetXaxis()->SetTitle("M_{K#pi} [GeV/c^{2}]");
	Unlike->GetYaxis()->SetTitle("Raw Yield (/0.01 GeV/c^{2})");
	Unlike->SetMarkerStyle(4);
	Unlike->GetYaxis()->SetTitleFont(42);
	Unlike->GetYaxis()->SetLabelFont(42);
	Unlike->Draw("");
//	MixedBg->SetMarkerStyle(7);
//	MixedBg->SetMarkerColor(2);
	MixedBg->SetLineColor(2);
	MixedBg->SetLineWidth(2);
//	MixedBg->Draw("Csame");
	LikeBg->SetMarkerStyle(7);
	LikeBg->SetMarkerColor(4);
	LikeBg->SetLineColor(4);
	LikeBg->SetLineWidth(2);
	LikeBg->Draw("hist C same");
	RotateBg->SetMarkerStyle(7);
	RotateBg->SetMarkerColor(8);
	RotateBg->SetLineColor(8);
	RotateBg->SetLineWidth(2);
	RotateBg->Draw("hist C same");	
	leg->AddEntry(Unlike,"Unlike sign","pf");
//	leg->AddEntry(MixedBg,"Mixed","l");
	leg->AddEntry(LikeBg,"Like sign","l");
	leg->AddEntry(RotateBg,"Rotated","l");
	
	
	
	// **************************************************************************************************************
	// ************************************************ Mixed Event *************************************************
	// **************************************************************************************************************
	x1 = 1.7;
	x2 = 2.1;	
	nBinsInvMass = MixedBg->GetNbinsX();
	lowEdge = MixedBg->GetBinLowEdge(1);
	highEdge = MixedBg->GetBinLowEdge(nBinsInvMass+1);

	TCanvas *c2 = new TCanvas("c2");
	c2->cd();
	c2->SetRightMargin(0.05);

	TH1D *D0_mixed = new TH1D("D0_mixed","N_{K^{-}#pi^{+}}+N_{K^{+}#pi^{-}} - (N_{K^{-}#pi^{+}}+N_{K^{+}#pi^{-}})_{buffer event}",nBinsInvMass,lowEdge,highEdge);
	D0_mixed->Sumw2();
	for(int i=0; i<D0_mixed->GetNbinsX(); i++){
		D0_mixed->SetBinContent(i+1,Unlike->GetBinContent(1+i) - MixedBg->GetBinContent(1+i));
		D0_mixed->SetBinError(i+1,TMath::Sqrt(Unlike->GetBinError(1+i)*Unlike->GetBinError(1+i)+
											  MixedBg->GetBinError(1+i)*MixedBg->GetBinError(1+i)));
	}	
	D0_mixed->Draw("bar");
	D0_mixed->GetXaxis()->SetTitle("m_{K#pi} [GeV/c^{2}]");
	D0_mixed->GetYaxis()->SetTitle("counts");

	//D0_mixed->SetMarkerStyle(7);
	D0_mixed->SetMarkerColor(2);
	D0_mixed->SetFillColor(2);
	
	TF1 *fOutskirt = new TF1("fOutskirt",bgfitfunction,x1,x2,4);
	fOutskirt->FixParameter(3,0.);
	fOutskirt->FixParameter(2,0.);
	D0_mixed->Fit("fOutskirt","N","",x1,x2);
	cout<<"*****************************************************************************"<<endl;
	cout<<"Mixed residual bg chisquare/ndf="<<fOutskirt->GetChisquare()<<"/"<<fOutskirt->GetNDF()<<endl;
	cout<<"*****************************************************************************"<<endl;
	

	TPad *zoom1Mixed = new TPad("zoom1Mixed","zoom",0.5,0.25,0.94,0.89);
	zoom1Mixed->SetFillColor(10);
	zoom1Mixed->Draw();
	zoom1Mixed->cd();
	zoom1Mixed->SetRightMargin(0.05);
	D0_mixed->Clone("mixedClone");
	TH1D *mixClon = (TH1D*)gDirectory->Get("mixedClone");
	mixClon->SetTitle("");
	mixClon->Draw();
	mixClon->GetXaxis()->SetTitle("m_{K#pi} [GeV/c^{2}]");
	mixClon->SetMarkerStyle(8);
	mixClon->GetXaxis()->SetRangeUser(x1,x2);
	PlotLine(1.8,1.8,mixClon->GetMinimum(),mixClon->GetMaximum(),1,3,2);
	PlotLine(1.92,1.92,mixClon->GetMinimum(),mixClon->GetMaximum(),1,3,2);
	
	int x1Bin = D0_mixed->FindBin(x1);
	int x2Bin = D0_mixed->FindBin(x2);
	int nBins = x2Bin - x1Bin;
	



	// ******************************* Mixed Residual bg Subtraction *******************************
	fOutskirt->GetParameters(par);
	TF1 *fSignalMixed = new TF1("fSignalMixed",fitFunction,x1,x2,7);
	par[4]=40; par[5]=1.864; par[6]=0.011;
	fSignalMixed->SetParameters(par);
	fSignalMixed->SetParLimits(5,1.84,1.89);
	fSignalMixed->SetParLimits(6,0.003,0.025);
	fSignalMixed->SetNpx(1000);
//	fSignalMixed->FixParameter(4,0);
//	fSignalMixed->FixParameter(5,1.865);
//	fSignalMixed->FixParameter(6,0.0111);
	fSignalMixed->FixParameter(0,par[0]);
	fSignalMixed->FixParameter(1,par[1]);
	fSignalMixed->FixParameter(2,par[2]);
	fSignalMixed->FixParameter(3,par[3]);
	mixClon->Fit("fSignalMixed","","N",xs1,xs2);
	fSignalMixed->SetLineColor(2);
	fSignalMixed->Draw("same");
	double signal = fSignalMixed->GetParameter(4)/mixClon->GetBinWidth(1);
	double signalerror = fSignalMixed->GetParError(4)/mixClon->GetBinWidth(1);
	sprintf(parstring,"RwYld = %.0f #pm %.0f",signal,signalerror);
	c2->cd();
	TPad *infoMixed = new TPad("infoMixed","info",0.3,0.3,0.5,0.8);
	infoMixed->SetFillColor(10);
	infoMixed->Draw();
	infoMixed->cd();	
	TLatex tl;
	tl.SetTextSize(0.1);
	tl.SetTextColor(2);
	tl.DrawLatex(0.05,0.85,parstring);
	sprintf(parstring,"#chi^{2}/ndf = %.3f/%i",fSignalMixed->GetChisquare(),fSignalMixed->GetNDF());	
	tl.DrawLatex(0.05,0.7,parstring);
	sprintf(parstring,"#mu = %.3f #pm %.3f",fSignalMixed->GetParameter(5),fSignalMixed->GetParError(5));	
	tl.DrawLatex(0.05,0.55,parstring);
	sprintf(parstring,"#sigma = %.4f #pm %.4f",fSignalMixed->GetParameter(6),fSignalMixed->GetParError(6));	
	tl.DrawLatex(0.05,0.4,parstring);
	tl.SetTextColor(1);
	sprintf(parstring,"p_{T} = %.1f:%.1f [GeV/c]",pT1,pT2);
	tl.DrawLatex(0.05,0.25,parstring);
	sprintf(parstring,"dm_{K#pi} = %.3f [GeV/c^{2}]",mixClon->GetBinWidth(1));
	tl.DrawLatex(0.05,0.1,parstring);
	// **************************************************************************************************************
	// ********************************************* Rotated momentum ***********************************************
	// **************************************************************************************************************	
	x1 = 1.72;
	x2 = 2.1;
//    xs1 = 1.8;
//    xs2 = 1.92;
	TCanvas *c3 = new TCanvas("c3");
	c3->cd();
	c3->SetRightMargin(0.05);
		
	TH1D *D0_Rotate = new TH1D("D0_rotate","N_{K^{-}#pi^{+}}+N_{K^{+}#pi^{-}} - (N_{K^{-}(rot)#pi^{+}}+N_{K^{+}(rot)#pi^{-}})",nBinsInvMass,lowEdge,highEdge);
	D0_Rotate->Sumw2();
	for(int i=0; i<D0_Rotate->GetNbinsX(); i++){
		D0_Rotate->SetBinContent(i+1,Unlike->GetBinContent(1+i) - RotateBg->GetBinContent(1+i));
		D0_Rotate->SetBinError(i+1,TMath::Sqrt(Unlike->GetBinError(1+i)*Unlike->GetBinError(1+i)+
											  RotateBg->GetBinError(1+i)*RotateBg->GetBinError(1+i)));
	}	
	D0_Rotate->Draw("bar");
	D0_Rotate->GetXaxis()->SetTitle("m_{K#pi} [GeV/c^{2}]");
	D0_Rotate->GetYaxis()->SetTitle("counts");
	//D0_rotate->SetMarkerStyle(7);
	D0_Rotate->SetMarkerColor(kGreen+2);
	D0_Rotate->SetFillColor(kGreen+2);
	D0_Rotate->Fit("fOutskirt","N","",x1,x2);
	cout<<"*****************************************************************************"<<endl;
	cout<<"Rotate residual bg chisquare/ndf="<<fOutskirt->GetChisquare()<<"/"<<fOutskirt->GetNDF()<<endl;
	cout<<"*****************************************************************************"<<endl;
	
	TPad *zoom1Rotate = new TPad("zoom1Rotate","zoom",0.5,0.25,0.94,0.89);
	zoom1Rotate->SetFillColor(10);
	zoom1Rotate->Draw();
	zoom1Rotate->cd();
	zoom1Rotate->SetRightMargin(0.05);
	D0_Rotate->Clone("RotateClone");
	D0_Rotate->Clone("CloneR");
	TH1D *rotClon = (TH1D*)gDirectory->Get("RotateClone");
	int lowborder = rotClon->FindBin(1.83);
	int highborder = rotClon->FindBin(1.9);
//	for (int i=0; i<lowborder; i=i+2) {
//		float content1 = rotClon->GetBinContent(i+1);
//		float content2 = rotClon->GetBinContent(i+2);
//		float error1 = rotClon->GetBinError(i+1);
//		float error2 = rotClon->GetBinError(i+2);
//		rotClon->SetBinContent(i+1,(float)(content1+content2)/2);
//		rotClon->SetBinContent(i+2,0);
//		rotClon->SetBinError(i+1,TMath::Sqrt(1/(1/(error1*error1)+1/(error2*error2))));
//		rotClon->SetBinError(i+2,0);
//	}
//	for (int i=highborder; i<rotClon->GetNbinsX(); i=i+2) {
//		float content1 = rotClon->GetBinContent(i+1);
//		float content2 = rotClon->GetBinContent(i+2);
//		float error1 = rotClon->GetBinError(i+1);
//		float error2 = rotClon->GetBinError(i+2);
//		rotClon->SetBinContent(i+1,(float)(content1+content2)/2);
//		rotClon->SetBinContent(i+2,0);
//		rotClon->SetBinError(i+1,TMath::Sqrt(1/(1/(error1*error1)+1/(error2*error2))));
//		rotClon->SetBinError(i+2,0);
//	}
	
	rotClon->SetTitle("");
	rotClon->Draw();
	rotClon->GetXaxis()->SetTitle("m_{K#pi} [GeV/c^{2}]");
	rotClon->SetMarkerStyle(8);
	rotClon->GetXaxis()->SetRangeUser(x1,x2);
	PlotLine(1.8,1.8,rotClon->GetMinimum(),rotClon->GetMaximum(),1,3,2);
	PlotLine(1.92,1.92,rotClon->GetMinimum(),rotClon->GetMaximum(),1,3,2);
	
	x1Bin = D0_Rotate->FindBin(x1);
	x2Bin = D0_Rotate->FindBin(x2);
	nBins = x2Bin - x1Bin;
	
	// ******************************* Rotate Residual bg Subtraction *******************************
	fOutskirt->GetParameters(par);
	TF1 *fSignalRotate = new TF1("fSignalRotate",fitFunction,x1,x2,7);
	par[4]=40; par[5]=1.864; par[6]=0.011;
	fSignalRotate->SetParameters(par);
	fSignalRotate->SetParLimits(5,1.84,1.89);
	fSignalRotate->SetParLimits(6,0.003,0.04);
	fSignalRotate->SetNpx(1000);
//	fSignalRotate->FixParameter(4,0); // zero signal hypothesis
//	fSignalRotate->FixParameter(5,1.864);
//	fSignalRotate->FixParameter(6,0.10);
//	fSignalRotate->FixParameter(0,par[0]);
//	fSignalRotate->FixParameter(1,par[1]);
	fSignalRotate->FixParameter(2,par[2]);
	fSignalRotate->FixParameter(3,par[3]);
	rotClon->Fit("fSignalRotate","","N",xs1,xs2);
	fSignalRotate->SetLineColor(kGreen+2);
	fSignalRotate->Draw("same");
	signal = fSignalRotate->GetParameter(4)/rotClon->GetBinWidth(1);
	signalerror = fSignalRotate->GetParError(4)/rotClon->GetBinWidth(1);
	sprintf(parstring,"RwYld = %.0f #pm %.0f",signal,signalerror);
	c3->cd();
	TPad *infoRotate = new TPad("infoRotate","info",0.3,0.3,0.5,0.8);
	infoRotate->SetFillColor(10);
	infoRotate->Draw();
	infoRotate->cd();
	tl.SetTextSize(0.1);
	tl.SetTextColor(kGreen+2);
    tl.SetTextFont(42);
	tl.DrawLatex(0.05,0.85,parstring);
	sprintf(parstring,"#chi^{2}/ndf = %.3f/%i",fSignalRotate->GetChisquare(),fSignalRotate->GetNDF());	
	tl.DrawLatex(0.05,0.7,parstring);
	sprintf(parstring,"#mu = %.3f #pm %.3f",fSignalRotate->GetParameter(5),fSignalRotate->GetParError(5));	
	tl.DrawLatex(0.05,0.55,parstring);
	sprintf(parstring,"#sigma = %.4f #pm %.4f",fSignalRotate->GetParameter(6),fSignalRotate->GetParError(6));	
	tl.DrawLatex(0.05,0.4,parstring);
	tl.SetTextColor(1);
	sprintf(parstring,"p_{T} = %.1f:%.1f [GeV/c]",pT1,pT2);
	tl.DrawLatex(0.05,0.25,parstring);
	sprintf(parstring,"dm_{K#pi} = %.3f [GeV/c^{2}]",rotClon->GetBinWidth(1));
	tl.DrawLatex(0.05,0.1,parstring);
	
	
	
	// **************************************************************************************************************
	// ********************************************* Like sign ***********************************************
	// **************************************************************************************************************	
	x1 = 1.72;
	x2 = 2.1;
//    xs1=1.8;
//    xs2=1.92;
	TCanvas *c4 = new TCanvas("c4");
	c4->cd();
	c4->SetRightMargin(0.05);
  TH1D *D0_Like;
	if(GeometricMean) 	D0_Like = new TH1D("D0_Like","N_{K^{-}#pi^{+}}+N_{K^{+}#pi^{-}} - 2 #sqrt{N_{K^{-}#pi^{-}}N_{K^{+}#pi^{+}}}",nBinsInvMass,lowEdge,highEdge);
	else D0_Like = new TH1D("D0_Like","N_{K^{-}#pi^{+}}+N_{K^{+}#pi^{-}} - (N_{K^{-}#pi^{-}}+N_{K^{+}#pi^{+}})",nBinsInvMass,lowEdge,highEdge);

	D0_Like->Sumw2();
	for(int i=0; i<D0_Like->GetNbinsX(); i++){
		D0_Like->SetBinContent(i+1,Unlike->GetBinContent(1+i) - LikeBg->GetBinContent(1+i));
		D0_Like->SetBinError(i+1,TMath::Sqrt(Unlike->GetBinError(1+i)*Unlike->GetBinError(1+i)+
											   LikeBg->GetBinError(1+i)*LikeBg->GetBinError(1+i)));
	}	
	D0_Like->Draw("bar");
	D0_Like->GetXaxis()->SetTitle("m_{K#pi} [GeV/c^{2}]");
	D0_Like->GetYaxis()->SetTitle("counts");
	D0_Like->SetFillColor(4);
//	D0_Like->SetMarkerStyle(7);
	D0_Like->SetMarkerColor(4);
//	D0_Like->Fit("fLinearOutskirt","N","",x1,x2);
	D0_Like->Fit("fOutskirt","N","",x1,x2);
	cout<<"*****************************************************************************"<<endl;
	cout<<"LikeSign residual bg chisquare/ndf="<<fOutskirt->GetChisquare()<<"/"<<fOutskirt->GetNDF()<<endl;
	cout<<"*****************************************************************************"<<endl;

	
	TPad *zoom1Like = new TPad("zoom1Like","zoom",0.5,0.25,0.94,0.89);
	zoom1Like->SetFillColor(10);
	zoom1Like->Draw();
	zoom1Like->cd();
	zoom1Like->SetRightMargin(0.05);
	D0_Like->Clone("LikeClone");
	D0_Like->Clone("CloneL");

	TH1D *likeClon = (TH1D*)gDirectory->Get("LikeClone");
//	for (int i=0; i<lowborder; i=i+2) {
//		float content1 = likeClon->GetBinContent(i+1);
//		float content2 = likeClon->GetBinContent(i+2);
//		float error1 = likeClon->GetBinError(i+1);
//		float error2 = likeClon->GetBinError(i+2);
//		likeClon->SetBinContent(i+1,(float)(content1+content2)/2);
//		likeClon->SetBinContent(i+2,0);
//		likeClon->SetBinError(i+1,TMath::Sqrt(1/(1/(error1*error1)+1/(error2*error2))));
//		likeClon->SetBinError(i+2,0);
//	}
//	for (int i=highborder; i<likeClon->GetNbinsX(); i=i+2) {
//		float content1 = likeClon->GetBinContent(i+1);
//		float content2 = likeClon->GetBinContent(i+2);
//		float error1 = likeClon->GetBinError(i+1);
//		float error2 = likeClon->GetBinError(i+2);
//		likeClon->SetBinContent(i+1,(float)(content1+content2)/2);
//		likeClon->SetBinContent(i+2,0);
//		likeClon->SetBinError(i+1,TMath::Sqrt(1/(1/(error1*error1)+1/(error2*error2))));
//		likeClon->SetBinError(i+2,0);
//	}
	likeClon->SetTitle("");
	likeClon->Draw();
	likeClon->GetXaxis()->SetTitle("m_{K#pi} [GeV/c^{2}]");
	likeClon->SetMarkerStyle(8);
	likeClon->GetXaxis()->SetRangeUser(x1,x2);
	PlotLine(1.8,1.8,likeClon->GetMinimum(),likeClon->GetMaximum(),1,3,2);
	PlotLine(1.92,1.92,likeClon->GetMinimum(),likeClon->GetMaximum(),1,3,2);
	
	x1Bin = D0_Like->FindBin(x1);
	x2Bin = D0_Like->FindBin(x2);
	nBins = x2Bin - x1Bin;
	
	// ******************************* Like Residual bg Subtraction *******************************
	fOutskirt->GetParameters(par);
	TF1 *fSignalLike = new TF1("fSignalLike",fitFunction,x1,x2,7);
	par[4]=40; par[5]=1.864; par[6]=0.011;
	fSignalLike->SetParameters(par);
	fSignalLike->SetParLimits(4,0,1000);
	fSignalLike->SetParLimits(5,1.84,1.89);
	fSignalLike->SetParLimits(6,0.005,0.025);
	fSignalLike->SetNpx(1000);
//	fSignalLike->FixParameter(4,0);	
//	fSignalLike->FixParameter(5,1.864);
//	fSignalLike->FixParameter(6,0.011);
//	fSignalLike->FixParameter(0,par[0]);
//	fSignalLike->FixParameter(1,par[1]);
	fSignalLike->FixParameter(2,par[2]);
	fSignalLike->FixParameter(3,par[3]);
	likeClon->Fit("fSignalLike","","N",xs1,xs2);
	fSignalLike->SetLineColor(4);
	fSignalLike->Draw("same");
	signal = fSignalLike->GetParameter(4)/likeClon->GetBinWidth(1);
	signalerror = fSignalLike->GetParError(4)/likeClon->GetBinWidth(1);
	sprintf(parstring,"RwYld = %.0f #pm %.0f",signal,signalerror);
	c4->cd();
	TPad *infoLike = new TPad("infoLike","info",0.3,0.3,0.5,0.8);
	infoLike->SetFillColor(10);
	infoLike->Draw();
	infoLike->cd();	
	tl.SetTextSize(0.1);
    tl.SetTextFont(42);
	tl.SetTextColor(4);
	tl.DrawLatex(0.05,0.85,parstring);
	sprintf(parstring,"#chi^{2}/ndf = %.3f/%i",fSignalLike->GetChisquare(),fSignalLike->GetNDF());	
	tl.DrawLatex(0.05,0.7,parstring);
	sprintf(parstring,"#mu = %.3f #pm %.3f",fSignalLike->GetParameter(5),fSignalLike->GetParError(5));	
	tl.DrawLatex(0.05,0.55,parstring);
	sprintf(parstring,"#sigma = %.4f #pm %.4f",fSignalLike->GetParameter(6),fSignalLike->GetParError(6));	
	tl.DrawLatex(0.05,0.4,parstring);
	tl.SetTextColor(1);
	sprintf(parstring,"p_{T} = %.1f:%.1f [GeV/c]",pT1,pT2);
	tl.DrawLatex(0.05,0.25,parstring);
	sprintf(parstring,"dm_{K#pi} = %.3f [GeV/c^{2}]",likeClon->GetBinWidth(1));
	tl.DrawLatex(0.05,0.1,parstring);

	TCanvas *c20 = new TCanvas("c20");
	c20->cd();
	mixClon->DrawCopy();
	likeClon->DrawCopy("same");
	rotClon->DrawCopy("same");
	
	TLegend *resultsLegend = new TLegend(0.2,0.82,0.5,0.98,"","brNDC");
	resultsLegend->AddEntry(mixClon,"mixed event","p");
	resultsLegend->AddEntry(rotClon,"rotated momentum","p");
	resultsLegend->AddEntry(likeClon,"like sign","p");
	resultsLegend->Draw();
	infoMixed->DrawClone();
	infoRotate->DrawClone();
	infoLike->DrawClone();
	
//	TCanvas *c30 = new TCanvas("c30");
//	c30->cd();
//	
//	D0_Like->SetMarkerStyle(20);
//	D0_rotate->SetMarkerStyle(20);
//	D0_mixed->SetMarkerStyle(20);
//	
//	D0_Like->Draw();
//	D0_rotate->Draw("same");
//	D0_mixed->Draw("same");
//	
//	TLegend *resultsLegend = new TLegend(0.2,0.82,0.5,0.98,"","brNDC");
//	resultsLegend->AddEntry(D0_Like,"mixed event","p");
//	resultsLegend->AddEntry(D0_rotate,"rotated momentum","p");
//	resultsLegend->AddEntry(D0_mixed,"like sign","p");
//	resultsLegend->Draw();
	
	
	//********************************** Draw and Save Zooms ************************************
	TCanvas *c10 = new TCanvas("c10","c10",400,400);
	c10->SetBottomMargin(0.15);
	c10->SetLeftMargin(0.15);
	c10->SetTickx();
	c10->SetTicky();
//	hSignal_Mixed->SetTitle("Mixed event");
//	hSignal_Mixed->GetXaxis()->SetRangeUser(1.73,2.05);
//	hSignal_Mixed->GetYaxis()->SetRangeUser(-800,2500);
//	hSignal_Mixed->GetYaxis()->SetTitle("Counts");
//	hSignal_Mixed->GetYaxis()->SetTitleOffset(1.7);
//	hSignal_Mixed->GetYaxis()->SetLabelSize(0.04);
////	hSignal_Mixed->SetMarkerStyle(8);
//	hSignal_Mixed->SetMarkerSize(0.7);
//	hSignal_Mixed->SetMarkerColor(2);
//	hSignal_Mixed->Draw();
//	infoMixed->SetPad(0.6,0.5,1.,1.);
//	infoMixed->SetFillColor(10);
//	infoMixed->Draw("same");
//	c10->cd(2);
//	gPad->SetBottomMargin(0.2);
//	gPad->SetLeftMargin(0.15);
	rotClon->SetTitle("After Rotated bg. subtracted");
	rotClon->GetXaxis()->SetRangeUser(1.73,2.05);
	rotClon->GetYaxis()->SetRangeUser(-1200,3500);
	rotClon->GetYaxis()->SetTitle("Raw Yield (/0.01 GeV/c^{2})");
	rotClon->GetYaxis()->SetTitleOffset(1.4);
	rotClon->GetYaxis()->SetLabelSize(0.05);
	rotClon->GetXaxis()->SetLabelSize(0.05);
	rotClon->GetYaxis()->SetTitleSize(0.05);
	rotClon->GetXaxis()->SetTitleSize(0.05);
	rotClon->GetYaxis()->SetLabelFont(42);
	rotClon->GetXaxis()->SetLabelFont(42);
	rotClon->GetYaxis()->SetTitleFont(42);
	rotClon->GetXaxis()->SetTitleFont(42);
	rotClon->SetMarkerSize(2);
	rotClon->SetLineWidth(2);
	//	rotClon->SetMarkerStyle(8);
	rotClon->SetMarkerColor(kGreen+1);
	rotClon->Draw();
	infoRotate->SetPad(0.6,0.5,1.,0.9);
	infoRotate->SetFillColor(10);
	infoRotate->Draw("same");
	TCanvas *c11 = new TCanvas("c11","c11",400,400);
	c11->SetBottomMargin(0.15);
	c11->SetTickx();
	c11->SetTicky();
	c11->SetLeftMargin(0.15);
	likeClon->SetTitle("After Like-sign bg. subtracted");
	likeClon->GetXaxis()->SetRangeUser(1.73,2.05);
	likeClon->GetYaxis()->SetRangeUser(-1200,3500);
	likeClon->GetYaxis()->SetTitle("Raw Yield (/0.01 GeV/c^{2})");
	likeClon->GetYaxis()->SetTitleOffset(1.4);
	likeClon->GetYaxis()->SetLabelSize(0.05);
	likeClon->GetXaxis()->SetLabelSize(0.05);
	likeClon->GetYaxis()->SetTitleSize(0.05);
	likeClon->GetXaxis()->SetTitleSize(0.05);
    likeClon->GetYaxis()->SetLabelFont(42);
	likeClon->GetXaxis()->SetLabelFont(42);
	likeClon->GetYaxis()->SetTitleFont(42);
	likeClon->GetXaxis()->SetTitleFont(42);
	//	likeClon->SetMarkerStyle(8);
	likeClon->SetMarkerColor(4);
	likeClon->SetMarkerSize(2);
	likeClon->SetLineWidth(2);
	likeClon->Draw();
	infoLike->SetPad(0.6,0.5,1.,0.9);
	infoLike->SetFillColor(10);
	infoLike->Draw("same");
//	
//	TF1 *signalFcn = new TF1("signalFcn",lorentzianPeak,0.8,1.05,3);
//	c4->cd();
//	signalFcn->SetParameters(1000,1,0.892);
//	D0_Like->Fit("signalFcn","","",0.8,0.95);
	


	//********************************** Draw a more plots into one **************************************

	c1->cd();
	TH1D *CloneR = (TH1D*)gDirectory->Get("CloneR");
	TH1D *CloneL = (TH1D*)gDirectory->Get("CloneL");
	CloneR->Scale(2);
	CloneL->Scale(2);
	CloneR->SetMarkerColor(8);
	CloneR->SetLineColor(8);
	CloneR->SetMarkerStyle(8);
	CloneL->SetMarkerColor(4);
	CloneL->SetLineColor(4);
	CloneL->SetMarkerStyle(8);

	CloneR->Draw("same");
	CloneL->Draw("same");

	TGaxis* RawYieldAxis = new TGaxis(2.5,0,2.5,14e4,0,7e4,510,"+L");
	//+ : draw on positive side
	//L : left adjusted
	RawYieldAxis->SetName("RawYieldAxis");
	RawYieldAxis->SetLineColor(8);
	RawYieldAxis->SetTextColor(4);
	// RawYieldAxis->SetTitle("Raw Yield (/0.1 GeV/c^{2})");
	RawYieldAxis->SetLabelColor(4);
	RawYieldAxis->SetLabelFont(42);
	RawYieldAxis->Draw();

	leg->AddEntry(CloneL,"2*(Unlike - Like)","p");
	leg->AddEntry(CloneR,"2*(Unlike - Rotated)","p");

	leg->Draw("same");



	
}

void PlotLine(Double_t x1_val, Double_t x2_val, Double_t y1_val,
			  Double_t y2_val, Int_t Line_Col, Int_t LineWidth, Int_t LineStyle)
{
	TLine* Zero_line = new TLine();
	Zero_line -> SetX1(x1_val);
	Zero_line -> SetX2(x2_val);
	Zero_line -> SetY1(y1_val);
	Zero_line -> SetY2(y2_val);
	Zero_line -> SetLineWidth(LineWidth);
	Zero_line -> SetLineStyle(LineStyle);
	Zero_line -> SetLineColor(Line_Col);
	Zero_line -> Draw();
	//delete Zero_line;
}

//PlotLine(0.896,0.896,0,hist_max_x_val,1,2,2); // x1_val, x2_val, y1_val, y2_val, Line_Col, LineWidth, LineStyle

