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

void RawYieldPt_sg_nomixedE(float pT1, float pT2, int rebin) {
	
	float x1,x2,lowEdge,highEdge,normEdge1,normEdge2, normEdge3, normEdge4, ptBin1,ptBin2;
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
	TH1D *MixedBg = new TH1D("MixedBg","",nBinsInvMass,lowEdge,highEdge);
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
		for(int i=0; i<=nBinsInvMass+1; i++){
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
				LikeBg->SetBinError(i,TMath::Sqrt(a*b)*TMath::Sqrt(((da/a)*(da/a)) + ((db/b)*(db/b))));
			}
		}
	}else LikeBg->Add(LikeBgD_py,LikeBgDbar_py,1,1);	
	RotateBg->Add(DRotate_py,DbarRotate_py,1,1);
	
	Double_t Scalefactorright;
	Double_t Scalefactorleft;
	Double_t Scalefactor;

	normEdge1 = Unlike->FindBin(1.7);
	normEdge2 = Unlike->FindBin(1.8);

	normEdge3 = Unlike->FindBin(1.9);
	normEdge4 = Unlike->FindBin(2.0);

	Scalefactorleft = Unlike->Integral(normEdge1,normEdge2)/MixedBg->Integral(normEdge1,normEdge2);
	Scalefactorright = Unlike->Integral(normEdge3,normEdge4)/MixedBg->Integral(normEdge3,normEdge4); 

	Scalefactor = (Scalefactorleft + Scalefactorright)/2.0;

	/////MixedBg->Scale(Scalefactor);
	cout<<"Mixed event Scale factor: "<<Scalefactor<<endl;
	Scalefactorleft = Unlike->Integral(normEdge1,normEdge2)/RotateBg->Integral(normEdge1,normEdge2);
	Scalefactorright = Unlike->Integral(normEdge3,normEdge4)/RotateBg->Integral(normEdge3,normEdge4);
	Scalefactor = (Scalefactorleft + Scalefactorright)/2.0;

	///////Scalefactor = 0.9999;
 	RotateBg->Scale(Scalefactorleft);
  cout<<"Rotated momentum scale factor: "<<Scalefactor<<endl;
	Scalefactorleft = Unlike->Integral(normEdge1,normEdge2)/LikeBg->Integral(normEdge1,normEdge2);
	Scalefactorright = Unlike->Integral(normEdge3,normEdge4)/LikeBg->Integral(normEdge3,normEdge4);
	Scalefactor = (Scalefactorleft + Scalefactorright)/2.0;
	LikeBg->Scale(Scalefactorleft);
  cout<<"Like sign scale factor: "<<Scalefactor<<endl;

	//********************************************* Draw Invariant mass **********************************************
	TCanvas *c1 = new TCanvas("c1");
	c1->cd();
	c1->SetRightMargin(0.05);
	Unlike->GetXaxis()->SetTitle("M_{K#pi} [GeV/c^{2}]");
	Unlike->GetYaxis()->SetTitle("Entries (/0.01 GeV/c^{2})");
	Unlike->SetMarkerStyle(4);
	Unlike->GetYaxis()->SetTitleFont(42);
	Unlike->GetYaxis()->SetLabelFont(42);
	//////Unlike->GetYaxis()->SetRangeUser(-0.04e6, 0.5e6);
	Unlike->Draw("");
	MixedBg->SetMarkerStyle(7);
	MixedBg->SetMarkerColor(kCyan+3);
	MixedBg->SetLineColor(kCyan+3);
	MixedBg->SetLineWidth(2);
	//////MixedBg->Draw("hist C same");
    
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
		
	leg->AddEntry(Unlike,"Unlike Sign","pf");
	leg->AddEntry(LikeBg,"Like Sign","l");
	leg->AddEntry(RotateBg,"Rotated Momentum","l");
	//////leg->AddEntry(MixedBg,"Mixed Event","l");
		
	// **************************************************************************************************************
	// ********************************************* Like sign ***********************************************
	// **************************************************************************************************************	
	x1 = 1.7;
	x2 = 2.1;
//    x1=1.8;
//    x2=1.92;
	TCanvas *c4 = new TCanvas("c4");
	c4->cd();
	c4->SetRightMargin(0.05);
  TH1D *D0_Like;
	if(GeometricMean) 	D0_Like = new TH1D("D0_Like","N_{K^{-}#pi^{+}}+N_{K^{+}#pi^{-}} - 2 #sqrt{N_{K^{-}#pi^{-}}N_{K^{+}#pi^{+}}}",nBinsInvMass,lowEdge,highEdge);
	else D0_Like = new TH1D("D0_Like","N_{K^{-}#pi^{+}}+N_{K^{+}#pi^{-}} - (N_{K^{-}#pi^{-}}+N_{K^{+}#pi^{+}})",nBinsInvMass,lowEdge,highEdge);

	D0_Like->Sumw2();
	for(int i=0; i<= D0_Like->GetNbinsX()+1; i++){
		D0_Like->SetBinContent(i,Unlike->GetBinContent(i) - LikeBg->GetBinContent(i));
		D0_Like->SetBinError(i,TMath::Sqrt(((Unlike->GetBinError(i)*Unlike->GetBinError(i)) +
											   (LikeBg->GetBinError(i)*LikeBg->GetBinError(i)))));
	}	
	D0_Like->Draw("bar");
	D0_Like->GetXaxis()->SetTitle("m_{K#pi} [GeV/c^{2}]");
	D0_Like->GetYaxis()->SetTitle("counts");
	D0_Like->SetFillColor(4);
//	D0_Like->SetMarkerStyle(7);
    D0_Like->SetLineColor(4);
	D0_Like->SetMarkerColor(4);

	TF1 *fOutskirt = new TF1("fOutskirt",bgfitfunction,x1,x2,4);
	////fOutskirt->SetParameter(2,0);
	////fOutskirt->SetParameter(1,4.33488e+04);
	////fOutskirt->SetParameter(0,-4.22559e+04);
    fOutskirt->FixParameter(2,0.);
	fOutskirt->FixParameter(3,0.);
//	D0_Like->Fit("fLinearOutskirt","","",x1,x2);
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
	
	int x1Bin = D0_Like->FindBin(x1);
	int x2Bin = D0_Like->FindBin(x2);
	int nBins = x2Bin - x1Bin;
	
	// ******************************* Like Residual bg Subtraction *******************************
	fOutskirt->GetParameters(par);
	TF1 *fSignalLike = new TF1("fSignalLike",fitFunction,x1,x2,7);
	fSignalLike->SetLineColor(kRed);
	fSignalLike->SetLineWidth(4);
	par[4]=40; par[5]=1.864; par[6]=0.011;
	fSignalLike->SetParameters(par);
	//////fSignalLike->SetParLimits(4,0,1000);
	fSignalLike->SetParLimits(5,1.850,1.890);
	fSignalLike->SetParLimits(6,0.006,0.018);
	////fSignalLike->SetNpx(1000);
//	fSignalLike->FixParameter(4,0);	
//	fSignalLike->FixParameter(5,1.864);
//	fSignalLike->FixParameter(6,0.011);
//	fSignalLike->FixParameter(0,par[0]);
//	fSignalLike->FixParameter(1,par[1]);
    fSignalLike->FixParameter(2,par[2]);
	fSignalLike->FixParameter(3,par[3]);
	likeClon->Fit("fSignalLike","","",x1,x2);
	fSignalLike->Draw("same");


	TF1 *fsignalonlyLike = new TF1("fsignalonlylike",normalDistribution,x1,x2,3);
	fsignalonlyLike->FixParameter(0,fSignalLike->GetParameter(4));
	fsignalonlyLike->FixParameter(1,fSignalLike->GetParameter(5));
	fsignalonlyLike->FixParameter(2,fSignalLike->GetParameter(6));

	fsignalonlyLike->SetLineColor(kBlack);
	fsignalonlyLike->SetLineWidth(2);
	fsignalonlyLike->SetLineStyle(10);
	fsignalonlyLike->Draw("same");

	TF1 *fresbgonlyLike = new TF1("fresbgonlyLike",parabola,x1,x2,4);
	fresbgonlyLike->FixParameter(0,fSignalLike->GetParameter(0));
	fresbgonlyLike->FixParameter(1,fSignalLike->GetParameter(1));
	fresbgonlyLike->FixParameter(2,fSignalLike->GetParameter(2));
	fresbgonlyLike->FixParameter(3,fSignalLike->GetParameter(3));

	fresbgonlyLike->SetLineColor(kRed);
	fresbgonlyLike->SetLineWidth(4);
	fresbgonlyLike->SetLineStyle(6);
	fresbgonlyLike->Draw("same");
	
	float signal = fSignalLike->GetParameter(4)/likeClon->GetBinWidth(1);
	float signalerror = fSignalLike->GetParError(4)/likeClon->GetBinWidth(1);
	sprintf(parstring,"RwYld = %.0f #pm %.0f",signal,signalerror);
	c4->cd();
	TPad *infoLike = new TPad("infoLike","info",0.3,0.3,0.5,0.8);
	infoLike->SetFillColor(10);
	infoLike->Draw();
	infoLike->cd();	

	TLatex tl;
	tl.SetTextSize(0.12);
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

	float x1a = (fSignalLike->GetParameter(5) - (3*fSignalLike->GetParameter(6)));
	float x2a = (fSignalLike->GetParameter(5) + (3*fSignalLike->GetParameter(6)));

	cout << x1a << "     " << x2a << endl;

	int x1aBin = Unlike->FindBin(x1a);
	int x2aBin = Unlike->FindBin(x2a);

	cout << x1aBin << "     " << x2aBin << endl;

	float x1aBinLowEdge = Unlike->GetBinLowEdge(x1aBin);
	float x2aBinHighEdge = Unlike->GetBinLowEdge(x2aBin) + Unlike->GetBinWidth(x2aBin);

	cout << x1aBinLowEdge << "    "  << x2aBinHighEdge << endl;

	float signalgauss = fsignalonlyLike->Integral(x1aBinLowEdge,x2aBinHighEdge)/likeClon->GetBinWidth(1);
	cout << signalgauss << "     " << sqrt(Unlike->Integral(x1aBin, x2aBin)) << endl;

	float sg = signalgauss/sqrt(Unlike->Integral(x1aBin, x2aBin));
	///sg = signalgauss/sqrt(signalgauss + 1*(fresbgonlyLike->Integral(x1a,x2a)/likeClon->GetBinWidth(1)) + 1*(LikeBg->Integral(x1aBin, x2aBin)));
    
	cout<<"*****************************************************************************"<<endl;
	cout << "significance from like-sign  =    " <<  sg << endl;
	cout<<"*****************************************************************************"<<endl;

	////sprintf(parstring,"dm_{K#pi} = %.3f [GeV/c^{2}]",likeClon->GetBinWidth(1));
	////tl.DrawLatex(0.05,0.1,parstring);


	// **************************************************************************************************************
	// ********************************************* Rotated momentum ***********************************************
	// **************************************************************************************************************	
	x1 = 1.70;
	x2 = 2.1;
//    x1 = 1.8;
//    x2 = 1.92;
	TCanvas *c3 = new TCanvas("c3");
	c3->cd();
	c3->SetRightMargin(0.05);
		
	TH1D *D0_Rotate = new TH1D("D0_rotate","N_{K^{-}#pi^{+}}+N_{K^{+}#pi^{-}} - (N_{K^{-}(rot)#pi^{+}}+N_{K^{+}(rot)#pi^{-}})",nBinsInvMass,lowEdge,highEdge);
	D0_Rotate->Sumw2();
	for(int i=0; i<= D0_Rotate->GetNbinsX()+1; i++){
		D0_Rotate->SetBinContent(i,Unlike->GetBinContent(i) - RotateBg->GetBinContent(i));
		D0_Rotate->SetBinError(i,TMath::Sqrt(Unlike->GetBinError(i)*Unlike->GetBinError(i) +
											  RotateBg->GetBinError(i)*RotateBg->GetBinError(i)));
	}	
	D0_Rotate->Draw("bar");
	D0_Rotate->GetXaxis()->SetTitle("m_{K#pi} [GeV/c^{2}]");
	D0_Rotate->GetYaxis()->SetTitle("counts");
	//D0_rotate->SetMarkerStyle(7);
	D0_Rotate->SetMarkerColor(kGreen+2);
	D0_Rotate->SetLineColor(kGreen+2);
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
	fSignalRotate->SetLineColor(kRed);
	fSignalRotate->SetLineWidth(4);
	fSignalRotate->SetParameters(par);
	/////fSignalRotate->SetParLimits(4,0,1000);
	/////cout << fSignalRotate->GetParameter(6) << "     " << fSignalLike->GetParameter(6) << endl;
	fSignalRotate->SetParLimits(5,1.850,fSignalLike->GetParameter(5));
	fSignalRotate->SetParLimits(6,0.007,fSignalLike->GetParameter(6));
	////fSignalRotate->SetNpx(1000);
//	fSignalRotate->FixParameter(4,0); // zero signal hypothesis
//	fSignalRotate->FixParameter(5,1.864);
//	fSignalRotate->FixParameter(6,0.10);
//	fSignalRotate->FixParameter(0,par[0]);
    ////fSignalRotate->FixParameter(4,par[4]);
	fSignalRotate->FixParameter(2,par[2]);
	fSignalRotate->FixParameter(3,par[3]);
	rotClon->Fit("fSignalRotate","","",x1,x2);
	fSignalRotate->Draw("same");

	TF1 *fsignalonlyrot = new TF1("fsignalonlyrot",normalDistribution,x1,x2,3);
	fsignalonlyrot->FixParameter(0,fSignalRotate->GetParameter(4));
	fsignalonlyrot->FixParameter(1,fSignalRotate->GetParameter(5));
	fsignalonlyrot->FixParameter(2,fSignalRotate->GetParameter(6));

	fsignalonlyrot->SetLineColor(kBlack);
	fsignalonlyrot->SetLineStyle(10);
	fsignalonlyrot->SetLineWidth(2);
	fsignalonlyrot->Draw("same");

	TF1 *fresbgonlyrot = new TF1("fresbgonlyrot",parabola,x1,x2,4);
	fresbgonlyrot->FixParameter(0,fSignalRotate->GetParameter(0));
	fresbgonlyrot->FixParameter(1,fSignalRotate->GetParameter(1));
	fresbgonlyrot->FixParameter(2,fSignalRotate->GetParameter(2));
	fresbgonlyrot->FixParameter(3,fSignalRotate->GetParameter(3));

	fresbgonlyrot->SetLineColor(kRed);
	fresbgonlyrot->SetLineStyle(6);
	fresbgonlyrot->SetLineWidth(4);
	fresbgonlyrot->Draw("same");
	

	signal = fSignalRotate->GetParameter(4)/rotClon->GetBinWidth(1);
	signalerror = fSignalRotate->GetParError(4)/rotClon->GetBinWidth(1);
	sprintf(parstring,"RwYld = %.0f #pm %.0f",signal,signalerror);
	c3->cd();
	TPad *infoRotate = new TPad("infoRotate","info",0.3,0.3,0.5,0.8);
	infoRotate->SetFillColor(10);
	infoRotate->Draw();
	infoRotate->cd();
	tl.SetTextSize(0.12);
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

	x1a = (fSignalRotate->GetParameter(5) - (3*fSignalRotate->GetParameter(6)));
	x2a = (fSignalRotate->GetParameter(5) + (3*fSignalRotate->GetParameter(6)));

	cout << x1a << "     " << x2a << endl;

	x1aBin = Unlike->FindBin(x1a);
	x2aBin = Unlike->FindBin(x2a);

	cout << x1aBin << "     " << x2aBin << endl;

	x1aBinLowEdge = Unlike->GetBinLowEdge(x1aBin);
	x2aBinHighEdge = Unlike->GetBinLowEdge(x2aBin) + Unlike->GetBinWidth(x2aBin);

	cout << x1aBinLowEdge << "    "  << x2aBinHighEdge << endl;

    signalgauss = fsignalonlyrot->Integral(x1aBinLowEdge,x2aBinHighEdge)/rotClon->GetBinWidth(1);

	cout << signalgauss << "     " << sqrt(Unlike->Integral(x1aBin, x2aBin)) << endl;

	sg = signalgauss/sqrt(Unlike->Integral(x1aBin, x2aBin));
	///sg = signalgauss/sqrt(signalgauss + 1*(fresbgonlyrot->Integral(x1a,x2a)/rotClon->GetBinWidth(1)) + 1*(RotateBg->Integral(x1aBin, x2aBin)));
    cout<<"*****************************************************************************"<<endl;
	cout << "significance from track-rotation  =    " <<  sg << endl;
	////sprintf(parstring,"sg = %.3f",sg);
	////tl.DrawLatex(0.05,0.1,parstring);

	cout<<"*****************************************************************************"<<endl;


	// **************************************************************************************************************
	// ************************************************ Mixed Event *************************************************
	// **************************************************************************************************************
	x1 = 1.70;
	x2 = 2.1;	
	nBinsInvMass = MixedBg->GetNbinsX();
	lowEdge = MixedBg->GetBinLowEdge(1);
	highEdge = MixedBg->GetBinLowEdge(nBinsInvMass+1);

	TCanvas *c2 = new TCanvas("c2");
	c2->cd();
	c2->SetRightMargin(0.05);

	TH1D *D0_mixed = new TH1D("D0_mixed","N_{K^{-}#pi^{+}}+N_{K^{+}#pi^{-}} - (N_{K^{-}#pi^{+}}+N_{K^{+}#pi^{-}})_{buffer event}",nBinsInvMass,lowEdge,highEdge);
	D0_mixed->Sumw2();
	for(int i=0; i<= D0_mixed->GetNbinsX()+1; i++){
		D0_mixed->SetBinContent(i,Unlike->GetBinContent(i) - MixedBg->GetBinContent(i));
		D0_mixed->SetBinError(i,TMath::Sqrt(((Unlike->GetBinError(i)*Unlike->GetBinError(i)) +
											  (MixedBg->GetBinError(i)*MixedBg->GetBinError(i)))));
	}	
	/////D0_mixed->Draw("bar");
	D0_mixed->GetXaxis()->SetTitle("m_{K#pi} [GeV/c^{2}]");
	D0_mixed->GetYaxis()->SetTitle("counts");

	//D0_mixed->SetMarkerStyle(7);
	D0_mixed->SetLineColor(kCyan+3);
	D0_mixed->SetMarkerColor(kCyan+3);
	D0_mixed->SetFillColor(kCyan+3);
	
	
	D0_mixed->Fit("fOutskirt","N","",x1,x2);
	cout<<"*****************************************************************************"<<endl;
	cout<<"Mixed residual bg chisquare/ndf="<<fOutskirt->GetChisquare()<<"/"<<fOutskirt->GetNDF()<<endl;
	cout<<"*****************************************************************************"<<endl;
	

	TPad *zoom1Mixed = new TPad("zoom1Mixed","zoom",0.5,0.25,0.94,0.89);
	zoom1Mixed->SetFillColor(10);
	////zoom1Mixed->Draw();
	zoom1Mixed->cd();
	zoom1Mixed->SetRightMargin(0.05);
	D0_mixed->Clone("mixedClone");
	D0_mixed->Clone("CloneM");


	TH1D *mixClon = (TH1D*)gDirectory->Get("mixedClone");
	mixClon->SetTitle("");
	mixClon->Draw();
	mixClon->GetXaxis()->SetTitle("m_{K#pi} [GeV/c^{2}]");
	mixClon->SetMarkerStyle(8);
	mixClon->GetXaxis()->SetRangeUser(x1,x2);
	PlotLine(1.8,1.8,mixClon->GetMinimum(),mixClon->GetMaximum(),1,3,2);
	PlotLine(1.92,1.92,mixClon->GetMinimum(),mixClon->GetMaximum(),1,3,2);
	
	x1Bin = D0_mixed->FindBin(x1);
	x2Bin = D0_mixed->FindBin(x2);
	nBins = x2Bin - x1Bin;
	



	// ******************************* Mixed Residual bg Subtraction *******************************
	fOutskirt->GetParameters(par);
	TF1 *fSignalMixed = new TF1("fSignalMixed",fitFunction,x1,x2,7);
	par[4]=40; par[5]=1.864; par[6]=0.011;
	fSignalMixed->SetLineColor(kRed);
	fSignalMixed->SetLineWidth(4);
	fSignalMixed->SetParameters(par);
	////////fSignalMixed->SetParLimits(4,0,1000);
	fSignalMixed->SetParLimits(5,1.850,fSignalLike->GetParameter(5));
	fSignalMixed->SetParLimits(6,0.007,fSignalLike->GetParameter(6));
	////fSignalMixed->SetNpx(1000);
//	fSignalMixed->FixParameter(4,0);
//	fSignalMixed->FixParameter(5,1.865);
//	fSignalMixed->FixParameter(6,0.0111);
	////fSignalMixed->FixParameter(0,par[0]);
	////fSignalMixed->FixParameter(1,par[1]);
	fSignalMixed->FixParameter(2,par[2]);
	fSignalMixed->FixParameter(3,par[3]);
	mixClon->Fit("fSignalMixed","","",x1,x2);
	/////fSignalMixed->Draw("same");

	TF1 *fsignalonlymix = new TF1("fsignalonlymix",normalDistribution,x1,x2,3);
	fsignalonlymix->FixParameter(0,fSignalMixed->GetParameter(4));
	fsignalonlymix->FixParameter(1,fSignalMixed->GetParameter(5));
	fsignalonlymix->FixParameter(2,fSignalMixed->GetParameter(6));

	fsignalonlymix->SetLineColor(kBlack);
	fsignalonlymix->SetLineWidth(4);
	fsignalonlymix->Draw("same");

	TF1 *fresbgonlymix = new TF1("fresbgonlymix",parabola,x1,x2,4);
	fresbgonlymix->FixParameter(0,fSignalMixed->GetParameter(0));
	fresbgonlymix->FixParameter(1,fSignalMixed->GetParameter(1));
	fresbgonlymix->FixParameter(2,fSignalMixed->GetParameter(2));
	fresbgonlymix->FixParameter(3,fSignalMixed->GetParameter(3));

	fresbgonlymix->SetLineColor(kRed);
	fresbgonlymix->SetLineWidth(4);
	fresbgonlymix->SetLineStyle(6);
	/////fresbgonlymix->Draw("same");

	signal = fSignalMixed->GetParameter(4)/mixClon->GetBinWidth(1);
	signalerror = fSignalMixed->GetParError(4)/mixClon->GetBinWidth(1);
	sprintf(parstring,"RwYld = %.0f #pm %.0f",signal,signalerror);
	c2->cd();
	TPad *infoMixed = new TPad("infoMixed","info",0.3,0.3,0.5,0.8);
	infoMixed->SetFillColor(10);
	/////infoMixed->Draw();
	infoMixed->cd();	
	
	tl.SetTextSize(0.12);
	tl.SetTextFont(42);
	tl.SetTextColor(kCyan+3);
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

	x1a = (fSignalMixed->GetParameter(5) - (3*fSignalMixed->GetParameter(6)));
	x2a = (fSignalMixed->GetParameter(5) + (3*fSignalMixed->GetParameter(6)));

	cout << x1a << "     " << x2a << endl;

	x1aBin = Unlike->FindBin(x1a);
	x2aBin = Unlike->FindBin(x2a);

	cout << x1aBin << "     " << x2aBin << endl;

	x1aBinLowEdge = Unlike->GetBinLowEdge(x1aBin);
	x2aBinHighEdge = Unlike->GetBinLowEdge(x2aBin) + Unlike->GetBinWidth(x2aBin);
	cout << x1aBinLowEdge << "    "  << x2aBinHighEdge << endl;

	signalgauss = (fsignalonlymix->Integral(x1aBinLowEdge,x2aBinHighEdge))/mixClon->GetBinWidth(1);
	cout << signalgauss << "     " << sqrt(Unlike->Integral(x1aBin, x2aBin)) << endl;

	sg = signalgauss/sqrt(Unlike->Integral(x1aBin, x2aBin));
    ////float sg = signalgauss/sqrt(signalgauss + 1*(fresbgonlymix->Integral(x1a,x2a)/mixClon->GetBinWidth(1)) + 1*(MixedBg->Integral(x1aBin, x2aBin)));
    cout<<"*****************************************************************************"<<endl;
	cout << "significance from mixed-event  =    " <<  sg << endl;
	////sprintf(parstring,"sg = %.3f",sg);
	////tl.DrawLatex(0.05,0.1,parstring);

	cout<<"*****************************************************************************"<<endl;
	
	////sprintf(parstring,"dm_{K#pi} = %.3f [GeV/c^{2}]",mixClon->GetBinWidth(1));
	////tl.DrawLatex(0.05,0.1,parstring);

	TCanvas *c20 = new TCanvas("c20");
	c20->cd();
	/////mixClon->DrawCopy();
	likeClon->DrawCopy("same");
	rotClon->DrawCopy("same");
	
	TLegend *resultsLegend = new TLegend(0.2,0.82,0.5,0.98,"","brNDC");
	resultsLegend->AddEntry(mixClon,"mixed event","p");
	resultsLegend->AddEntry(rotClon,"rotated momentum","p");
	resultsLegend->AddEntry(likeClon,"like sign","p");
	resultsLegend->Draw();
	//////infoMixed->DrawClone();
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
//	hSignal_Mixed->GetXaxis()->SetRangeUser(1.75,2.1);
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
	rotClon->GetYaxis()->SetTitle("Entries (/0.01 GeV/c^{2})");
	rotClon->GetYaxis()->SetTitleOffset(1.4);
	rotClon->GetYaxis()->SetLabelSize(0.05);
	rotClon->GetXaxis()->SetLabelSize(0.05);
	rotClon->GetYaxis()->SetTitleSize(0.05);
	rotClon->GetXaxis()->SetTitleSize(0.05);
	rotClon->GetYaxis()->SetLabelFont(42);
	rotClon->GetXaxis()->SetLabelFont(42);
	rotClon->GetYaxis()->SetTitleFont(42);
	rotClon->GetXaxis()->SetTitleFont(42);
	rotClon->SetMarkerSize(1);
	rotClon->SetLineWidth(2);
	//	rotClon->SetMarkerStyle(8);
	rotClon->SetMarkerColor(kGreen+1);
	///rotClon->Scale(2.);
	rotClon->Draw();
	infoRotate->SetPad(0.6,0.5,1.,0.9);
	infoRotate->SetFillColor(10);
	infoRotate->Draw("same");
	//////fsignalonlyrot->Draw("same");
	fresbgonlyrot->Draw("same");
	TCanvas *c11 = new TCanvas("c11","c11",400,400);
	c11->SetBottomMargin(0.15);
	c11->SetTickx();
	c11->SetTicky();
	c11->SetLeftMargin(0.15);
	likeClon->SetTitle("After Like-sign bg. subtracted");
	likeClon->GetXaxis()->SetRangeUser(1.73,2.05);
	likeClon->GetYaxis()->SetRangeUser(-1200,3500);
	likeClon->GetYaxis()->SetTitle("Entries (/0.01 GeV/c^{2})");
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
	likeClon->SetMarkerSize(1);
	likeClon->SetLineWidth(2);
	likeClon->Draw();
	infoLike->SetPad(0.6,0.5,1.,0.9);
	infoLike->SetFillColor(10);
	infoLike->Draw("same");
	//////fsignalonlyLike->Draw("same");
	fresbgonlyLike->Draw("same");


	TCanvas *c12 = new TCanvas("c12","c12",400,400);
	c12->SetBottomMargin(0.15);
	c12->SetLeftMargin(0.15);
	c12->SetTickx();
	c12->SetTicky();
//	hSignal_Mixed->SetTitle("Mixed event");
//	hSignal_Mixed->GetXaxis()->SetRangeUser(1.75,2.1);
//	hSignal_Mixed->GetYaxis()->SetRangeUser(-1200,2500);
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
//	c12->cd(2);
//	gPad->SetBottomMargin(0.2);
//	gPad->SetLeftMargin(0.15);
	mixClon->SetTitle("After Mixed-Event bg. subtracted");
	mixClon->GetXaxis()->SetRangeUser(1.73,2.05);
	mixClon->GetYaxis()->SetRangeUser(-750,2500);
	mixClon->GetYaxis()->SetTitle("Entries (/0.01 GeV/c^{2})");
	mixClon->GetYaxis()->SetTitleOffset(1.4);
	mixClon->GetYaxis()->SetLabelSize(0.05);
	mixClon->GetXaxis()->SetLabelSize(0.05);
	mixClon->GetYaxis()->SetTitleSize(0.05);
	mixClon->GetXaxis()->SetTitleSize(0.05);
	mixClon->GetYaxis()->SetLabelFont(42);
	mixClon->GetXaxis()->SetLabelFont(42);
	mixClon->GetYaxis()->SetTitleFont(42);
	mixClon->GetXaxis()->SetTitleFont(42);
	//	mixClon->SetMarkerStyle(8);
	mixClon->SetMarkerColor(kCyan+3);
	mixClon->SetMarkerSize(1);
	mixClon->SetLineWidth(2);
	///mixClon->Scale(2.);
	/////mixClon->Draw();
	infoMixed->SetPad(0.6,0.5,1.,0.9);
	infoMixed->SetFillColor(10);
	/////infoMixed->Draw("same");
	////fsignalonlymix->Draw("same");
	////fresbgonlymix->Draw("same");
//	
//	TF1 *signalFcn = new TF1("signalFcn",lorentzianPeak,0.8,1.05,3);
//	c4->cd();
//	signalFcn->SetParameters(1000,1,0.892);
//	D0_Like->Fit("signalFcn","","",0.8,0.95);
	


	//********************************** Draw a more plots into one **************************************

	c1->cd();
	TH1D *CloneR = (TH1D*)gDirectory->Get("CloneR");
	TH1D *CloneL = (TH1D*)gDirectory->Get("CloneL");
	TH1D *CloneM = (TH1D*)gDirectory->Get("CloneM");
	CloneR->Scale(2);
	CloneL->Scale(2);
	CloneM->Scale(2);
	CloneR->SetMarkerColor(8);
	CloneR->SetLineColor(8);
	CloneR->SetMarkerStyle(8);
	CloneL->SetMarkerColor(4);
	CloneL->SetLineColor(4);
	CloneL->SetMarkerStyle(8);
	CloneM->SetMarkerColor(kCyan+3);
	CloneM->SetLineColor(kCyan+3);
	CloneM->SetMarkerStyle(8);

	CloneL->Draw("same");
	CloneR->Draw("same");
	////CloneM->Draw("same");

	TGaxis* RawYieldAxis = new TGaxis(2.5,0,2.5,1e6,0,0.5e6,510,"+L");
	

	//+ : draw on positive side
	//L : left adjusted
	RawYieldAxis->SetName("RawYieldAxis");
	RawYieldAxis->SetLineColor(8);
	RawYieldAxis->SetTextColor(4);
	// RawYieldAxis->SetTitle("Raw Yield (/0.1 GeV/c^{2})");
	RawYieldAxis->SetLabelColor(4);
	RawYieldAxis->SetLabelFont(42);
	RawYieldAxis->Draw();

    ////leg->AddEntry(CloneM,"2*(Unlike Sign - Mixed Event)","p");
	leg->AddEntry(CloneR,"2*(Unlike Sign - Rotated Momentum)","p");
	leg->AddEntry(CloneL,"2*(Unlike Sign - Like Sign)","p");
	
	

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

