#include "style.C"
Double_t bgfitfunction(Double_t *x, Double_t *par)
	{
		if (x[0] > 0.1442 && x[0] < 0.1481) {
			////if (x[0] > 0.14 && x[0] < 0.15) {
			TF1::RejectPoint();
			return 0;
		}
		return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0];
//		return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
	}
Double_t linearbgfitfunction(Double_t *x, Double_t *par)
{
	if (x[0] > 0.1442 && x[0] < 0.1481) {
	////if (x[0] > 0.14 && x[0] < 0.15) {
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

void RawYieldPt_sg_dstar(float pT1, float pT2, int rebin) {
	style();
	float x1,x2,lowEdge,highEdge,normEdge1,normEdge2, normEdge3, normEdge4, ptBin1,ptBin2;
	int nBinsInvMass;
	bool GeometricMean = false;
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
	
	
	int rebinWrongSign = 1;
	int rebinSideBand = 1;
	int rebinMixed = 1;
	
    ptBin1 = DstarPlusRightSign->GetXaxis()->FindBin(pT1);
	ptBin2 = DstarPlusRightSign->GetXaxis()->FindBin(pT2);

	DstarPlusRightSign->ProjectionY("_py",ptBin1,ptBin2);
	DstarMinusRightSign->ProjectionY("_py",ptBin1,ptBin2);
	DstarPlusWrongSign->ProjectionY("_py",ptBin1,ptBin2);
	DstarMinusWrongSign->ProjectionY("_py",ptBin1,ptBin2);
	DstarPlusSideBand->ProjectionY("_py",ptBin1,ptBin2);
	DstarMinusSideBand->ProjectionY("_py",ptBin1,ptBin2);
	MixedDstarPlusRightSign->ProjectionY("_py",ptBin1,ptBin2);
	MixedDstarMinusRightSign->ProjectionY("_py",ptBin1,ptBin2);
	
  TH1D *DstarPlusRightSign_py = (TH1D*)gDirectory->Get("DstarPlusRightSign_py");
  TH1D *DstarMinusRightSign_py = (TH1D*)gDirectory->Get("DstarMinusRightSign_py");
  TH1D *DstarPlusWrongSign_py = (TH1D*)gDirectory->Get("DstarPlusWrongSign_py");
  TH1D *DstarMinusWrongSign_py = (TH1D*)gDirectory->Get("DstarMinusWrongSign_py");
  TH1D *DstarPlusSideBand_py = (TH1D*)gDirectory->Get("DstarPlusSideBand_py");
  TH1D *DstarMinusSideBand_py = (TH1D*)gDirectory->Get("DstarMinusSideBand_py");
  TH1D *MixedDstarPlusRightSign_py = (TH1D*)gDirectory->Get("MixedDstarPlusRightSign_py");
  TH1D *MixedDstarMinusRightSign_py = (TH1D*)gDirectory->Get("MixedDstarMinusRightSign_py");


    DstarPlusRightSign_py->Sumw2();
	DstarMinusRightSign_py->Sumw2();

	DstarPlusWrongSign_py->Sumw2();
	DstarMinusWrongSign_py->Sumw2();
	DstarPlusSideBand_py->Sumw2();
	DstarMinusSideBand_py->Sumw2();
	MixedDstarPlusRightSign_py->Sumw2();
	MixedDstarMinusRightSign_py->Sumw2();
	
	DstarPlusRightSign_py->Rebin(rebin);
	DstarMinusRightSign_py->Rebin(rebin);

	DstarPlusWrongSign_py->Rebin(rebin);
	DstarMinusWrongSign_py->Rebin(rebin);
	DstarPlusSideBand_py->Rebin(rebin);
	DstarMinusSideBand_py->Rebin(rebin);
	MixedDstarPlusRightSign_py->Rebin(rebin);
	MixedDstarMinusRightSign_py->Rebin(rebin);
	
	nBinsInvMass = DstarPlusRightSign_py->GetNbinsX();
	lowEdge = DstarPlusRightSign_py->GetBinLowEdge(1);
	highEdge = DstarPlusRightSign_py->GetBinLowEdge(nBinsInvMass+1);
	
	cout<<nBinsInvMass<<"   "<<lowEdge<<"    "<<highEdge<<endl;
	
	TLegend *leg = new TLegend(0.6,0.7,0.88,0.88,"");
	leg->SetFillColor(10);
	leg->SetTextFont(42);
    TH1D *RightSign = new TH1D("RightSign","",nBinsInvMass,lowEdge,highEdge);
	
	TH1D *WrongSignBg= new TH1D("WrongSignBg","",nBinsInvMass,lowEdge,highEdge);
	TH1D *SideBandBg= new TH1D("SideBandBg","",nBinsInvMass,lowEdge,highEdge);
	TH1D *MixedBg = new TH1D("MixedBg","",nBinsInvMass,lowEdge,highEdge);
	
	RightSign->Sumw2();
	
	WrongSignBg->Sumw2();
	SideBandBg->Sumw2();
	MixedBg->Sumw2();

    RightSign->Add(DstarPlusRightSign_py,DstarMinusRightSign_py,1,1);
	
	if(GeometricMean){
		float a,b,da,db;
		for(int i=0; i<=nBinsInvMass+1; i++){
			a = DstarPlusWrongSign_py->GetBinContent(i);
			b = DstarMinusWrongSign_py->GetBinContent(i);
			da = DstarPlusWrongSign_py->GetBinError(i);
			db = DstarMinusWrongSign_py->GetBinError(i);	
			if(a==0 || b==0) 
			{
				WrongSignBg->SetBinContent(i,2*TMath::Sqrt(a*b));
				WrongSignBg->SetBinError(i,0);
			}else{
				WrongSignBg->SetBinContent(i,2*TMath::Sqrt(a*b));
				WrongSignBg->SetBinError(i,TMath::Sqrt(a*b)*TMath::Sqrt(((da/a)*(da/a)) + ((db/b)*(db/b))));
			}
		}
	} else WrongSignBg->Add(DstarPlusWrongSign_py,DstarMinusWrongSign_py,1,1);
	SideBandBg->Add(DstarPlusSideBand_py,DstarMinusSideBand_py,1,1);
	MixedBg->Add(MixedDstarPlusRightSign_py,MixedDstarMinusRightSign_py,1,1);
	
	Double_t Scalefactorright;
	Double_t Scalefactorleft;
	Double_t Scalefactor;
	
	normEdge1 = RightSign->FindBin(0.138);
	normEdge2 = RightSign->FindBin(0.142);
	
	normEdge3 = RightSign->FindBin(0.15);
	normEdge4 = RightSign->FindBin(0.16);

	Scalefactorleft = RightSign->Integral(normEdge1,normEdge2)/WrongSignBg->Integral(normEdge1,normEdge2);
	Scalefactorright = RightSign->Integral(normEdge3,normEdge4)/WrongSignBg->Integral(normEdge3,normEdge4);
	Scalefactor = (Scalefactorleft + Scalefactorright)/2.0;
	WrongSignBg->Scale(Scalefactorright);
    cout<<"WrongSign scale factor: "<<Scalefactor<<endl;


	Scalefactorleft = RightSign->Integral(normEdge1,normEdge2)/SideBandBg->Integral(normEdge1,normEdge2);
	Scalefactorright = RightSign->Integral(normEdge3,normEdge4)/SideBandBg->Integral(normEdge3,normEdge4);
	Scalefactor = (Scalefactorleft + Scalefactorright)/2.0; 
	
	///////Scalefactor = 0.9999;
	SideBandBg->Scale(Scalefactorright);
  cout<<"SideBand scale factor: "<<Scalefactor<<endl;
	
	Scalefactorleft = RightSign->Integral(normEdge1,normEdge2)/MixedBg->Integral(normEdge1,normEdge2);
	Scalefactorright = RightSign->Integral(normEdge3,normEdge4)/MixedBg->Integral(normEdge3,normEdge4);
	
	Scalefactor = (Scalefactorleft + Scalefactorright)/2.0;
	
	///////MixedBg->Scale(Scalefactor);
	cout<<"Mixed event Scale factor: "<<Scalefactor<<endl;
	
	
	//********************************************* Draw Invariant mass **********************************************
	TCanvas *c1 = new TCanvas("c1");
	c1->cd();
	c1->SetRightMargin(0.05);
	RightSign->GetXaxis()->SetTitle("M_{K#pi#pi} - M_{K#pi} [GeV/c^{2}]");
	RightSign->GetYaxis()->SetTitle(Form("Entries (/%.4f GeV/c^{2})", 0.0002*rebin));
	RightSign->SetMarkerStyle(4);
	RightSign->GetYaxis()->SetTitleFont(42);
	RightSign->GetYaxis()->SetLabelFont(42);
	//////RightSign->GetYaxis()->SetRangeUser(-100, 300);
	RightSign->Draw("");



	WrongSignBg->SetMarkerStyle(7);
	WrongSignBg->SetMarkerColor(4);
	WrongSignBg->SetLineColor(4);
	WrongSignBg->SetLineWidth(2);
	WrongSignBg->Draw("hist C same");
	SideBandBg->SetMarkerStyle(7);
	SideBandBg->SetMarkerColor(8);
	SideBandBg->SetLineColor(8);
	SideBandBg->SetLineWidth(2);
	SideBandBg->Draw("hist C same");
	
	MixedBg->SetMarkerStyle(7);
	MixedBg->SetMarkerColor(kCyan+3);
	MixedBg->SetLineColor(kCyan+3);
	MixedBg->SetLineWidth(2);
	//////////MixedBg->Draw("hist C same");
    
	leg->AddEntry(RightSign,"Right Sign","pf");
	leg->AddEntry(WrongSignBg,"Wrong Sign","l");
	leg->AddEntry(SideBandBg,"Side Band","l");
	//////////leg->AddEntry(MixedBg,"Mixed Event","l");

	// **************************************************************************************************************
	// ********************************************* Wrong sign ***********************************************
	// **************************************************************************************************************	
	x1 = 0.138;
	x2 = 0.154;
//    x1=1.14;
//    x2=1.15;
	TCanvas *c4 = new TCanvas("c4");
	c4->cd();
	c4->SetRightMargin(0.05);
  TH1D *D0_WrongSign;
	if(GeometricMean) 	D0_WrongSign = new TH1D("D0_WrongSign","N_{K^{-}#pi^{+}#pi^{+}}+N_{K^{+}#pi^{-}#pi^{-}} - 2 #sqrt{N_{K^{-}#pi^{+}#pi^{-}}N_{K^{+}#pi^{-}#pi^{+}}}",nBinsInvMass,lowEdge,highEdge);
	else D0_WrongSign = new TH1D("D0_WrongSign","N_{K^{-}#pi^{+}#pi^{+}}+N_{K^{+}#pi^{-}#pi^{-}} - (N_{K^{-}#pi^{+}#pi^{-}}+N_{K^{+}#pi^{-}#pi^{+}})",nBinsInvMass,lowEdge,highEdge);

	D0_WrongSign->Sumw2();
	for(int i=0; i<= D0_WrongSign->GetNbinsX()+1; i++){
		D0_WrongSign->SetBinContent(i,RightSign->GetBinContent(i) - WrongSignBg->GetBinContent(i));
		D0_WrongSign->SetBinError(i,TMath::Sqrt(((RightSign->GetBinError(i)*RightSign->GetBinError(i)) +
											   (WrongSignBg->GetBinError(i)*WrongSignBg->GetBinError(i)))));
	}	
	D0_WrongSign->Draw("bar");
	D0_WrongSign->GetXaxis()->SetTitle("M_{K#pi#pi} - M_{K#pi} [GeV/c^{2}]");
	D0_WrongSign->GetYaxis()->SetTitle("counts");
	D0_WrongSign->SetFillColor(4);
//	D0_WrongSign->SetMarkerStyle(7);
    D0_WrongSign->SetLineColor(4);
	D0_WrongSign->SetMarkerColor(4);

    TF1 *fOutskirt = new TF1("fOutskirt",bgfitfunction,x1,x2,4);
	/////fOutskirt->SetParameter(2,0);
	/////fOutskirt->SetParameter(1,-1.77653e+03);
	////fOutskirt->SetParameter(0,2.48347e+02);
    fOutskirt->FixParameter(2,0.);
	fOutskirt->FixParameter(3,0.);

	D0_WrongSign->Fit("fOutskirt","N","",x1,x2);
	cout<<"*****************************************************************************"<<endl;
	cout<<"WrongSign residual bg chisquare/ndf="<<fOutskirt->GetChisquare()<<"/"<<fOutskirt->GetNDF()<<endl;
	cout<<"*****************************************************************************"<<endl;

	
	TPad *zoom1WrongSign = new TPad("zoom1WrongSign","zoom",0.5,0.25,0.94,0.89);
	zoom1WrongSign->SetFillColor(10);
	zoom1WrongSign->Draw();
	zoom1WrongSign->cd();
	zoom1WrongSign->SetRightMargin(0.05);
	D0_WrongSign->Clone("WrongSignClone");
	D0_WrongSign->Clone("CloneWS");

	TH1D *WrongSignClon = (TH1D*)gDirectory->Get("WrongSignClone");
//	for (int i=0; i<lowborder; i=i+2) {
//		float content1 = WrongSignClon->GetBinContent(i+1);
//		float content2 = WrongSignClon->GetBinContent(i+2);
//		float error1 = WrongSignClon->GetBinError(i+1);
//		float error2 = WrongSignClon->GetBinError(i+2);
//		WrongSignClon->SetBinContent(i+1,(float)(content1+content2)/2);
//		WrongSignClon->SetBinContent(i+2,0);
//		WrongSignClon->SetBinError(i+1,TMath::Sqrt(1/(1/(error1*error1)+1/(error2*error2))));
//		WrongSignClon->SetBinError(i+2,0);
//	}
//	for (int i=highborder; i<WrongSignClon->GetNbinsX(); i=i+2) {
//		float content1 = WrongSignClon->GetBinContent(i+1);
//		float content2 = WrongSignClon->GetBinContent(i+2);
//		float error1 = WrongSignClon->GetBinError(i+1);
//		float error2 = WrongSignClon->GetBinError(i+2);
//		WrongSignClon->SetBinContent(i+1,(float)(content1+content2)/2);
//		WrongSignClon->SetBinContent(i+2,0);
//		WrongSignClon->SetBinError(i+1,TMath::Sqrt(1/(1/(error1*error1)+1/(error2*error2))));
//		WrongSignClon->SetBinError(i+2,0);
//	}
	WrongSignClon->SetTitle("");
	WrongSignClon->Draw();
	WrongSignClon->GetXaxis()->SetTitle("M_{K#pi#pi} - M_{K#pi} [GeV/c^{2}]");
	WrongSignClon->SetMarkerStyle(8);
	WrongSignClon->GetXaxis()->SetRangeUser(x1,x2);
	PlotLine(0.14,0.14,WrongSignClon->GetMinimum(),WrongSignClon->GetMaximum(),1,3,2);
	PlotLine(0.15,0.15,WrongSignClon->GetMinimum(),WrongSignClon->GetMaximum(),1,3,2);
	
	int x1Bin = D0_WrongSign->FindBin(x1);
	int x2Bin = D0_WrongSign->FindBin(x2);
	int nBins = x2Bin - x1Bin;
	
	// ******************************* WrongSign Residual bg Subtraction *******************************
	fOutskirt->GetParameters(par);
	TF1 *fSignalWrongSign = new TF1("fSignalWrongSign",fitFunction,x1,x2,7);
	fSignalWrongSign->SetLineColor(kRed);
	fSignalWrongSign->SetLineWidth(6);
	par[4]=0.304228; par[5]=0.145558; par[6]=0.000809652;
	////////par[4]=2.52872e-01; par[5]=1.45468e-01; par[6]=4.06005e-04;
	fSignalWrongSign->SetParameters(par);
	///////fSignalWrongSign->SetParLimits(4,0,1000);
	fSignalWrongSign->SetParLimits(5,0.143,0.147);
	fSignalWrongSign->SetParLimits(6,0.0001,0.0015);
	/*fSignalWrongSign->SetParLimits(5,0.1448,0.1459);
	fSignalWrongSign->SetParLimits(6,0.00015,0.002);*/
	fSignalWrongSign->SetNpx(1000);
//	fSignalWrongSign->FixParameter(4,0);	
//	fSignalWrongSign->FixParameter(5,0.145);
//	fSignalWrongSign->FixParameter(6,0.0003);
//	fSignalWrongSign->FixParameter(0,par[0]);
//	fSignalWrongSign->FixParameter(1,par[1]);
    fSignalWrongSign->FixParameter(2,par[2]);
	fSignalWrongSign->FixParameter(3,par[3]);
	WrongSignClon->Fit("fSignalWrongSign","","",x1,x2);
	fSignalWrongSign->Draw("same");


	TF1 *fsignalonlyWrongSign = new TF1("fsignalonlyWrongSign",normalDistribution,x1,x2,3);
	fsignalonlyWrongSign->FixParameter(0,fSignalWrongSign->GetParameter(4));
	fsignalonlyWrongSign->FixParameter(1,fSignalWrongSign->GetParameter(5));
	fsignalonlyWrongSign->FixParameter(2,fSignalWrongSign->GetParameter(6));

	fsignalonlyWrongSign->SetLineColor(kBlack);
	fsignalonlyWrongSign->SetLineWidth(6);
	fsignalonlyWrongSign->SetLineStyle(10);
	fsignalonlyWrongSign->Draw("same");

	TF1 *fresbgonlyWrongSign = new TF1("fresbgonlyWrongSign",parabola,x1,x2,4);
	fresbgonlyWrongSign->FixParameter(0,fSignalWrongSign->GetParameter(0));
	fresbgonlyWrongSign->FixParameter(1,fSignalWrongSign->GetParameter(1));
	fresbgonlyWrongSign->FixParameter(2,fSignalWrongSign->GetParameter(2));
	fresbgonlyWrongSign->FixParameter(3,fSignalWrongSign->GetParameter(3));

	fresbgonlyWrongSign->SetLineColor(kRed);
	fresbgonlyWrongSign->SetLineWidth(6);
	fresbgonlyWrongSign->SetLineStyle(6);
	fresbgonlyWrongSign->Draw("same");
	
	float signal = fSignalWrongSign->GetParameter(4)/WrongSignClon->GetBinWidth(1);
	float signalerror = fSignalWrongSign->GetParError(4)/WrongSignClon->GetBinWidth(1);
	sprintf(parstring,"RwYld = %.0f #pm %.0f",signal,signalerror);
	c4->cd();
	TPad *infoWrongSign = new TPad("infoWrongSign","info",0.3,0.3,0.5,0.8);
	infoWrongSign->SetFillColor(10);
	infoWrongSign->SetFillStyle(4000);
	infoWrongSign->Draw();
	infoWrongSign->cd();

	TLatex tl;	
	tl.SetTextSize(0.12);
    tl.SetTextFont(42);
	tl.SetTextColor(4);
	tl.DrawLatex(0.05,0.85,parstring);
	sprintf(parstring,"#chi^{2}/ndf = %.3f/%i",fSignalWrongSign->GetChisquare(),fSignalWrongSign->GetNDF());	
	tl.DrawLatex(0.05,0.7,parstring);
	/*sprintf(parstring,"#mu = %.3f #pm %.4f",fSignalWrongSign->GetParameter(5),fSignalWrongSign->GetParError(5));	
	tl.DrawLatex(0.05,0.55,parstring);
	sprintf(parstring,"#sigma = %.4f #pm %.4f",fSignalWrongSign->GetParameter(6),fSignalWrongSign->GetParError(6));	
	tl.DrawLatex(0.05,0.4,parstring);*/
	tl.SetTextColor(1);
	sprintf(parstring,"p_{T} #in [%.1f,%.1f] (GeV/c)",pT1,pT2);
	tl.DrawLatex(0.05,0.55,parstring);

	float x1a = (fSignalWrongSign->GetParameter(5) - (3*fSignalWrongSign->GetParameter(6)));
	float x2a = (fSignalWrongSign->GetParameter(5) + (3*fSignalWrongSign->GetParameter(6)));

	cout << x1a << "     " << x2a << endl;

	int x1aBin = RightSign->FindBin(x1a);
	int x2aBin = RightSign->FindBin(x2a);

	cout << x1aBin << "     " << x2aBin << endl;

	float x1aBinLowEdge = RightSign->GetBinLowEdge(x1aBin);
	float x2aBinHighEdge = RightSign->GetBinLowEdge(x2aBin) + RightSign->GetBinWidth(x2aBin);

	cout << x1aBinLowEdge << "    "  << x2aBinHighEdge << endl;

	float signalgauss = (fsignalonlyWrongSign->Integral(x1aBinLowEdge,x2aBinHighEdge))/WrongSignClon->GetBinWidth(1);
	cout << signalgauss << "     " << sqrt(RightSign->Integral(x1aBin, x2aBin)) << endl;

	float sg = signalgauss/sqrt(RightSign->Integral(x1aBin, x2aBin));
	///sg = signalgauss/sqrt(signalgauss + 1*(fresbgonlyWrongSign->Integral(x1a,x2a)/WrongSignClon->GetBinWidth(1)) + 1*(WrongSignBg->Integral(x1aBin, x2aBin)));
	sprintf(parstring,"Significance = %.1f",sg);
	tl.DrawLatex(0.05,0.4,parstring);
    
	cout<<"*****************************************************************************"<<endl;
	cout << "significance from WrongSign  =    " <<  sg << endl;
	cout<<"*****************************************************************************"<<endl;

	const char* inputFile = "wrongsign_table.tex";
    const char* tempFile = "temp_wrongsign_table.txt";

	std::ifstream infile(inputFile);
    if (!infile.is_open()) {
        std::cerr << "Error: could not open file " << inputFile << std::endl;
        return;
    }

	std::ofstream outfile(tempFile);
    if (!outfile.is_open()) {
        std::cerr << "Error: could not open file " << tempFile << std::endl;
        return;
    }

	std::string line;
    int lineNumber = 0;


    while (getline(infile, line)) {
        lineNumber++;

        // Modify line 5 by appending "pT"
        if (lineNumber == 5) {
            line += Form(" & %.1f  $ < p_{T} < $  %.1f", pT1, pT2);
        }

		if (lineNumber == 7) {
			line += Form(" & %.1f/%i", fSignalWrongSign->GetChisquare(), fSignalWrongSign->GetNDF());
		}

		if (lineNumber == 8) {
			line += Form(" & %.0f $\\pm$ %.0f", signal, signalerror);
		}

		if (lineNumber == 9) {
			line += Form(" & %.1f $\\pm$ %.1f", 1000*fSignalWrongSign->GetParameter(5), 1000*fSignalWrongSign->GetParError(5));
		}

		if (lineNumber == 10) {
			line += Form(" & %.2f $\\pm$ %.2f", 1000*fSignalWrongSign->GetParameter(6), 1000*fSignalWrongSign->GetParError(6));
		}

		if (lineNumber == 11) {
			line += Form(" & %.1f", sg);
		}

		// Write the modified (or unmodified) line to the temporary file
		outfile << line << std::endl;
	}

	// Close files
	infile.close();
	outfile.close();

	// Replace the original file with the modified one
	std::remove(inputFile);            // Delete the original file
	std::rename(tempFile, inputFile); // Rename the temporary file to the original name
	

	////sprintf(parstring,"dm_{K#pi} = %.3f [GeV/c^{2}]",WrongSignClon->GetBinWidth(1));
	////tl.DrawLatex(0.05,0.1,parstring);
	
	
	// **************************************************************************************************************
	// ********************************************* Side Band ***********************************************
	// **************************************************************************************************************	
	x1 = 0.138;
	x2 = 0.154;
//    x1 = 1.8;
//    x2 = 1.92;
	TCanvas *c3 = new TCanvas("c3");
	c3->cd();
	c3->SetRightMargin(0.05);
		
	TH1D *Dstar_SideBand = new TH1D("Dstar_SideBand","RightSign - SideBand",nBinsInvMass,lowEdge,highEdge);
	Dstar_SideBand->Sumw2();
	for(int i=0; i<= Dstar_SideBand->GetNbinsX()+1; i++){
		Dstar_SideBand->SetBinContent(i,RightSign->GetBinContent(i) - SideBandBg->GetBinContent(i));
		Dstar_SideBand->SetBinError(i,TMath::Sqrt(RightSign->GetBinError(i)*RightSign->GetBinError(i) +
											  SideBandBg->GetBinError(i)*SideBandBg->GetBinError(i)));
	}	
	Dstar_SideBand->Draw("bar");
	Dstar_SideBand->GetXaxis()->SetTitle("M_{K#pi#pi} - M_{K#pi} [GeV/c^{2}]");
	Dstar_SideBand->GetYaxis()->SetTitle("counts");
	//Dstar_SideBand->SetMarkerStyle(7);
	Dstar_SideBand->SetMarkerColor(kGreen+2);
	Dstar_SideBand->SetLineColor(kGreen+2);
	Dstar_SideBand->SetFillColor(kGreen+2);
	Dstar_SideBand->Fit("fOutskirt","N","",x1,x2);
	cout<<"*****************************************************************************"<<endl;
	cout<<"SideBand residual bg chisquare/ndf="<<fOutskirt->GetChisquare()<<"/"<<fOutskirt->GetNDF()<<endl;
	cout<<"*****************************************************************************"<<endl;
	
	TPad *zoom1SideBand = new TPad("zoom1SideBand","zoom",0.5,0.25,0.94,0.89);
	zoom1SideBand->SetFillColor(10);
	zoom1SideBand->Draw();
	zoom1SideBand->cd();
	zoom1SideBand->SetRightMargin(0.05);
	Dstar_SideBand->Clone("SideBandClone");
	Dstar_SideBand->Clone("CloneSB");
	TH1D *SideBandClon = (TH1D*)gDirectory->Get("SideBandClone");
	int lowborder = SideBandClon->FindBin(0.142);
	int highborder = SideBandClon->FindBin(0.148);
//	for (int i=0; i<lowborder; i=i+2) {
//		float content1 = SideBandClon->GetBinContent(i+1);
//		float content2 = SideBandClon->GetBinContent(i+2);
//		float error1 = SideBandClon->GetBinError(i+1);
//		float error2 = SideBandClon->GetBinError(i+2);
//		SideBandClon->SetBinContent(i+1,(float)(content1+content2)/2);
//		SideBandClon->SetBinContent(i+2,0);
//		SideBandClon->SetBinError(i+1,TMath::Sqrt(1/(1/(error1*error1)+1/(error2*error2))));
//		SideBandClon->SetBinError(i+2,0);
//	}
//	for (int i=highborder; i<SideBandClon->GetNbinsX(); i=i+2) {
//		float content1 = SideBandClon->GetBinContent(i+1);
//		float content2 = SideBandClon->GetBinContent(i+2);
//		float error1 = SideBandClon->GetBinError(i+1);
//		float error2 = SideBandClon->GetBinError(i+2);
//		SideBandClon->SetBinContent(i+1,(float)(content1+content2)/2);
//		SideBandClon->SetBinContent(i+2,0);
//		SideBandClon->SetBinError(i+1,TMath::Sqrt(1/(1/(error1*error1)+1/(error2*error2))));
//		SideBandClon->SetBinError(i+2,0);
//	}
	
	SideBandClon->SetTitle("");
	SideBandClon->Draw();
	SideBandClon->GetXaxis()->SetTitle("M_{K#pi#pi} - M_{K#pi} [GeV/c^{2}]");
	SideBandClon->SetMarkerStyle(8);
	SideBandClon->GetXaxis()->SetRangeUser(x1,x2);
	PlotLine(0.140,0.140,SideBandClon->GetMinimum(),SideBandClon->GetMaximum(),1,3,2);
	PlotLine(0.148,0.148,SideBandClon->GetMinimum(),SideBandClon->GetMaximum(),1,3,2);
	
	x1Bin = Dstar_SideBand->FindBin(x1);
	x2Bin = Dstar_SideBand->FindBin(x2);
	nBins = x2Bin - x1Bin;
	
	// ******************************* SideBand Residual bg Subtraction *******************************
	fOutskirt->GetParameters(par);
	TF1 *fSignalSideBand = new TF1("fSignalSideBand",fitFunction,x1,x2,7);
	par[4]=0.304228; par[5]=0.145558; par[6]=0.000809652;
	fSignalSideBand->SetLineColor(kRed);
	fSignalSideBand->SetLineWidth(6);
	fSignalSideBand->SetParameters(par);
	/////fSignalSideBand->SetParLimits(4,0,1000);
	/////cout << fSignalSideBand->GetParameter(6) << "     " << fSignalWrongSign->GetParameter(6) << endl;
	fSignalSideBand->SetParLimits(5,0.143,0.147);
	fSignalSideBand->SetParLimits(6,0.0001,0.0015); ////0.0001,0.0015
	fSignalSideBand->SetNpx(1000);
//	fSignalSideBand->FixParameter(4,0); // zero signal hypothesis
//	fSignalSideBand->FixParameter(5,1.864);
//	fSignalSideBand->FixParameter(6,0.10);
//	fSignalSideBand->FixParameter(0,par[0]);
    ////fSignalSideBand->FixParameter(4,par[4]);
	fSignalSideBand->FixParameter(2,par[2]);
	fSignalSideBand->FixParameter(3,par[3]);
	SideBandClon->Fit("fSignalSideBand","","",x1,x2);
	fSignalSideBand->Draw("same");

	TF1 *fsignalonlySideBand = new TF1("fsignalonlySideBand",normalDistribution,x1,x2,3);
	fsignalonlySideBand->FixParameter(0,fSignalSideBand->GetParameter(4));
	fsignalonlySideBand->FixParameter(1,fSignalSideBand->GetParameter(5));
	fsignalonlySideBand->FixParameter(2,fSignalSideBand->GetParameter(6));

	fsignalonlySideBand->SetLineColor(kBlack);
	fsignalonlySideBand->SetLineStyle(10);
	fsignalonlySideBand->SetLineWidth(6);
	fsignalonlySideBand->Draw("same");

	TF1 *fresbgonlySideBand = new TF1("fresbgonlySideBand",parabola,x1,x2,4);
	fresbgonlySideBand->FixParameter(0,fSignalSideBand->GetParameter(0));
	fresbgonlySideBand->FixParameter(1,fSignalSideBand->GetParameter(1));
	fresbgonlySideBand->FixParameter(2,fSignalSideBand->GetParameter(2));
	fresbgonlySideBand->FixParameter(3,fSignalSideBand->GetParameter(3));

	fresbgonlySideBand->SetLineColor(kRed);
	fresbgonlySideBand->SetLineWidth(6);
	fresbgonlySideBand->SetLineStyle(6);
	///////fresbgonlySideBand->Draw("same");
	

	signal = fSignalSideBand->GetParameter(4)/SideBandClon->GetBinWidth(1);
	signalerror = fSignalSideBand->GetParError(4)/SideBandClon->GetBinWidth(1);
	sprintf(parstring,"RwYld = %.0f #pm %.0f",signal,signalerror);
	c3->cd();
	TPad *infoSideBand = new TPad("infoSideBand","info",0.3,0.3,0.5,0.8);
	infoSideBand->SetFillColor(10);
	infoSideBand->SetFillStyle(4000);
	infoSideBand->Draw();
	infoSideBand->cd();
	tl.SetTextSize(0.12);
	tl.SetTextColor(kGreen+2);
    tl.SetTextFont(42);
	tl.DrawLatex(0.05,0.85,parstring);
	sprintf(parstring,"#chi^{2}/ndf = %.3f/%i",fSignalSideBand->GetChisquare(),fSignalSideBand->GetNDF());	
	tl.DrawLatex(0.05,0.7,parstring);
	/*sprintf(parstring,"#mu = %.3f #pm %.4f",fSignalSideBand->GetParameter(5),fSignalSideBand->GetParError(5));	
	tl.DrawLatex(0.05,0.55,parstring);
	sprintf(parstring,"#sigma = %.4f #pm %.4f",fSignalSideBand->GetParameter(6),fSignalSideBand->GetParError(6));	
	tl.DrawLatex(0.05,0.4,parstring);*/
	tl.SetTextColor(1);
	sprintf(parstring,"p_{T} #in [%.1f,%.1f] (GeV/c)",pT1,pT2);
	////tl.DrawLatex(0.05,0.25,parstring);
	tl.DrawLatex(0.05,0.55,parstring);

	x1a = (fSignalSideBand->GetParameter(5) - (3*fSignalSideBand->GetParameter(6)));
	x2a = (fSignalSideBand->GetParameter(5) + (3*fSignalSideBand->GetParameter(6)));

	cout << x1a << "     " << x2a << endl;

	x1aBin = RightSign->FindBin(x1a);
	x2aBin = RightSign->FindBin(x2a);

	cout << x1aBin << "     " << x2aBin << endl;

	x1aBinLowEdge = RightSign->GetBinLowEdge(x1aBin);
	x2aBinHighEdge = RightSign->GetBinLowEdge(x2aBin) + RightSign->GetBinWidth(x2aBin);

	cout << x1aBinLowEdge << "    "  << x2aBinHighEdge << endl;

    signalgauss = fsignalonlySideBand->Integral(x1aBinLowEdge,x2aBinHighEdge)/SideBandClon->GetBinWidth(1);

	cout << signalgauss << "     " << sqrt(RightSign->Integral(x1aBin, x2aBin)) << endl;

	sg = signalgauss/sqrt(RightSign->Integral(x1aBin, x2aBin));
	sprintf(parstring,"Significance = %.1f",sg);
	////tl.DrawLatex(0.05,0.10,parstring);
	tl.DrawLatex(0.05,0.40,parstring);
	///sg = signalgauss/sqrt(signalgauss + 1*(fresbgonlySideBand->Integral(x1a,x2a)/SideBandClon->GetBinWidth(1)) + 1*(SideBandBg->Integral(x1aBin, x2aBin)));
    cout<<"*****************************************************************************"<<endl;
	cout << "significance from SideBand  =    " <<  sg << endl;
	////sprintf(parstring,"sg = %.3f",sg);
	////tl.DrawLatex(0.05,0.1,parstring);

	cout<<"*****************************************************************************"<<endl;

	inputFile = "sideband_table.tex";
	tempFile = "temp_sideband_table.txt";

	infile.open(inputFile);
	if (!infile.is_open()) {
		std::cerr << "Error: could not open file " << inputFile << std::endl;
		return;
	}

	outfile.open(tempFile);
	if (!outfile.is_open()) {
		std::cerr << "Error: could not open file " << tempFile << std::endl;
		return;
	}

	lineNumber = 0;

	while (getline(infile, line)) {
        lineNumber++;

        // Modify line 5 by appending "pT"
        if (lineNumber == 5) {
            line += Form(" & %.1f  $ < p_{T} < $  %.1f", pT1, pT2);
        }

		if (lineNumber == 7) {
			line += Form(" & %.1f/%i", fSignalSideBand->GetChisquare(), fSignalSideBand->GetNDF());
		}

		if (lineNumber == 8) {
			line += Form(" & %.0f $\\pm$ %.0f", signal, signalerror);
		}

		if (lineNumber == 9) {
			line += Form(" & %.1f $\\pm$ %.1f", 1000*fSignalSideBand->GetParameter(5), 1000*fSignalSideBand->GetParError(5));
		}

		if (lineNumber == 10) {
			line += Form(" & %.2f $\\pm$ %.2f", 1000*fSignalSideBand->GetParameter(6), 1000*fSignalSideBand->GetParError(6));
		}

		if (lineNumber == 11) {
			line += Form(" & %.1f", sg);
		}

		// Write the modified (or unmodified) line to the temporary file
		outfile << line << std::endl;
	}

	// Close files
	infile.close();
	outfile.close();

	// Replace the original file with the modified one
	std::remove(inputFile);            // Delete the original file
	std::rename(tempFile, inputFile); // Rename the temporary file to the original name
	
	
	
	// **************************************************************************************************************
	// ************************************************ Mixed Event *************************************************
	// **************************************************************************************************************
	x1 = 0.138;
	x2 = 0.154;	
	nBinsInvMass = MixedBg->GetNbinsX();
	lowEdge = MixedBg->GetBinLowEdge(1);
	highEdge = MixedBg->GetBinLowEdge(nBinsInvMass+1);

	TCanvas *c2 = new TCanvas("c2");
	c2->cd();
	c2->SetRightMargin(0.05);

	TH1D *D0_mixed = new TH1D("D0_mixed","N_{K^{-}#pi^{+}}+N_{K^{+}#pi^{-}} - (N_{K^{-}#pi^{+}}+N_{K^{+}#pi^{-}})_{buffer event}",nBinsInvMass,lowEdge,highEdge);
	D0_mixed->Sumw2();
	for(int i=0; i<= D0_mixed->GetNbinsX()+1; i++){
		D0_mixed->SetBinContent(i,RightSign->GetBinContent(i) - MixedBg->GetBinContent(i));
		D0_mixed->SetBinError(i,TMath::Sqrt(((RightSign->GetBinError(i)*RightSign->GetBinError(i)) +
											  (MixedBg->GetBinError(i)*MixedBg->GetBinError(i)))));
	}	
	//////////D0_mixed->Draw("bar");
	D0_mixed->GetXaxis()->SetTitle("M_{K#pi#pi} - M_{K#pi} [GeV/c^{2}]");
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
	//////////zoom1Mixed->Draw();
	zoom1Mixed->cd();
	zoom1Mixed->SetRightMargin(0.05);
	D0_mixed->Clone("mixedClone");
	D0_mixed->Clone("CloneM");


	TH1D *mixClon = (TH1D*)gDirectory->Get("mixedClone");
	mixClon->SetTitle("");
	mixClon->Draw();
	mixClon->GetXaxis()->SetTitle("M_{K#pi#pi} - M_{K#pi} [GeV/c^{2}]");
	mixClon->SetMarkerStyle(8);
	mixClon->GetXaxis()->SetRangeUser(x1,x2);
	PlotLine(0.140,0.140,mixClon->GetMinimum(),mixClon->GetMaximum(),1,3,2);
	PlotLine(0.148,0.148,mixClon->GetMinimum(),mixClon->GetMaximum(),1,3,2);
	
	x1Bin = D0_mixed->FindBin(x1);
	x2Bin = D0_mixed->FindBin(x2);
	nBins = x2Bin - x1Bin;
	



	// ******************************* Mixed Residual bg Subtraction *******************************
	fOutskirt->GetParameters(par);
	TF1 *fSignalMixed = new TF1("fSignalMixed",fitFunction,x1,x2,7);
	par[4]=0.304228; par[5]=0.145558; par[6]=0.000809652;
	fSignalMixed->SetLineColor(kRed);
	fSignalMixed->SetLineWidth(6);
	fSignalMixed->SetParameters(par);
	/////////fSignalMixed->SetParLimits(4,0,1000);
	fSignalMixed->SetParLimits(5,0.143,0.147);
	fSignalMixed->SetParLimits(6,0.0001,0.0015);
	fSignalMixed->SetNpx(1000);
//	fSignalMixed->FixParameter(4,0);
//	fSignalMixed->FixParameter(5,1.865);
//	fSignalMixed->FixParameter(6,0.0111);
	////fSignalMixed->FixParameter(0,par[0]);
	////fSignalMixed->FixParameter(1,par[1]);
	fSignalMixed->FixParameter(2,par[2]);
	fSignalMixed->FixParameter(3,par[3]);
	mixClon->Fit("fSignalMixed","","",x1,x2);
	//////////fSignalMixed->Draw("same");

	TF1 *fsignalonlymix = new TF1("fsignalonlymix",normalDistribution,x1,x2,3);
	fsignalonlymix->FixParameter(0,fSignalMixed->GetParameter(4));
	fsignalonlymix->FixParameter(1,fSignalMixed->GetParameter(5));
	fsignalonlymix->FixParameter(2,fSignalMixed->GetParameter(6));

	fsignalonlymix->SetLineColor(kBlack);
	fsignalonlymix->SetLineWidth(6);
	fsignalonlymix->Draw("same");

	TF1 *fresbgonlymix = new TF1("fresbgonlymix",parabola,x1,x2,4);
	fresbgonlymix->FixParameter(0,fSignalMixed->GetParameter(0));
	fresbgonlymix->FixParameter(1,fSignalMixed->GetParameter(1));
	fresbgonlymix->FixParameter(2,fSignalMixed->GetParameter(2));
	fresbgonlymix->FixParameter(3,fSignalMixed->GetParameter(3));

	fresbgonlymix->SetLineColor(kRed);
	fresbgonlymix->SetLineWidth(6);
	fresbgonlymix->SetLineStyle(6);
	//////////fresbgonlymix->Draw("same");

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

	x1aBin = RightSign->FindBin(x1a);
	x2aBin = RightSign->FindBin(x2a);

	cout << x1aBin << "     " << x2aBin << endl;

	x1aBinLowEdge = RightSign->GetBinLowEdge(x1aBin);
	x2aBinHighEdge = RightSign->GetBinLowEdge(x2aBin) + RightSign->GetBinWidth(x2aBin);
	cout << x1aBinLowEdge << "    "  << x2aBinHighEdge << endl;

	signalgauss = (fsignalonlymix->Integral(x1aBinLowEdge,x2aBinHighEdge))/mixClon->GetBinWidth(1);
	cout << signalgauss << "     " << sqrt(RightSign->Integral(x1aBin, x2aBin)) << endl;

	sg = signalgauss/sqrt(RightSign->Integral(x1aBin, x2aBin));
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
	
	WrongSignClon->DrawCopy("same");
	SideBandClon->DrawCopy("same");
	/////mixClon->DrawCopy();
	
	TLegend *resultsLegend = new TLegend(0.2,0.82,0.5,0.98,"","brNDC");
	resultsLegend->AddEntry(SideBandClon,"Side Band","p");
	resultsLegend->AddEntry(WrongSignClon,"Wrong Sign","p");
	//////////resultsLegend->AddEntry(mixClon,"mixed event","p");
	resultsLegend->Draw();
	infoSideBand->DrawClone();
	infoWrongSign->DrawClone();
	//////////infoMixed->DrawClone();


	
//	TCanvas *c30 = new TCanvas("c30");
//	c30->cd();
//	
//	D0_WrongSign->SetMarkerStyle(20);
//	Dstar_SideBand->SetMarkerStyle(20);
//	D0_mixed->SetMarkerStyle(20);
//	
//	D0_WrongSign->Draw();
//	Dstar_SideBand->Draw("same");
//	D0_mixed->Draw("same");
//	
//	TLegend *resultsLegend = new TLegend(0.2,0.82,0.5,0.98,"","brNDC");
//	resultsLegend->AddEntry(D0_WrongSign,"mixed event","p");
//	resultsLegend->AddEntry(Dstar_SideBand,"Side Band","p");
//	resultsLegend->AddEntry(D0_mixed,"Wrong sign","p");
//	resultsLegend->Draw();
	
	
	
	
	//********************************** Draw and Save Zooms ************************************
	
	TCanvas *c11 = new TCanvas("c11","c11",400,400);
	c11->SetBottomMargin(0.15);
	c11->SetTickx();
	c11->SetTicky();
	c11->SetLeftMargin(0.15);
	WrongSignClon->SetTitle("After WrongSign bg. subtracted");
	WrongSignClon->GetXaxis()->SetRangeUser(0.139,0.161);
	WrongSignClon->GetYaxis()->SetRangeUser(-210, 510);
	if (pT1 > 2.9 && pT2 < 4.3) {
		WrongSignClon->GetYaxis()->SetRangeUser(-80, 260);
	}
	if (pT1 > 4.1 && pT2 < 6.5) {
		WrongSignClon->GetYaxis()->SetRangeUser(-30, 110);
	}
	WrongSignClon->GetYaxis()->SetTitle(Form("Entries (/%.4f GeV/c^{2})", 0.0002*rebin));
	WrongSignClon->GetYaxis()->SetTitleOffset(1.4);
	WrongSignClon->GetYaxis()->SetLabelSize(0.05);
	WrongSignClon->GetXaxis()->SetLabelSize(0.05);
	WrongSignClon->GetYaxis()->SetTitleSize(0.05);
	WrongSignClon->GetXaxis()->SetTitleSize(0.05);
    WrongSignClon->GetYaxis()->SetLabelFont(42);
	WrongSignClon->GetXaxis()->SetLabelFont(42);
	WrongSignClon->GetYaxis()->SetTitleFont(42);
	WrongSignClon->GetXaxis()->SetTitleFont(42);
	//	WrongSignClon->SetMarkerStyle(8);
	WrongSignClon->SetMarkerColor(4);
	WrongSignClon->SetMarkerSize(1);
	WrongSignClon->SetLineWidth(2);
	WrongSignClon->Draw();
	infoWrongSign->SetPad(0.6,0.5,1.,0.9);
	infoWrongSign->SetFillColor(10);
	infoWrongSign->Draw("same");
	//////fsignalonlyWrongSign->Draw("same");
	//////fresbgonlyWrongSign->Draw("same");
	TLegend *legLike = new TLegend(0.1934673,0.6093333,0.4924623,0.8093333,NULL,"brNDC");
	legLike->AddEntry(WrongSignClon, "Data", "lep");
	legLike->AddEntry(fSignalWrongSign, "Fit", "l");
	legLike->Draw();

	TCanvas *c10 = new TCanvas("c10","c10",400,400);
	c10->SetBottomMargin(0.15);
	c10->SetLeftMargin(0.15);
	c10->SetTickx();
	c10->SetTicky();

	SideBandClon->SetTitle("After SideBand bg. subtracted");
	SideBandClon->GetXaxis()->SetRangeUser(0.139,0.161);
	SideBandClon->GetYaxis()->SetRangeUser(-210, 510);
	if (pT1 > 2.9 && pT2 <= 4.3) {
		SideBandClon->GetYaxis()->SetRangeUser(-80, 260);
	}
	if (pT1 > 4.1 && pT2 <= 6.5) {
		SideBandClon->GetYaxis()->SetRangeUser(-30, 110);
	}
	SideBandClon->GetYaxis()->SetTitle(Form("Entries (/%.4f GeV/c^{2})", 0.0002*rebin));
	SideBandClon->GetYaxis()->SetTitleOffset(1.4);
	SideBandClon->GetYaxis()->SetLabelSize(0.05);
	SideBandClon->GetXaxis()->SetLabelSize(0.05);
	SideBandClon->GetYaxis()->SetTitleSize(0.05);
	SideBandClon->GetXaxis()->SetTitleSize(0.05);
	SideBandClon->GetYaxis()->SetLabelFont(42);
	SideBandClon->GetXaxis()->SetLabelFont(42);
	SideBandClon->GetYaxis()->SetTitleFont(42);
	SideBandClon->GetXaxis()->SetTitleFont(42);
	SideBandClon->SetMarkerSize(1);
	SideBandClon->SetLineWidth(2);
	//	SideBandClon->SetMarkerStyle(8);
	SideBandClon->SetMarkerColor(kGreen+1);
	///SideBandClon->Scale(2.);
	SideBandClon->Draw();
	infoSideBand->SetPad(0.6,0.5,1.,0.9);
	infoSideBand->SetFillColor(10);
	infoSideBand->Draw("same");
	//////fsignalonlySideBand->Draw("same");
	//////fresbgonlySideBand->Draw("same");TLegend *legRot = new TLegend(0.1934673,0.6093333,0.4924623,0.8093333,NULL,"brNDC");
	TLegend *legRot = new TLegend(0.1934673,0.6093333,0.4924623,0.8093333,NULL,"brNDC");
	legRot->AddEntry(SideBandClon, "Data", "lep");
	legRot->AddEntry(fSignalSideBand, "Fit", "l");
	legRot->Draw();

	
	TCanvas *c12 = new TCanvas("c12","c12",400,400);
	c12->SetBottomMargin(0.15);
	c12->SetLeftMargin(0.15);
	c12->SetTickx();
	c12->SetTicky();
//	hSignal_Mixed->SetTitle("Mixed event");
//	hSignal_Mixed->GetXaxis()->SetRangeUser(0.139,0.161);
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
//	c12->cd(2);
//	gPad->SetBottomMargin(0.2);
//	gPad->SetLeftMargin(0.15);
	mixClon->SetTitle("After Mixed-Event bg. subtracted");
	mixClon->GetXaxis()->SetRangeUser(0.139,0.161);
	mixClon->GetYaxis()->SetRangeUser(-210, 510);
	mixClon->GetYaxis()->SetTitle(Form("Entries (/%.4f GeV/c^{2})", 0.0002*rebin));
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
	/*//////////mixClon->Draw();
	infoMixed->SetPad(0.6,0.5,1.,0.9);
	infoMixed->SetFillColor(10);
	infoMixed->Draw("same");
	fsignalonlymix->Draw("same");
	fresbgonlymix->Draw("same");//////////*/
//	
//	TF1 *signalFcn = new TF1("signalFcn",lorentzianPeak,0.8,1.05,3);
//	c4->cd();
//	signalFcn->SetParameters(1000,1,0.892);
//	D0_WrongSign->Fit("signalFcn","","",0.8,0.95);
	


	//********************************** Draw a more plots into one **************************************

	c1->cd();
	TH1D *CloneSB = (TH1D*)gDirectory->Get("CloneSB");
	TH1D *CloneWS = (TH1D*)gDirectory->Get("CloneWS");
	TH1D *CloneM = (TH1D*)gDirectory->Get("CloneM");
	CloneSB->Scale(2);
	CloneWS->Scale(2);
	CloneM->Scale(2);
	CloneSB->SetMarkerColor(8);
	CloneSB->SetLineColor(8);
	CloneSB->SetMarkerStyle(8);
	CloneWS->SetMarkerColor(4);
	CloneWS->SetLineColor(4);
	CloneWS->SetMarkerStyle(8);
	CloneM->SetMarkerColor(kCyan+3);
	CloneM->SetLineColor(kCyan+3);
	CloneM->SetMarkerStyle(8);

	CloneWS->Draw("same");
	CloneSB->Draw("same");
	//////////CloneM->Draw("same");

	TGaxis* RawYieldAxis = new TGaxis(2.5,0,2.5,30e4,0,15e4,510,"+L");
	
	
	//+ : draw on positive side
	//L : left adjusted
	RawYieldAxis->SetName("RawYieldAxis");
	RawYieldAxis->SetLineColor(8);
	RawYieldAxis->SetTextColor(4);
	// RawYieldAxis->SetTitle("Raw Yield (/0.1 GeV/c^{2})");
	RawYieldAxis->SetLabelColor(4);
	RawYieldAxis->SetLabelFont(42);
	RawYieldAxis->Draw();

    
	
	leg->AddEntry(CloneSB,"2*(Right Sign - Side Band)","p");
	leg->AddEntry(CloneWS,"2*(Right Sign - Wrong sign)","p");
	//////////leg->AddEntry(CloneM,"2*(RightSign Sign - Mixed Event)","p");
	
	

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

