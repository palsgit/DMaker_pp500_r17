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
#include "TNtuple.h"

using namespace std;

Double_t fitf(Double_t *x, Double_t *par)
{ 
    return par[0]+(par[1]/(pow((x[0]+par[2]),par[3])));   
}

Double_t normalDistribution(Double_t *x, Double_t *par) {
	return par[0]/(TMath::Sqrt(2*TMath::Pi()*par[2]*par[2]))*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2*par[2]*par[2]));
	
}

void nSigmaslicesK(TString h_name)
{

    TH2D *h_nSigma = (TH2D*) picoD0AnaMaker->FindObject(h_name);
    TString title = h_name + "; p [GeV]; " + h_name + "; Entries";
    h_nSigma->SetTitle(title);
    
    int nSlices = 101;
////float sliceWidth = 0.02;

float fitrange[nSlices][2];

for (Int_t i = 0; i < 6; i++)
{
    fitrange[i][0] = -8;
    fitrange[i][1] = 8;
}

for (Int_t i = 6; i < 8; i++)
{
    fitrange[i][0] = -6;
    fitrange[i][1] = 6;
}

for (Int_t i = 8; i < 9; i++)
{
    fitrange[i][0] = -4;
    fitrange[i][1] = 4;
}

for (Int_t i = 9; i < 10; i++)
{
    fitrange[i][0] = -3;
    fitrange[i][1] = 3;
}

for (Int_t i = 10; i < 25; i++)
{
    fitrange[i][0] = -2.5;
    fitrange[i][1] = 2.5;
}

for (Int_t i = 25; i < 33; i++)
{
    fitrange[i][0] = -2.0;
    fitrange[i][1] = 2.0;
}

for (Int_t i = 33; i < 37; i++)
{
    fitrange[i][0] = -1.5;
    fitrange[i][1] = 1.5;
}

for (Int_t i = 37; i < nSlices; i++)
{
    fitrange[i][0] = -1.0;
    fitrange[i][1] = 1.0;
}


// Create a canvas to plot the slices
TFile* dataRes = new TFile("slice_Projections"+ h_name +".root","RECREATE");
TCanvas *canvas = new TCanvas("canvas", "Fitted Slices", 1200, 800);
canvas->Divide(5, 20);

/*TH2D *mean = new TH2D("mean", "mean; p [GeV]; #mu_{n#sigma}; Entries", 350, 0.0, 3.5, 400, -0.5, 3.5); 
TH2D *sigma = new TH2D("sigma", "sigma; p [GeV]; #sigma_{n#sigma}; Entries", 350, 0.0, 3.5, 600, 0.0, 6.0);*/


/*vector<double> meanvect;
vector<double> meanerrvect;
vector<double> sigmavect;
vector<double> sigmaerrvect;
vector<double> pvect;
vector<double> perrvect;*/

TGraphErrors *mean = new TGraphErrors();
mean->SetName("mean");
mean->SetTitle("mean; p [GeV]; #mu_{n#sigma}");
mean->SetMarkerStyle(20);
mean->SetMarkerSize(1);
mean->SetMarkerColor(kRed);
mean->SetLineColor(kRed);
mean->SetLineWidth(2);

TGraphErrors *sigma = new TGraphErrors();
sigma->SetName("sigma");
sigma->SetTitle("sigma; p [GeV]; #sigma_{n#sigma}");
sigma->SetMarkerStyle(20);
sigma->SetMarkerSize(1);
sigma->SetMarkerColor(kRed);
sigma->SetLineColor(kRed);
sigma->SetLineWidth(2);


int firstslicegpad = 11;

// Fit each slice with a Gaussian and plot the result
int numBinsX = h_nSigma->GetNbinsX();
float sliceWidth = h_nSigma->GetXaxis()->GetBinWidth(1);
for (int i = 0; i < nSlices; i++) {
    float sliceMin = h_nSigma->GetXaxis()->GetBinLowEdge(i);
    float sliceMax = sliceMin + sliceWidth;

    if (sliceMin < 0.16) continue;

    // Create a new pad to plot the slice
    
    /*gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);*/

    // Project the 2D histogram onto the y-axis within the slice range
    TH1D *slice = h_nSigma->ProjectionY(Form("slice_%d", i), i, i);
    slice->SetMarkerStyle(27);
    slice->SetMarkerSize(1);


    cout << "Slice " << i << ": " << i << " < p < " << i+1 << endl;
    cout << "Slice " << i << ": " << sliceMin << " < p < " << sliceMax << endl;
    

    // Fit the slice with a Gaussian
    TF1 *fitFunc = new TF1(Form("fitFunc_%d", i), "gaus", sliceMin, sliceMax);
    
    ////canvas->cd(i+1-firstslicegpad);
    canvas->cd(i);
// Set the size dimensions of the gPad
    gPad->SetCanvasSize(1200,4800);
    slice->Fit(fitFunc, "ME", "", fitrange[i][0], fitrange[i][1]);

    // Get the fit parameters
    float fitMean = fitFunc->GetParameter(1);
    float fitMeanErr = fitFunc->GetParError(1);
    float fitSigma = fitFunc->GetParameter(2);
    float fitSigmaErr = fitFunc->GetParError(2);
    float fitChi2 = fitFunc->GetChisquare();
    float NDF = fitFunc->GetNDF();
    
    ////meanvect.pushback
    /*mean->Fill((sliceMin+sliceMax)/2.0, fitMean);
    sigma->Fill((sliceMin+sliceMax)/2.0, fitSigma);*/
    mean->SetPoint(i, (sliceMin+sliceMax)/2.0, fitMean);
    mean->SetPointError(i, 0.0, fitMeanErr);
    sigma->SetPoint(i, (sliceMin+sliceMax)/2.0, fitSigma);
    sigma->SetPointError(i, 0.0, fitSigmaErr);
    printf("Slice %d: mean = %.5f, sigma = %.5f, p = %.5f GeV,  chi2/NDF = %.5f\n", i, fitMean, fitSigma, (sliceMin+sliceMax)/2.0, fitChi2/NDF);

    // Plot the slice and the fit function
    slice->GetXaxis()->SetRangeUser(-10, 10);
    

    // Set the axis labels and titles
    slice->GetXaxis()->SetTitle("n#sigma");
    slice->GetYaxis()->SetTitle("Entries");
    slice->SetTitle(Form("Slice %d: Mean=%.5f, Sigma=%.5f, p = %.5f GeV, Chi2/NDF=%.5f", i, fitMean, fitSigma, (sliceMin+sliceMax)/2.0, fitChi2/NDF));

    slice->Draw();
    fitFunc->Draw("same");

    // Delete the slice and fit function objects
    slice->Write();
    fitFunc->Write();
    ///delete slice;
    ////delete fitFunc;
}

TCanvas *meancanvas = new TCanvas("meancanvas", "mean", 1200, 800);
meancanvas->cd();
TF1 *fitmean = new TF1("fitmean",fitf,0.20,2.0,4);
fitmean->SetParameter(0,0.0188267);
fitmean->SetParameter(1,0.00558778);
fitmean->SetParameter(2,-0.125211);
fitmean->SetParameter(3,2.78366);

fitmean->SetLineColor(kBlack);
///TProfile *prof_mean;
///prof_mean->Sumw2();
///prof_mean = mean->ProfileX("prof_mean", 1, -1);
//prof_mean->Fit(fitmean,"L","LN",0.20,1.);
mean->Fit(fitmean,"MEN","",0.2, 2.0);


 double r0 = fitmean->GetParameter(0);
 double r1 = fitmean->GetParameter(1);
 double r2 = fitmean->GetParameter(2);
 double r3 = fitmean->GetParameter(3);

 cout << "Parametr [0] (mean ) je " << r0 <<endl;
 cout << "Parametr [1] (mean ) je " << r1 <<endl;
 cout << "Parametr [2] (mean ) je " << r2 <<endl;
 cout << "Parametr [3] (mean ) je " << r3 <<endl;

 cout << "chi2/ndf = " << fitmean->GetChisquare() << "/   " << fitmean->GetNDF() << endl;

/////mean->GetYaxis()->SetRangeUser(-0.5, 4000.5);
mean->SetMarkerStyle(4);
mean->SetMarkerColor(kRed);
mean->Draw("ap");
fitmean->Draw("same");
///prof_mean->Draw("same");
mean->Write();
//prof_mean->Write();

TCanvas *sigmacanvas = new TCanvas("sigmacanvas", "", 1200, 800);
sigmacanvas->cd();
TF1 *fitsigma = new TF1("fitsigma",fitf,0.20,2.0,4);
fitsigma->SetParameter(0,1.32732);
fitsigma->SetParameter(1,0.0257344);
fitsigma->SetParameter(2,-0.0213994);
fitsigma->SetParameter(3, 2.99865);

fitsigma->SetLineColor(kBlack);
////TProfile *prof_sigma;
////prof_sigma->Sumw2();
////prof_sigma = sigma->ProfileX("prof_sigma", 1, -1);
////prof_sigma->Fit(fitsigma,"L","LN",0.20,1.);
sigma->Fit(fitsigma,"MEN","",0.2, 2.0);

double p0 = fitsigma->GetParameter(0);
double p1 = fitsigma->GetParameter(1);
double p2 = fitsigma->GetParameter(2);
double p3 = fitsigma->GetParameter(3);

cout << "Parametr [0] (sigma ) je " << p0 <<endl;
cout << "Parametr [1] (sigma ) je " << p1 <<endl;
cout << "Parametr [2] (sigma ) je " << p2 <<endl;
cout << "Parametr [3] (sigma ) je " << p3 <<endl;

cout << "chi2/ndf = " << fitsigma->GetChisquare() << "/   " << fitsigma->GetNDF() << endl;
sigma->SetMarkerStyle(4);
sigma->SetMarkerColor(kRed);
sigma->Draw("ap");

fitsigma->Draw("same");
sigma->Write();
////prof_sigma->Write();
// Save the canvas as a PDF file
///canvas->SaveAs("fitted_slices.pdf");

// Delete the canvas object

  ////TFile* dataRes = new TFile("Projections.root","RECREATE");
     
    h_nSigma->Write();


    ////dataRes->Close();
    ///file1->Close();

    ////delete canvas;
    ////delete sigmacanvas;
    ////delete meancanvas;
    
}

