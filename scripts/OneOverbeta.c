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

void OneOverbeta()
{
    TFile *file1 = TFile::Open("output_all.root", "READ");
    TList *list1 = (TList*) file1->Get("picoD0AnaMaker");
    TH2D *h_OneOverBetaDiffK_TPC = (TH2D*) list1->FindObject("h_OneOverBetaDiffK_TPC");
    h_OneOverBetaDiffK_TPC->SetTitle("nSigmaTOF_Kaon; p [GeV]; n#sigma_{K}^{TOF}; Counts");
    TH2D *h_OneOverBetaDiffPi_TPC = (TH2D*) list1->FindObject("h_OneOverBetaDiffPi_TPC");
    h_OneOverBetaDiffPi_TPC->SetTitle("nSigmaTOF_Pion; p [GeV]; n#sigma_{#pi}^{TOF}; Counts");

    h_OneOverBetaDiffK_TPC->GetXaxis()->SetRangeUser(0.2, 3.5);

    int nSlices = 60;
float sliceWidth = 0.02;

float fitrange[nSlices][2];

for (Int_t i = 0; i < 12; i++)
{
    fitrange[i][0] = -8;
    fitrange[i][1] = 8;
}

for (Int_t i = 12; i < 15; i++)
{
    fitrange[i][0] = -6;
    fitrange[i][1] = 6;
}

for (Int_t i = 15; i < 17; i++)
{
    fitrange[i][0] = -4;
    fitrange[i][1] = 4;
}

for (Int_t i = 17; i < 19; i++)
{
    fitrange[i][0] = -3;
    fitrange[i][1] = 3;
}

for (Int_t i = 19; i < 50; i++)
{
    fitrange[i][0] = -2.5;
    fitrange[i][1] = 2.5;
}

for (Int_t i = 50; i < 60; i++)
{
    fitrange[i][0] = -2.0;
    fitrange[i][1] = 2.0;
}


// Create a canvas to plot the slices
TFile* dataRes = new TFile("Projections.root","RECREATE");
TCanvas *canvas = new TCanvas("canvas", "Fitted Slices", 1200, 800);
canvas->Divide(5, 10);

TH2D *meanKaon = new TH2D("meankaon", "meankaon; p [GeV]; mean_{n#sigma^{TOF}; Entries", 100, 0.0, 2.0, 350, -0.5, 3.0); 
TH2D *sigmaKaon = new TH2D("sigmakaon", "sigmakaon; p [GeV]; #sigma_{n#sigma^{TOF}; Entries", 100, 0.0, 2.0, 500, 0.0, 5.0);
meanKaon->Sumw2();
sigmaKaon->Sumw2();
/*Double_t meanKaonarr[50];
Double_t sigmaKaonarr[50];
Double_t pKaonarr[50];
Double_t pKaonerrarr[50];
*/
vector<double> meanKaonvect;
vector<double> sigmaKaonvect;
vector<double> pKaonvect;
vector<double> pKaonerrvect;

// Fit each slice with a Gaussian and plot the result
for (int i = 0; i < nSlices; i++) {
    float sliceMin = i * sliceWidth;
    float sliceMax = (i + 1) * sliceWidth;

    // Create a new pad to plot the slice
    canvas->cd(i+1);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);

    // Project the 2D histogram onto the y-axis within the slice range
    TH1D *slice = h_OneOverBetaDiffK_TPC->ProjectionY(Form("slice_%d", i), h_OneOverBetaDiffK_TPC->GetXaxis()->FindBin(sliceMin), h_OneOverBetaDiffK_TPC->GetXaxis()->FindBin(sliceMax));
    slice->SetMarkerStyle(27);
    slice->SetMarkerSize(1);

    // Fit the slice with a Gaussian
    TF1 *fitFunc = new TF1(Form("fitFunc_%d", i), "gaus", sliceMin, sliceMax);
    
    
    if (sliceMin < 0.18) continue;
    slice->Fit(fitFunc, "Q", "", fitrange[i][0], fitrange[i][1]);

    // Get the fit parameters
    float fitMean = fitFunc->GetParameter(1);
    float fitSigma = fitFunc->GetParameter(2);
    float fitChi2 = fitFunc->GetChisquare();
    float NDF = fitFunc->GetNDF();
    
    ////meanKaonvect.pushback
    meanKaon->Fill((sliceMin+sliceMax)/2.0, fitMean);
    sigmaKaon->Fill((sliceMin+sliceMax)/2.0, fitSigma);
    printf("Slice %d: mean = %.2f, sigma = %.2f, chi2/NDF = %.2f\n", i, fitMean, fitSigma, fitChi2/NDF);

    // Plot the slice and the fit function
    slice->Draw("Hist");
    fitFunc->Draw("same");

    // Set the axis labels and titles
    slice->GetXaxis()->SetTitle("Y");
    slice->GetYaxis()->SetTitle("Counts");
    slice->SetTitle(Form("Slice %d: Mean=%.2f, Sigma=%.2f, Chi2/NDF=%.2f", i, fitMean, fitSigma, fitChi2/NDF));

    // Delete the slice and fit function objects
    slice->Write();
    fitFunc->Write();
    delete slice;
    delete fitFunc;
}

TF1 *fitmean = new TF1("fitmean",fitf,0.2,1.,4);
fitmean->SetParameter(0,0.03);
fitmean->SetParameter(1,0.0014);
fitmean->SetParameter(2,0.1);
fitmean->SetParameter(3,6.9);
TProfile *prof_meanKaon;
///prof_meanKaon->Sumw2();
prof_meanKaon = meanKaon->ProfileX("prof_meanKaon", 1, -1);
prof_meanKaon->Fit(fitmean,"L","LN",0.2,1.);

 double r0 = fitmean->GetParameter(0);
 double r1 = fitmean->GetParameter(1);
 double r2 = fitmean->GetParameter(2);
 double r3 = fitmean->GetParameter(3);

 cout << "Parametr [0] (mean kaon) je " << r0 <<endl;
 cout << "Parametr [1] (mean kaon) je " << r1 <<endl;
 cout << "Parametr [2] (mean kaon) je " << r2 <<endl;
 cout << "Parametr [3] (mean kaon) je " << r3 <<endl;

prof_meanKaon->Write();

TF1 *fitsigma = new TF1("fitsigma",fitf,0.2,1.,4);
fitsigma->SetParameter(0,0.9);
fitsigma->SetParameter(1,0.02);
fitsigma->SetParameter(2,0.08);
fitsigma->SetParameter(3,4.23);
TProfile *prof_sigmaKaon;
////prof_sigmaKaon->Sumw2();
prof_sigmaKaon = sigmaKaon->ProfileX("prof_sigmaKaon", 1, -1);
prof_sigmaKaon->Fit(fitsigma,"L","LN",0.2,1.);

double p0 = fitsigma->GetParameter(0);
double p1 = fitsigma->GetParameter(1);
double p2 = fitsigma->GetParameter(2);
double p3 = fitsigma->GetParameter(3);

cout << "Parametr [0] (sigma kaon) je " << p0 <<endl;
cout << "Parametr [1] (sigma kaon) je " << p1 <<endl;
cout << "Parametr [2] (sigma kaon) je " << p2 <<endl;
cout << "Parametr [3] (sigma kaon) je " << p3 <<endl;

prof_sigmaKaon->Write();
// Save the canvas as a PDF file
canvas->SaveAs("fitted_slices.pdf");

// Delete the canvas object

  ////TFile* dataRes = new TFile("Projections.root","RECREATE");
     
    h_OneOverBetaDiffK_TPC->Write();
    h_OneOverBetaDiffPi_TPC->Write();
    dataRes->Close();
    file1->Close();

    delete canvas;
    
}

