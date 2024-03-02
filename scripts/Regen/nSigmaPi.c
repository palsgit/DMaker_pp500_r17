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
Double_t nSigmaTOFPionTrack_upper(Double_t x) { 
	double f_nsigmaTOFPionTrack_res =  0.907 + (0.002/(pow((x - 0.164), 1.990)));
    double f_nsigmaTOFPionTrack_pos =  0.010 + (0.005/(pow((x - 0.080), 2.544))); 

    /*float mSigma = 3.18;
    if (x > 2.0) mSigma = 2.41;*/
    /////float mSigma = 2.5; //weighted average
    float mSigma = 2.74;//midstep
	double Pion_higher = mSigma*f_nsigmaTOFPionTrack_res + f_nsigmaTOFPionTrack_pos;
    return Pion_higher;  
	}
Double_t nSigmaTOFPionTrack_lower(Double_t x) {
	double f_nsigmaTOFPionTrack_res =  0.907 + (0.002/(pow((x - 0.164), 1.990)));
    double f_nsigmaTOFPionTrack_pos =  0.010 + (0.005/(pow((x - 0.080), 2.544))); 
    /*float mSigma = -2.31;
    if (x > 1.25 && x < 1.50) mSigma = -1.56;
    if (x > 1.50 && x < 1.65) mSigma = -1.18;
    if (x > 1.65) mSigma = -0.8;*/
    float mSigma = -4.42; //midstep
    /////float mSigma = -1.0; //weighted average
	 double Pion_lower = mSigma*f_nsigmaTOFPionTrack_res + f_nsigmaTOFPionTrack_pos;
	 return Pion_lower; 
	 }

void nSigmaPi()
{
    TFile *file1 = TFile::Open("output_all.root", "READ");
    TList *list1 = (TList*) file1->Get("picoD0AnaMaker");
    TH2D *h_nSigmaOneOverBetaK_tr = (TH2D*) list1->FindObject("h_nSigmaOneOverBetaK_tr");
    h_nSigmaOneOverBetaK_tr->SetTitle("n#sigma_{K}^{TOF}; p [GeV]; n#sigma_{K}^{TOF}; Entries");
    TH2D *h_nSigmaOneOverBetaPi_tr = (TH2D*) list1->FindObject("h_nSigmaOneOverBetaPi_tr");
    h_nSigmaOneOverBetaPi_tr->SetTitle("n#sigma_{#pi}^{TOF}; p [GeV]; n#sigma_{#pi}^{TOF}; Entries");

    TH2D *h_nSigmadEdxK_tr = (TH2D*) list1->FindObject("h_nSigmadEdxK_tr");
    h_nSigmadEdxK_tr->SetTitle("n#sigma_{K}^{TPC}; p [GeV]; n#sigma_{K}^{TPC}; Entries");
    TH2D *h_nSigmadEdxPi_tr = (TH2D*) list1->FindObject("h_nSigmadEdxPi_tr");
    h_nSigmadEdxPi_tr->SetTitle("n#sigma_{#pi}^{TPC}; p [GeV]; n#sigma_{#pi}^{TPC}; Entries");

    TH2D *hDedx_tr = (TH2D*) list1->FindObject("hDedx_tr");
    hDedx_tr->SetTitle("dE/dx after track cut; p*charge [GeV/c]; dE/dx [keV/cm]; Entries");

    TH2D *hBetavsP_tr = (TH2D*) list1->FindObject("hBetavsP_tr");
    hBetavsP_tr->SetTitle("1/#beta after track cut; p [GeV/c]; 1/#beta; Entries");

    TH2D *h_nSigmaOneOverBetaK_TOF_TPC = (TH2D*) list1->FindObject("h_nSigmaOneOverBetaK_TOF_TPC");
    TH2D *h_nSigmaOneOverBetaPi_TOF_TPC = (TH2D*) list1->FindObject("h_nSigmaOneOverBetaPi_TOF_TPC");

    ////h_nSigmaOneOverBetaK_tr->GetXaxis()->SetRangeUser(0.18, 3.5);
    ////h_nSigmaOneOverBetaPi_tr->GetXaxis()->SetRangeUser(0.18, 3.5);

    int nSlices = 61;
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

for (Int_t i = 50; i < nSlices; i++)
{
    fitrange[i][0] = -2.0;
    fitrange[i][1] = 2.0;
}

/*for (Int_t i = 75; i < nSlices; i++)
{
    fitrange[i][0] = -1.0;
    fitrange[i][1] = 1.0;
}*/


// Create a canvas to plot the slices
TFile* dataRes = new TFile("Projections_pi.root","RECREATE");
TCanvas *canvas = new TCanvas("canvas", "Fitted Slices", 3600, 3600);
canvas->Divide(5, 10);

TH2D *meanPion = new TH2D("meanPion", "meanPion; p [GeV]; mean_{n#sigma^{TOF}}; Entries", 350, 0.0, 3.5, 400, -0.5, 3.5); 
TH2D *sigmaPion = new TH2D("sigmaPion", "sigmaPion; p [GeV]; #sigma_{n#sigma^{TOF}}; Entries", 350, 0.0, 3.5, 600, 0.0, 6.0);
meanPion->Sumw2();
sigmaPion->Sumw2();
/*Double_t meanPionarr[50];
Double_t sigmaPionarr[50];
Double_t pPionarr[50];
Double_t pPionerrarr[50];
*/
vector<double> meanPionvect;
vector<double> sigmaPionvect;
vector<double> pPionvect;
vector<double> pPionerrvect;
int firstslicegpad = 11;
// Fit each slice with a Gaussian and plot the result
for (int i = 0; i < nSlices; i++) {
    float sliceMin = i * sliceWidth;
    float sliceMax = (i + 1) * sliceWidth;

    if (sliceMin < 0.20) continue;

    // Create a new pad to plot the slice
    
    /*gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.15);
    gPad->SetTopMargin(0.05);*/

    // Project the 2D histogram onto the y-axis within the slice range
    TH1D *slice = h_nSigmaOneOverBetaPi_TOF_TPC->ProjectionY(Form("slice_%d", i), h_nSigmaOneOverBetaPi_TOF_TPC->GetXaxis()->FindBin(sliceMin), h_nSigmaOneOverBetaPi_TOF_TPC->GetXaxis()->FindBin(sliceMax));
    slice->SetMarkerStyle(27);
    slice->SetMarkerSize(1);

    // Fit the slice with a Gaussian
    TF1 *fitFunc = new TF1(Form("fitFunc_%d", i), "gaus", sliceMin, sliceMax);
    
    canvas->cd(i+1-firstslicegpad);
// Set the size dimensions of the gPad
    gPad->SetCanvasSize(3600, 3600);
    slice->Fit(fitFunc, "", "", fitrange[i][0], fitrange[i][1]);

    // Get the fit parameters
    float fitMean = fitFunc->GetParameter(1);
    float fitSigma = fitFunc->GetParameter(2);
    float fitChi2 = fitFunc->GetChisquare();
    float NDF = fitFunc->GetNDF();
    
    ////meanPionvect.pushback
    meanPion->Fill((sliceMin+sliceMax)/2.0, fitMean);
    sigmaPion->Fill((sliceMin+sliceMax)/2.0, fitSigma);
    printf("Slice %d: mean = %.2f, sigma = %.2f, chi2/NDF = %.2f\n", i, fitMean, fitSigma, fitChi2/NDF);

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
TF1 *fitmean = new TF1("fitmean",fitf,0.2, 5.0,4);
fitmean->SetParameter(0,-0.0427681);
fitmean->SetParameter(1,0.0416044);
fitmean->SetParameter(2, -0.201444);
fitmean->SetParameter(3,1.21657);

fitmean->SetLineColor(kBlack);
///TProfile *prof_meanPion;
///prof_meanPion->Sumw2();
////prof_meanPion = meanPion->ProfileX("prof_meanPion", 1, -1);
meanPion->Fit(fitmean,"N","",0.20,1.2);

 double r0 = fitmean->GetParameter(0);
 double r1 = fitmean->GetParameter(1);
 double r2 = fitmean->GetParameter(2);
 double r3 = fitmean->GetParameter(3);

 cout << "Parametr [0] (mean Pion) je " << r0 <<endl;
 cout << "Parametr [1] (mean Pion) je " << r1 <<endl;
 cout << "Parametr [2] (mean Pion) je " << r2 <<endl;
 cout << "Parametr [3] (mean Pion) je " << r3 <<endl;

////meanPion->GetYaxis()->SetRangeUser(-0.5, 4000.5);
meanPion->SetMarkerStyle(4);
meanPion->SetMarkerColor(kRed);
meanPion->Draw();
fitmean->Draw("same");
///prof_mean->Draw("same");
meanPion->Write();
//prof_mean->Write();

TCanvas *sigmacanvas = new TCanvas("sigmacanvas", "", 1200, 800);
sigmacanvas->cd();
TF1 *fitsigma = new TF1("fitsigma",fitf,0.2, 5.0,4);
fitsigma->SetParameter(0,0.9);
fitsigma->SetParameter(1,0.02);
fitsigma->SetParameter(2,0.08);
fitsigma->SetParameter(3,4.23);

fitsigma->SetLineColor(kBlack);
////TProfile *prof_sigma;
////prof_sigma->Sumw2();
////prof_sigma = sigma->ProfileX("prof_sigma", 1, -1);
////prof_sigma->Fit(fitsigma,"L","LN",0.18,1.2);
sigmaPion->Fit(fitsigma,"N","",0.20,1.2);

double p0 = fitsigma->GetParameter(0);
double p1 = fitsigma->GetParameter(1);
double p2 = fitsigma->GetParameter(2);
double p3 = fitsigma->GetParameter(3);

cout << "Parametr [0] (sigma Pion) je " << p0 <<endl;
cout << "Parametr [1] (sigma Pion) je " << p1 <<endl;
cout << "Parametr [2] (sigma Pion) je " << p2 <<endl;
cout << "Parametr [3] (sigma Pion) je " << p3 <<endl;

sigmaPion->SetMarkerStyle(4);
sigmaPion->SetMarkerColor(kRed);
sigmaPion->Draw();
fitsigma->Draw("same");
sigmaPion->Write();
// Save the canvas as a PDF file
//canvas->SaveAs("fitted_slices.pdf");

TCanvas *canvas2 = new TCanvas("cut_check_canvas", "check cut", 3600, 3600);
canvas2->SetLogz();
canvas2->cd();
TF1 *fitsigmaPionTrackupper = new TF1("fitsigmaPionTrackupper","nSigmaTOFPionTrack_upper(x)",0.2, 5.25);
TF1 *fitsigmaPionTracklower = new TF1("fitsigmaPionTracklower","nSigmaTOFPionTrack_lower(x)",0.2, 5.25);

h_nSigmaOneOverBetaPi_TOF_TPC->GetXaxis()->SetRangeUser(0.20, 3.5);
h_nSigmaOneOverBetaPi_TOF_TPC->GetYaxis()->SetRangeUser(-5.0, 5.0);
h_nSigmaOneOverBetaPi_TOF_TPC->Draw("colz");
fitsigmaPionTrackupper->Draw("same");
fitsigmaPionTracklower->Draw("same");



// Delete the canvas object

  ////TFile* dataRes = new TFile("Projections.root","RECREATE");
     
    h_nSigmaOneOverBetaK_tr->Write();
    h_nSigmaOneOverBetaPi_tr->Write();
    h_nSigmadEdxK_tr->Write();
    h_nSigmadEdxPi_tr->Write();
    h_nSigmaOneOverBetaK_TOF_TPC->Write();
    h_nSigmaOneOverBetaPi_TOF_TPC->Write();


    hDedx_tr->Write();
    hBetavsP_tr->Write();

    /*dataRes->Close();
    file1->Close();*/

    ////delete canvas;
    ////delete sigmacanvas;
    ///delete meancanvas;
    
}

