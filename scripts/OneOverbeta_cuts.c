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


void OneOverbeta_cuts()
{

    TFile *file1 = TFile::Open("Projections.root", "read");
    TH2D *h_OneOverBetaDiffK_TPC = (TH2D*)file1->Get("h_OneOverBetaDiffK_TPC");
    h_OneOverBetaDiffK_TPC->Rebin2D(5, 1);
    //h_OneOverBetaDiffK_TPC->Rebin2D(2, 2);
    h_OneOverBetaDiffK_TPC->FitSlicesY();
    /*TH2F *sigmakaon = (TH2F*)file1->Get("h_OneOverBetaDiffK_TPC_2");  ///sigma

    TF1 *fitsigma = new TF1("fitsigma","[0]+([1]/(x+[2])^[3])",0.22,0.85);
    fitsigma->SetParameter(0,7.84921e-02);// 1.12043e+00
    fitsigma->SetParameter(1,4.22205e-04);// 1.05315e-01
    fitsigma->SetParameter(2,6.17544e-02);// -1.41250e-01
    fitsigma->SetParameter(3,6.95657e+00);// 1.37384e+00

    sigmakaon->Fit(fitsigma,"L","LN",0.22,0.85);

    double r0 = fitsigma->GetParameter(0);
    double r1 = fitsigma->GetParameter(1);
    double r2 = fitsigma->GetParameter(2);
    double r3 = fitsigma->GetParameter(3);

    fitsigma->SetLineWidth(0.1);
    fitsigma->SetLineColor(3);
    */

     TF1 *fitsigmaupper = new TF1("fitsigmaupper","2.0*(1.09889 +(0.0115498 /(x + 0.0612253)^4.39129)) + (0.0345587 +(0.513216 /(x + 0.682501)^14.3189))",0.22,3.5);
     TF1 *fitsigmalower = new TF1("fitsigmalower","-1.15*(1.09889 +(0.0115498 /(x + 0.0612253)^4.39129)) + (0.0345587 +(0.513216 /(x + 0.682501)^14.3189))",0.22,3.5);
    ////TF1 *fitsigmaupper = new TF1("fitsigmaupper","1.55*(1.12190e+00 +(1.04744e-01 /(x - 1.41968e-01)^1.37291e+00)) + (7.84921e-02 +(4.22205e-04 /(x + 6.17544e-02)^6.95657e+00))",0.22,3.5);
    ////TF1 *fitsigmalower = new TF1("fitsigmalower","-1.15*(1.12190e+00 +(1.04744e-01 /(x - 1.41968e-01)^1.37291e+00)) + (7.84921e-02 +(4.22205e-04 /(x + 6.17544e-02)^6.95657e+00))",0.22,3.5);
    h_OneOverBetaDiffK_TPC->Draw("colz");
    //////sigmakaon->Draw("Psame");
    ///////fitsigma->Draw("same");
    fitsigmaupper->Draw("same");
    fitsigmalower->Draw("same");


    ///cout << fitsigma->GetChisquare() << "   "  << fitsigma->GetNDF() << endl;
    //h_OneOverBetaDiffK_TPC_1->Draw();
}
