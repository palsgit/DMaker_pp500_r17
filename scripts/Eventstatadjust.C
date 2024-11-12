#include <TFile.h>
#include <TList.h>
#include <TH1.h>

void Eventstatadjust() {
    // Open the .root file
    TFile *file = TFile::Open("anaoutput_dstar_forCollabmeeting.root");

    // Get the TList from the file
    ///////TList *list = (TList*)file->Get("picoD0AnaMaker");

    // Get the histogram from the list
    /////TH1 *hEventStat0 = (TH1*)file->Get("hEventStat0");
    TH1 *hEventStat1 = (TH1*)file->Get("hEventStat1");
    TH1 *hPrimVtZ = (TH1*)file->Get("hPrimVtZ");
    TH1 *hVzVPDvsVzTPCDiff = (TH1*)file->Get("hVzVPDvsVzTPCDiff");
    TH2 *hPrimVtXY = (TH2*)file->Get("hPrimVtXY");
    

    /*TH2 *hBetavsP_tr = (TH2*)file->FindObject("hBetavsP_tr");

    // Remove bin contents of hBetavsP_tr where x <= 0.2
    for (int i = 1; i <= hBetavsP_tr->GetNbinsX(); ++i) {
        if (hBetavsP_tr->GetXaxis()->GetBinCenter(i) <= 0.2) {
            hBetavsP_tr->SetBinContent(i, 0);
        }
    }

    hBetavsP_tr->GetXaxis()->SetRangeUser(0.2, 3.5);
    hBetavsP_tr->GetYaxis()->SetRangeUser(0.0, 3.0);
    hBetavsP_tr->SetTitle("1/#beta vs p after track quality cuts");

    

    hBetavsP_tr->Draw("colz");
    double y_min = hBetavsP_tr->GetYaxis()->GetXmin();
    double y_max = hBetavsP_tr->GetYaxis()->GetXmax();
    TLine *line = new TLine(1.6, y_min, 1.6, 3.0);
    line->SetLineColor(kRed);  // Set line color to red
    line->SetLineStyle(10);     // Set line style to dashed
    line->Draw("same");        // Draw line on the same canvas
    */

   TH1 *hEventStat1_final = new TH1F("hEventStat1_final", "Event Statistics", 6, 0.5, 6.5);
    for (int i = 1; i <= 3; i++) {
        hEventStat1_final->GetXaxis()->SetBinLabel(i, hEventStat1->GetXaxis()->GetBinLabel(i));
        hEventStat1_final->SetBinContent(i, hEventStat1->GetBinContent(i));
    }

    cout << "Bin 3:" << hEventStat1_final->GetBinContent(3) << endl;

    int bin1 = hPrimVtZ->GetXaxis()->FindBin(-60.0);
    int bin2 = hPrimVtZ->GetXaxis()->FindBin(60.0);

    long integral = hPrimVtZ->Integral(bin1, bin2);

    cout << "Integral: " << integral << endl;

    hEventStat1_final->SetBinContent(4, (integral < hEventStat1_final->GetBinContent(3)) ? integral : hEventStat1_final->GetBinContent(3));
    hEventStat1_final->GetXaxis()->SetBinLabel(4, "|V_{z}^{TPC}| < 60 cm");

    cout << "Bin 4:" << hEventStat1_final->GetBinContent(4) << endl;

    bin1 = hVzVPDvsVzTPCDiff->GetXaxis()->FindBin(-4.0);
    bin2 = hVzVPDvsVzTPCDiff->GetXaxis()->FindBin(4.0);

    integral = hVzVPDvsVzTPCDiff->Integral(bin1, bin2);

    cout << "Integral: " << integral << endl;

    hEventStat1_final->SetBinContent(5, (integral < hEventStat1_final->GetBinContent(4)) ? integral : hEventStat1_final->GetBinContent(4));
    hEventStat1_final->GetXaxis()->SetBinLabel(5, "|V_{z}^{TPC} - V_{z}^{VPD}| < 4 cm");

    cout << "Bin 5:" << hEventStat1_final->GetBinContent(5) << endl;

    bin1 = hPrimVtXY->GetXaxis()->FindBin(-0.3);
    bin2 = hPrimVtXY->GetXaxis()->FindBin(0.14);
    int bin3 = hPrimVtXY->GetYaxis()->FindBin(-0.26);
    int bin4 = hPrimVtXY->GetYaxis()->FindBin(0.02);

    integral = hPrimVtXY->Integral(bin1, bin2, bin3, bin4);

    cout << "Integral: " << integral << endl;

    hEventStat1_final->SetBinContent(6, (integral < hEventStat1_final->GetBinContent(5)) ? integral : hEventStat1_final->GetBinContent(5));
    hEventStat1_final->GetXaxis()->SetBinLabel(6, "V_{xy} cut");

    cout << "Bin 6:" << hEventStat1_final->GetBinContent(6) << endl;

    hEventStat1_final->SetBinContent(7, hEventStat1_final->GetBinContent(6));
    hEventStat1_final->GetXaxis()->SetBinLabel(7, "Accepted");

    cout << "Bin 7:" << hEventStat1_final->GetBinContent(7) << endl;



    // Create a new .root file to write the histograms
    TFile *outfile = TFile::Open("output_histograms.root", "RECREATE");

    
    hEventStat1_final->Write();
    //////hBetavsP_tr->Write();

    // Close the output file
    outfile->Close();
    

    // Don't forget to close the file
    file->Close();
}