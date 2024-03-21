#include <TFile.h>
#include <TList.h>
#include <TH1.h>

void Eventstatadjust() {
    // Open the .root file
    TFile *file = TFile::Open("output_all.root");

    // Get the TList from the file
    TList *list = (TList*)file->Get("picoD0AnaMaker");

    // Get the histogram from the list
    TH1 *hEventStat0 = (TH1*)list->FindObject("hEventStat0");
    TH1 *hEventStat1 = (TH1*)list->FindObject("hEventStat1");
    TH1 *hVzVPDvsVzTPCDiff = (TH1*)list->FindObject("hVzVPDvsVzTPCDiff");

    TH2 *hBetavsP_tr = (TH2*)list->FindObject("hBetavsP_tr");

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


    


    // Find the bin numbers corresponding to the range
    int bin1 = hVzVPDvsVzTPCDiff->GetXaxis()->FindBin(-4.0);
    int bin2 = hVzVPDvsVzTPCDiff->GetXaxis()->FindBin(4.0);

    // Take the integral of the histogram in the range
    int integral = hVzVPDvsVzTPCDiff->Integral(bin1, bin2);

    cout << "Integral: " << integral << endl;




    // Modify the binx = 5 of the histogram with a different bincontent
    hEventStat0->SetBinContent(5, integral);
    hEventStat0->GetXaxis()->SetBinLabel(6, "V_{xy}");
    hEventStat0->GetXaxis()->SetBinLabel(5, "|V_{z}^{TPC} - V_{z}^{VPD}| < 4 cm");
    hEventStat0->GetXaxis()->SetBinLabel(4, "|V_{z}^{TPC}| < 60 cm");
    hEventStat0->GetXaxis()->SetBinLabel(3, "VPDMB trigger");
    /////hEventStat0->Draw();


    int binContent = hEventStat1->GetBinContent(4);
    hEventStat1->SetBinContent(5, (integral < binContent) ? integral : binContent);
    hEventStat1->SetBinContent(6, (hEventStat1->GetBinContent(6) < hEventStat1->GetBinContent(5)) ? hEventStat1->GetBinContent(6) : hEventStat1->GetBinContent(5));
    hEventStat1->SetBinContent(7, (hEventStat1->GetBinContent(7) < hEventStat1->GetBinContent(6)) ? hEventStat1->GetBinContent(7) : hEventStat1->GetBinContent(6));
    
    //////TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    hEventStat1->GetXaxis()->SetBinLabel(6, "V_{xy}");
    hEventStat1->GetXaxis()->SetBinLabel(5, "|V_{z}^{TPC} - V_{z}^{VPD}| < 4 cm");
    hEventStat1->GetXaxis()->SetBinLabel(4, "|V_{z}^{TPC}| < 60 cm");
    hEventStat1->GetXaxis()->SetBinLabel(3, "VPDMB trigger");
    //////hEventStat1->Draw();

    // Create a new .root file to write the histograms
    TFile *outfile = TFile::Open("output_histograms.root", "RECREATE");

    // Write the histograms to the file
    hEventStat0->Write();
    hEventStat1->Write();
    hBetavsP_tr->Write();

    // Close the output file
    outfile->Close();
    

    // Don't forget to close the file
    file->Close();
}