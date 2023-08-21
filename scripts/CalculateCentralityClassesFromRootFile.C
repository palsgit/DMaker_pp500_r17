#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TFile.h>
#include <TList.h>
#include <iostream>


void CalculateCentralityClasses(TH1D *multiplicityHist) {
    if (!multiplicityHist) {
        std::cerr << "Error: Invalid histogram provided." << std::endl;
        return;
    }

    double totalEvents = multiplicityHist->Integral(0, multiplicityHist->GetNbinsX());

    const int nCentralityClasses = 10;  // Dividing into 10% intervals
    /////double centralityPercentiles[nCentralityClasses + 1] = {0.0, 10.0, 20.0, 40.0, 60.0, 80.0, 100.0};
    double centralityPercentiles[nCentralityClasses + 1] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0};


    TH1D *centralityHist = new TH1D("centralityHist", "Centrality Classes", nCentralityClasses, 0, nCentralityClasses);
    centralityHist->GetXaxis()->SetTitle("Centrality (%)");
    centralityHist->GetYaxis()->SetTitle("Events");

    double sum = 0.0;
    int centralityBin = 1;

    for (int bin = 0; bin <= multiplicityHist->GetNbinsX(); ++bin) {
        sum += multiplicityHist->GetBinContent(bin);
        double currentPercent = (sum / totalEvents) * 100.0;
        if (currentPercent == 100) break;
        cout << currentPercent << "     "  <<  multiplicityHist->GetBinCenter(bin) << endl;

        if (currentPercent > centralityPercentiles[centralityBin]) {

            cout << "Multiplicity  " << multiplicityHist->GetBinCenter(bin) << "  centrality class " << centralityBin << endl;
            centralityHist->SetBinContent(centralityBin, sum);
            /////sum = 0.0;
            ++centralityBin;
        }

    }

}

void CalculateCentralityClassesFromRootFile(const char *fileName = "output_all.root", const char *histogramName = "hPrimVtZ_evcut") {
    TFile *file = TFile::Open(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file." << std::endl;
        return;
    }

    TList *histList = (TList*) file->Get("picoD0AnaMaker"); // Replace "YourListName" with the actual list name
    if (!histList) {
        std::cerr << "Error: List not found in the file." << std::endl;
        file->Close();
        return;
    }

    TH1D *multiplicityHist = (TH1D*) histList->FindObject(histogramName);
    if (!multiplicityHist) {
        std::cerr << "Error: Histogram not found in the list." << std::endl;
        file->Close();
        return;
    }

    TCanvas *c1 = new TCanvas("c1", "Centrality Classes", 800, 600);
    multiplicityHist->Draw();

    CalculateCentralityClasses(multiplicityHist);

    // Rest of the code remains the same as before...
    // ... CalculateCentralityClasses function and so on ...

    file->Close();
}



