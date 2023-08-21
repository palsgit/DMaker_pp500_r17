void DrawBinContents() {
    TFile *file1 = TFile::Open("output_all.root", "READ");
    TList *list1 = (TList*) file1->Get("picoD0AnaMaker");
    TH1D *hEventStat1 = (TH1D*) list1->FindObject("hEventStat1");

    
    TCanvas *canvas = new TCanvas("canvas", "hEventStat1ogram with Bin Contents", 800, 600);
    canvas->cd();
    
    hEventStat1->Draw("hist"); // Draw the hEventStat1ogram first
    
    // Loop over the bins and put bin contents as TLatex above each bin
    for (Int_t iBin = 1; iBin <= hEventStat1->GetNbinsX(); ++iBin) {
        Double_t binContent = hEventStat1->GetBinContent(iBin);
        Double_t binCenterX = hEventStat1->GetBinCenter(iBin);
        Double_t binCenterY = hEventStat1->GetBinContent(iBin);
        
        TLatex *latex = new TLatex(binCenterX, binCenterY, Form("%.2f", binContent));
        latex->SetTextAlign(22);
        latex->SetTextSize(0.03);
        latex->SetTextAngle(90);
        latex->Draw();
    }
    
    canvas->Update();
}
