#include "style.C"
void DrawBinContents() {
    style();
    TFile *file1 = TFile::Open("anaoutput_dstar_forCollabmeeting.root", "READ");
    TFile *file2 = TFile::Open("output_histograms.root", "READ");

    //////TList *list1 = (TList*) file1->Get("picoD0AnaMaker");
    /*TH1D *hEventStat1 = (TH1D*) file1->Get("hEventStat1");
    hEventStat1->SetStats(0);*/

    TH1D *hEventStat1 = (TH1D*) file2->Get("hEventStat1_final");

    
    TCanvas *canvas = new TCanvas("canvas", "hEventStat1ogram with Bin Contents", 800, 600);
    canvas->cd();
    
    hEventStat1->Draw("hist"); // Draw the hEventStat1ogram first
    
    // Loop over the bins and put bin contents as TLatex above each bin
    for (Int_t iBin = 1; iBin <= hEventStat1->GetNbinsX(); ++iBin) {
        long binContent = hEventStat1->GetBinContent(iBin);
        int binCenterX = hEventStat1->GetBinCenter(iBin);
        int binCenterY = hEventStat1->GetBinContent(iBin);
        
        TLatex *latex = new TLatex(binCenterX, binCenterY, Form("%.2e", static_cast<double>(binContent)));
        ////TLatex *latex = new TLatex(binCenterX, binCenterY, Form("%ld", binContent));
        latex->SetTextAlign(22);
        latex->SetTextSize(0.03);
        latex->SetTextAngle(90);
        latex->Draw();
    }
    
   

    TH1D *hPrimVtZ = (TH1D*) file1->Get("hPrimVtZ");
    TCanvas *canvas2 = new TCanvas("canvas2", "hPrimVtZ", 800, 600);
    canvas2->cd();
    hPrimVtZ->SetStats(0);
    hPrimVtZ->SetLineWidth(4);
    hPrimVtZ->Draw("hist");
    hPrimVtZ->GetXaxis()->SetTitle("V_{z}TPC [cm]");
    hPrimVtZ->GetYaxis()->SetTitle("Entries");

    // Draw red dashed lines at x = -60 and x = 60
    TLine *line1 = new TLine(-60, hPrimVtZ->GetMinimum(), -60, hPrimVtZ->GetMaximum());
    line1->SetLineColor(kRed);
    line1->SetLineStyle(10);
    line1->SetLineWidth(2);
    line1->Draw("same");

    TLine *line2 = new TLine(60, hPrimVtZ->GetMinimum(), 60, hPrimVtZ->GetMaximum());
    line2->SetLineColor(kRed);
    line2->SetLineStyle(10);
    line2->SetLineWidth(2);
    line2->Draw("same");

    TH1D *hVzVPDvsVzTPCDiff = (TH1D*) file1->Get("hVzVPDvsVzTPCDiff");
    TCanvas *canvas3 = new TCanvas("canvas3", "hVzVPDvsVzTPCDiff", 800, 600);
    canvas3->cd();
    canvas3->SetLogy();
    hVzVPDvsVzTPCDiff->SetStats(0);
    hVzVPDvsVzTPCDiff->SetLineWidth(4);
    hVzVPDvsVzTPCDiff->GetXaxis()->SetTitle("V_{z}TPC - V_{z}VPD [cm]");
    hVzVPDvsVzTPCDiff->GetYaxis()->SetTitle("Entries");
    hVzVPDvsVzTPCDiff->GetXaxis()->SetRangeUser(-40, 40);
    hVzVPDvsVzTPCDiff->Draw("hist");

    /////Draw red dashed lines at x = -4 and x = 4
    TLine *line3 = new TLine(-4, hVzVPDvsVzTPCDiff->GetMinimum(), -4, hVzVPDvsVzTPCDiff->GetMaximum());
    line3->SetLineColor(kRed);
    line3->SetLineStyle(10);
    line3->SetLineWidth(2);
    line3->Draw("same");

    TLine *line4 = new TLine(4, hVzVPDvsVzTPCDiff->GetMinimum(), 4, hVzVPDvsVzTPCDiff->GetMaximum());
    line4->SetLineColor(kRed);
    line4->SetLineStyle(10);
    line4->SetLineWidth(2);
    line4->Draw("same");
    

    TH2D *hPrimVtXY = (TH2D*) file1->Get("hPrimVtXY");
    TCanvas *canvas4 = new TCanvas("canvas4", "hPrimVtXY", 800, 600);
    canvas4->cd();
    canvas4->SetLogz();
    hPrimVtXY->SetStats(0);
    hPrimVtXY->GetXaxis()->SetTitle("V_{x}TPC [cm]");
    hPrimVtXY->GetXaxis()->SetRangeUser(-3, 3);
    hPrimVtXY->GetYaxis()->SetTitle("V_{y}TPC [cm]");
    hPrimVtXY->GetYaxis()->SetRangeUser(-3, 3);
    hPrimVtXY->GetZaxis()->SetTitle("Entries");
    hPrimVtXY->Draw("colz");

    /*// Draw dashed circle centered around (0,0) with radius 0.3
    TEllipse *circle = new TEllipse(0, 0, 0.3);
    circle->SetFillColorAlpha(0, 0); // Set fill color to transparent
    circle->SetFillStyle(0); // Set fill style to transparent
    circle->SetLineColor(kRed);
    circle->SetLineStyle(10);
    circle->SetLineWidth(6);
    circle->Draw("same");*/
    


// Draw a rectangle with the given edge points
TBox *rectangle = new TBox(-0.3, -0.26, 0.14, 0.02);
    rectangle->SetFillColorAlpha(0, 0); // Set fill color to transparent
    rectangle->SetFillStyle(0);
    rectangle->SetLineColor(kRed);
    rectangle->SetLineStyle(2);
    rectangle->SetLineWidth(2);
    rectangle->Draw("same");
    

    TH2D *hPrimVtXY_evcut = (TH2D*) file1->Get("hPrimVtXY_evcut");
    TCanvas *canvas5 = new TCanvas("canvas5", "hPrimVtXY_evcut", 800, 600);
    canvas5->cd();
    canvas5->SetLogz();
    hPrimVtXY_evcut->SetStats(0);
    hPrimVtXY_evcut->GetXaxis()->SetTitle("V_{x}TPC [cm]");
    hPrimVtXY_evcut->GetXaxis()->SetRangeUser(-3, 3);
    hPrimVtXY_evcut->GetYaxis()->SetTitle("V_{y}TPC [cm]");
    hPrimVtXY_evcut->GetYaxis()->SetRangeUser(-3, 3);
    hPrimVtXY_evcut->GetZaxis()->SetTitle("Entries");
    hPrimVtXY_evcut->Draw("colz");

    TH1D *hPrimVr = new TH1D("Vr", "Vr", hPrimVtXY->GetNbinsX(), -0.5, 5.0);
    TCanvas *canvas6 = new TCanvas("canvas6", "hPrimVr", 800, 600);
    canvas6->SetLogy();
    canvas6->cd();
    hPrimVr->SetStats(0);
    hPrimVr->GetXaxis()->SetTitle("V_{r}TPC [cm]");
    hPrimVr->GetXaxis()->SetRangeUser(-0.25, 3);
    hPrimVr->GetYaxis()->SetTitle("Entries");


  // Fill the TH1D with the sqrt(x^2 + y^2) of the TH2D
  for (int i = 1; i <=hPrimVtXY->GetNbinsX(); i++) {
    for (int j = 1; j <=hPrimVtXY->GetNbinsY(); j++) {
      ///if (hist2D->GetBinContent(i, j) == 0) continue;
      double x =hPrimVtXY->GetXaxis()->GetBinCenter(i);
      double y =hPrimVtXY->GetYaxis()->GetBinCenter(j);
      double value = 0;
      if (x != 0 || y != 0) value = sqrt(x*x + y*y);
      hPrimVr->Fill(value,hPrimVtXY->GetBinContent(i, j));
    }
  }

  hPrimVr->Draw();

    // Draw dashed line at x = 0.3
    TLine *line5 = new TLine(0.3, hPrimVr->GetMinimum(), 0.3, hPrimVr->GetMaximum());
    line5->SetLineColor(kRed);
    line5->SetLineStyle(10);
    line5->SetLineWidth(2);
    ///line5->Draw("same");

    TH2D *hVzVPDvsVzTPC = (TH2D*) file1->Get("hVzVPDvsVzTPC");
    TCanvas *canvas7 = new TCanvas("canvas7", "hVzVPDvsVzTPC", 800, 600);
    canvas7->cd();
    canvas7->SetLogz();
    hVzVPDvsVzTPC->SetStats(0);
    hVzVPDvsVzTPC->GetXaxis()->SetTitle("V_{z}TPC [cm]");
    hVzVPDvsVzTPC->GetXaxis()->SetRangeUser(-200, 200);
    hVzVPDvsVzTPC->GetYaxis()->SetTitle("V_{z}VPD [cm]");
    hVzVPDvsVzTPC->GetYaxis()->SetRangeUser(-500, 500);
    hVzVPDvsVzTPC->GetZaxis()->SetTitle("Entries");
    hVzVPDvsVzTPC->Draw("colz");

    TH2D *hVzVPDvsVzTPC_evcut = (TH2D*) file1->Get("hVzVPDvsVzTPC_evcut");
    TCanvas *canvas8 = new TCanvas("canvas8", "hVzVPDvsVzTPC_evcut", 800, 600);
    canvas8->cd();
    canvas8->SetLogz();
    hVzVPDvsVzTPC_evcut->SetStats(0);
    hVzVPDvsVzTPC_evcut->GetXaxis()->SetTitle("V_{z}TPC [cm]");
    hVzVPDvsVzTPC_evcut->GetXaxis()->SetRangeUser(-200, 200);
    hVzVPDvsVzTPC_evcut->GetYaxis()->SetTitle("V_{z}VPD [cm]");
    hVzVPDvsVzTPC_evcut->GetYaxis()->SetRangeUser(-500, 500);
    hVzVPDvsVzTPC_evcut->GetZaxis()->SetTitle("Entries");
    hVzVPDvsVzTPC_evcut->Draw("colz");

    TH1D *hNHitsFit = (TH1D*) file1->Get("hNHitsFit");
    TCanvas *canvas9 = new TCanvas("canvas9", "hNHitsFit", 800, 600);
    canvas9->cd();
    hNHitsFit->SetStats(0);
    hNHitsFit->SetLineWidth(4);
    hNHitsFit->GetXaxis()->SetTitle("nHitsFit");
    hNHitsFit->GetXaxis()->SetRangeUser(0, 50);
    hNHitsFit->GetYaxis()->SetTitle("Entries");
    hNHitsFit->Draw("hist");

    // Draw red dashed lines at x = 18
    TLine *line6 = new TLine(18, hNHitsFit->GetMinimum(), 15, hNHitsFit->GetMaximum());
    line6->SetLineColor(kRed);
    line6->SetLineStyle(10);
    line6->SetLineWidth(2);
    line6->Draw("same");

    TH1D *hNHitsFitnHitsMax = (TH1D*) file1->Get("hNHitsFitnHitsMax");
    TCanvas *canvas10 = new TCanvas("canvas10", "hNHitsFitnHitsMax", 800, 600);
    canvas10->cd();
    hNHitsFitnHitsMax->SetStats(0);
    hNHitsFitnHitsMax->SetLineWidth(4);
    hNHitsFitnHitsMax->GetXaxis()->SetTitle("nHitsFit/nHitsMax");
    hNHitsFitnHitsMax->GetXaxis()->SetRangeUser(0, 1.2);
    hNHitsFitnHitsMax->GetYaxis()->SetTitle("Entries");
    hNHitsFitnHitsMax->Draw("hist");

    // Draw red dashed lines at x = 0.52
    TLine *line7 = new TLine(0.52, hNHitsFitnHitsMax->GetMinimum(), 0.52, hNHitsFitnHitsMax->GetMaximum());
    line7->SetLineColor(kRed);
    line7->SetLineStyle(10);
    line7->SetLineWidth(2);
    line7->Draw("same");

    TH1D *hDca = (TH1D*) file1->Get("hDca");
    TCanvas *canvas11 = new TCanvas("canvas11", "hDca", 800, 600);
    canvas11->cd();
    hDca->SetStats(0);
    hDca->SetLineWidth(4);
    hDca->GetXaxis()->SetTitle("gDCA [cm]");
    hDca->GetXaxis()->SetRangeUser(0, 4);
    hDca->GetYaxis()->SetTitle("Entries");
    hDca->Draw("hist");

    // Draw red dashed lines at x = 1.5
    TLine *line8 = new TLine(1.5, hDca->GetMinimum(), 1.5, hDca->GetMaximum());
    line8->SetLineColor(kRed);
    line8->SetLineStyle(10);
    line8->SetLineWidth(2);
    line8->Draw("same");

    //Draw blue dashed lines at x = 3.0
    TLine *line12 = new TLine(3.0, hDca->GetMinimum(), 3.0, hDca->GetMaximum());
    line12->SetLineColor(kGreen);
    line12->SetLineStyle(10);
    line12->SetLineWidth(2);
    line12->Draw("same");

    // Create a legend
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0, 0); // Set fill color to transparent

    // Add entries to the legend
    legend->AddEntry(line8, "K and Soft #pi", "l");
    legend->AddEntry(line12, "D^{0} #pi matched to BHT1 Trigger", "l");

    // Draw the legend
    legend->Draw();


    TH1D *hpT_tr = (TH1D*) file1->Get("hpT_tr");
    TCanvas *canvas12 = new TCanvas("canvas12", "hpT_tr", 800, 600);
    canvas12->cd();
    hpT_tr->SetStats(0);
    hpT_tr->SetLineWidth(4);
    hpT_tr->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    hpT_tr->GetXaxis()->SetRangeUser(0, 3);
    hpT_tr->GetYaxis()->SetTitle("Entries");
    hpT_tr->Draw("hist");

    // Draw red dashed lines at x = 0.2
    TLine *line9 = new TLine(0.2, hpT_tr->GetMinimum(), 0.2, hpT_tr->GetMaximum());
    line9->SetLineColor(kGreen);
    line9->SetLineStyle(10);
    line9->SetLineWidth(2);
    line9->Draw("same");

    // Draw red dashed lines at x = 0.1
    TLine *line13 = new TLine(0.1, hpT_tr->GetMinimum(), 0.1, hpT_tr->GetMaximum());
    line13->SetLineColor(kRed);
    line13->SetLineStyle(10);
    line13->SetLineWidth(2);
    line13->Draw("same");

    //create a legend
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend1->SetBorderSize(0);
    legend1->SetFillColorAlpha(0, 0); // Set fill color to transparent

    // Add entries to the legend
    legend1->AddEntry(line9, "D^{0} daughters", "l");
    legend1->AddEntry(line13, "Soft #pi", "l");

    // Draw the legend
    legend1->Draw();

    

    TH1D *hEta_tr = (TH1D*) file1->Get("hEta_tr");
    TCanvas *canvas13 = new TCanvas("canvas13", "hEta_tr", 800, 600);
    canvas13->cd();
    hEta_tr->SetStats(0);
    hEta_tr->SetLineWidth(4);
    hEta_tr->GetXaxis()->SetTitle("#eta");
    hEta_tr->GetXaxis()->SetRangeUser(-2, 2);
    hEta_tr->GetYaxis()->SetTitle("Entries");
    hEta_tr->Draw("hist");

    // Draw red dashed lines at x = -1 and x = 1
    TLine *line10 = new TLine(-1, hEta_tr->GetMinimum(), -1, hEta_tr->GetMaximum());
    line10->SetLineColor(kRed);
    line10->SetLineStyle(10);
    line10->SetLineWidth(2);
    line10->Draw("same");

    TLine *line11 = new TLine(1, hEta_tr->GetMinimum(), 1, hEta_tr->GetMaximum());
    line11->SetLineColor(kRed);
    line11->SetLineStyle(10);
    line11->SetLineWidth(2);
    line11->Draw("same");
    





























        
        
    }

    



