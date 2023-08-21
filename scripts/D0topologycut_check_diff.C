#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TROOT.h>

void drawHistogramsFromRootFile(const char* inputfileName1, const char* inputfileName2, const char* outputfileName) {
    /////gROOT->SetBatch(); // Run in batch mode (no GUI)
    
    const int numCanvases = 13;
    ////const int numHistograms = 4;
    const int numHistograms = 2;
    
    TFile* inputFile1 = new TFile(inputfileName1);
    TFile* inputFile2 = new TFile(inputfileName2); // Open the input ROOT file
    
    if (inputFile1->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputfileName1 << std::endl;
        return;
    }

    if (inputFile2->IsZombie()) {
        std::cerr << "Error: Cannot open input file " << inputfileName2 << std::endl;
        return;
    }
   
    TFile* outputFile = new TFile(outputfileName, "RECREATE");

    TString histname; 
    /////TString spectraname = "unlike"; 
    /////Double_t scalefactor;
    /////Int_t nRebin;

    ////if (inputfileName1 == "pp500_D0_prereqgDCA1.5_all_y0.4eta1.0.root") scalefactor = 0.942514;
    /////if (inputfileName2 == "pp500_D0_prereqgDCA1.5_all_y0.4to1.0_eta1.0.root") scalefactor = 1.00253;
    
    // Create and configure canvases
    //////TCanvas* canvases[numCanvases];
    for (int i = 0; i < numCanvases; ++i) {
        /////canvases[i] = new TCanvas(Form("canvas%d", i), Form("Canvas %d", i), 800, 600);
        ////canvases[i]->Divide(1, 1); // Split canvas into 1x1 grid

        if (i == 0) histname = "dcaDaughters_";
        if (i == 1) histname = "D_theta_";
        if (i == 2) histname = "cosTheta_";
        if (i == 3) histname = "D_decayL_";
        if (i == 4) histname = "dcaD0ToPv_";
        if (i == 5) histname = "D_cosThetaStar_";
        if (i == 6) histname = "D_pt_";
        if (i == 7) histname = "D_mass_";
        if (i == 8) histname = "D_rapidity_";
        if (i == 9) histname = "D_phi_";
        if (i == 10) histname = "D_eta_";
        if (i == 11) histname = "pi_DCA_";
        if (i == 12) histname = "k_DCA_";

        /////canvases[i]->cd();

        

        TH1D* unlike_histogram = dynamic_cast<TH1D*>(inputFile1->Get((histname + "rot")));
        /*float nBinsInvMass = unlike_histogram->GetNbinsX();
        
        ///spectraname = "rot";
        TH1D* rot_histogram = dynamic_cast<TH1D*>(inputFile1->Get((histname + "rot")));
        rot_histogram->Scale(0.942514);
        
        TH1D *hist_topo = new TH1D(histname + "0_4",histname + "0_4",nBinsInvMass,(unlike_histogram->GetBinLowEdge(1)),(unlike_histogram->GetBinLowEdge(nBinsInvMass+1)));
	    hist_topo->Sumw2();
	    for(int i=0; i<hist_topo->GetNbinsX(); i++){
		hist_topo->SetBinContent(i+1,unlike_histogram->GetBinContent(1+i) - rot_histogram->GetBinContent(1+i));
		hist_topo->SetBinError(i+1,TMath::Sqrt(unlike_histogram->GetBinError(1+i)*unlike_histogram->GetBinError(1+i)+
											  rot_histogram->GetBinError(1+i)*rot_histogram->GetBinError(1+i)));
	   }
       */

       /*if (i == 0) hist_topo->Rebin(10);
       if (i == 1) hist_topo->Rebin(15);
       if (i == 2) hist_topo->Rebin(12);
       if (i == 3) hist_topo->Rebin(20);
       if (i == 4) hist_topo->Rebin(1);
       if (i == 5) hist_topo->Rebin(1);
       if (i == 8) hist_topo->Rebin(24);
       if (i == 9) hist_topo->Rebin(4);

       hist_topo->SetLineColor(1); // Set histogram line color
       hist_topo->SetMarkerColor(1); 
       hist_topo->SetMarkerStyle(8);
       hist_topo->Draw("ap");

       hist_topo->Write();*/
       unlike_histogram->SetName(histname + "rot" + "_0.4");
       unlike_histogram->Scale(1.0/unlike_histogram->GetEntries());
       unlike_histogram->SetLineColor(1); // Set histogram line color
       unlike_histogram->SetMarkerColor(1); 
       unlike_histogram->SetMarkerStyle(8);
       unlike_histogram->Draw("ap");

       unlike_histogram->Write();


        unlike_histogram = dynamic_cast<TH1D*>(inputFile2->Get((histname + "rot")));
        /////float nBinsInvMass = unlike_histogram->GetNbinsX();
        
        ///spectraname = "rot";
        /*rot_histogram = dynamic_cast<TH1D*>(inputFile2->Get((histname + "rot")));
        rot_histogram->Scale(0.963398);
        
        hist_topo = new TH1D(histname + "1_0",histname + "1_0",nBinsInvMass,(unlike_histogram->GetBinLowEdge(1)),(unlike_histogram->GetBinLowEdge(nBinsInvMass+1)));
	    hist_topo->Sumw2();
	    for(int i=0; i<hist_topo->GetNbinsX(); i++){
		hist_topo->SetBinContent(i+1,unlike_histogram->GetBinContent(1+i) - rot_histogram->GetBinContent(1+i));
		hist_topo->SetBinError(i+1,TMath::Sqrt(unlike_histogram->GetBinError(1+i)*unlike_histogram->GetBinError(1+i)+
											  rot_histogram->GetBinError(1+i)*rot_histogram->GetBinError(1+i)));
	   }
       */
       
       	
       /*if (i == 0) hist_topo->Rebin(10);
       if (i == 1) hist_topo->Rebin(15);
       if (i == 2) hist_topo->Rebin(12);
       if (i == 3) hist_topo->Rebin(20);
       if (i == 4) hist_topo->Rebin(1);
       if (i == 5) hist_topo->Rebin(1);
       if (i == 8) hist_topo->Rebin(24);
       if (i == 9) hist_topo->Rebin(4);*/

       /*hist_topo->SetLineColor(2); // Set histogram line color
       hist_topo->SetMarkerColor(2); 
       hist_topo->SetMarkerStyle(8);
       hist_topo->Draw("apsame");*/
       
       unlike_histogram->SetName(histname + "rot" + "_1.0");
       unlike_histogram->Scale(1.0/unlike_histogram->GetEntries());
       unlike_histogram->SetLineColor(2); // Set histogram line color
       unlike_histogram->SetMarkerColor(2); 
       unlike_histogram->SetMarkerStyle(8);
       unlike_histogram->Draw("apsame");

       TLegend *leg = new TLegend(0.6,0.7,0.88,0.88,"");
	   leg->SetFillColor(10);
	   leg->SetTextFont(42);


       unlike_histogram->Write();

        
       /* for (int j = 0; j < numHistograms; ++j) {
            canvases[i]->cd(); // Switch to the current canvas
            if (j == 0) spectraname = "unlike";
            if (j == 1) spectraname = "rot";
            /*if (j == 1) spectraname = "mixed";
            if (j == 2) spectraname = "rot";
            if (j == 3) spectraname = "like";
            
            TH1D* histogram = dynamic_cast<TH1D*>(inputFile->Get((histname + spectraname))); // Get histogram from file
            if (i == 1) histogram->Rebin(30);
            if (i == 2) histogram->Rebin(6);
            if (i == 3) histogram->Rebin(80);
            if (i == 4) histogram->Rebin(24);
            if (i == 5) histogram->Rebin(6);
            if (i == 8) histogram->Rebin(24);
            if (i == 9) histogram->Rebin(4);
            if (!histogram) {
                cout << histname + spectraname << endl;
               std::cerr << "Error: Invalid histogram provided." << std::endl;
                return;
              }

            if (histogram) {
                histogram->SetLineColor(j + 1); // Set histogram line color
                histogram->SetMarkerColor(j + 1); 
                histogram->SetMarkerStyle(8);
                histogram->Draw(j == 0 ? "p" : "psame"); // Draw histograms on top of each other
            }
            leg->AddEntry(histogram, histname + spectraname, "lp");
            leg->Draw();
        }
        */
        ////canvases[i]->Update();
        ////canvases[i]->SaveAs(histname+".png");
        ////canvases[i]->Write();
    }
    

    /////gApplication->Run();


    ////inputFile->Close(); // Close the input ROOT file
    outputFile->Close();
    outputFile->Save();
}

void D0topologycut_check_diff() {
    const char* inputfileName1 = "pp500_D0_prereqgDCA1.5_all_y0.4_rotonly.root";
    const char* inputfileName2 = "pp500_D0_prereqgDCA1.5_all_y1.0_rotonly.root"; // Replace with your input ROOT file name
    const char* outputFileName = "output_rot_only_difference_background_scaled.root";
    drawHistogramsFromRootFile(inputfileName1, inputfileName2, outputFileName);
}
