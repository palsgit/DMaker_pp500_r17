void DrawTProfiles() {
  // Open the ROOT files containing the TLists with the TProfiles
   TFile* file1 = TFile::Open("output_all.root");
   /////TFile* file2 = TFile::Open("output_2xfastmatched_events_gDCA_1.5cm_no_minDCA_cut.root");
   ///TFile* file3 = TFile::Open("output_DTlusty_thesis_2xfastmatched_event_matchedtracks.root");
  /////TFile* file4 = TFile::Open("output_no_matching_vzvpd-vztpc100.root");

  // Retrieve the TLists with the TProfiles
  TList* list1 = (TList*) file1->Get("picoD0AnaMaker");
  /////TList* list2 = (TList*) file2->Get("picoD0AnaMaker");
  ///TList* list3 = (TList*) file3->Get("picoD0AnaMaker");
  /////TList* list4 = (TList*) file4->Get("picoD0AnaMaker");

  // Create the canvas
  TCanvas* canvas = new TCanvas("canvas", "TProfiles on canvas", 800, 600);

  // Draw the first TProfile from the first list
  TProfile* profile1 = (TProfile*) list1->FindObject("prof_ntracksPionTPC_vs_BBCx");
  profile1->SetMarkerStyle(20);
  profile1->GetYaxis()->SetTitle("track multiplicity");
  profile1->GetXaxis()->SetRangeUser(1500,4500);
  profile1->GetYaxis()->SetRangeUser(0,18);
  profile1->SetMarkerColor(kRed);
  profile1->SetMarkerSize(1);
  profile1->SetLineColor(kRed);
  profile1->SetLineWidth(1);
  profile1->Draw();
  
  // Draw the remaining TProfiles from the first list on the same canvas
  
  
  // Draw the TProfiles from the second list on the same canvas
  TProfile* profile2 = (TProfile*) list1->FindObject("prof_ntracksKaonTPC_vs_BBCx");
  profile2->SetMarkerStyle(20);
  profile2->SetMarkerSize(1);
  profile2->SetMarkerColor(kBlue);
  profile2->SetLineColor(kBlue);
  profile2->SetLineWidth(1);
  profile2->Draw("same");
    
  
  
  /*TProfile* profile3 = (TProfile*) list3->FindObject("prof_ntracksKaon_vs_BBCx");
  profile3->SetMarkerStyle(20);
  profile3->SetMarkerSize(1);
  profile3->SetLineColor(kGreen);
  profile3->SetMarkerColor(kGreen);
  profile3->SetLineWidth(1);
  profile3->Draw("same");
  

  TProfile* profile4 = (TProfile*) list4->FindObject("prof_ntracksKaon_vs_BBCx");
  profile4->SetMarkerStyle(20);
  profile4->SetMarkerSize(1);
  profile4->SetLineColor(kBlack);
  profile4->SetMarkerColor(kBlack);
  profile4->SetLineWidth(1);
  profile4->Draw("same");
  */
  

   TLegend* legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend1->SetFillColor(0);
  legend1->SetBorderSize(0);
  //legend1->AddEntry(profile4, "no matching in fast detectors: |Vz_{vpd} -  Vz_{tpc}| < 100 cm.", "l");
  //legend1->AddEntry(profile1, "No matching in fast detectors", "l");
  legend1->AddEntry(profile1, "Pions", "p");
  legend1->AddEntry(profile2, "Kaons", "p");
  ///legend1->AddEntry(profile2, "Events with at least 2 tracks with BEMC or TOF matched hit", "l");
  ///legend1->AddEntry(profile3, "Tracks with TOF/BEMC matched hit in 2xfastmatched Events", "l");
  
  legend1->Draw();
  // Update the canvas
  canvas->Update();
  
  // Close the ROOT files
  file1->Close();
  ////file2->Close();
  //file3->Close();
  ///file4->Close();
}
