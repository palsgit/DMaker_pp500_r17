#ifndef __CINT__

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "StChain/StMaker.h"
#include "StChain/StChain.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoMixedEventMaker/StPicoMixedEventMaker.h"
#include "macros/loadSharedHFLibraries.C"
#include <iostream>
#include <ctime>
#include <cstdio>
#include "StPicoD0AnaMaker/StPicoD0AnaMaker.h"
#include "StPicoCutsBase/StPicoCutsBase.h"

using namespace std;

#else
class StChain;
#endif

void runPicoMixedEventLocal(
			const Char_t *inputFile="./picoLists/runs_local_test.list",
			const Char_t *outputFile="outputLocal",
			const Char_t *badRunListFileName = "./picoLists/picoList_bad.list") {

  string SL_version = "SL22b";
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) {
      cout<<"Environment Star Library does not match the requested library in runPicoHFMyAnaMaker.C. Exiting..."<<endl;
      exit(1);
  }

#ifdef __CINT__
    gROOT->LoadMacro("loadSharedHFLibraries.C");
  loadSharedHFLibraries();
#endif

  gROOT->LoadMacro("loadSharedHFLibraries.C");
  loadSharedHFLibraries();
  StChain *chain = new StChain();

  TString sInputFile(inputFile);

  if (!sInputFile.Contains(".list") && !sInputFile.Contains("picoDst.root")) {
    cout << "No input list or picoDst root file provided! Exiting..." << endl;
    exit(1);
  }

  StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");

    hfCuts->setBadRunListFileName(badRunListFileName);
    
    //event cuts
    hfCuts->addTriggerId(570001); //VPDMB-30
    hfCuts->showTriggers();
    hfCuts->setnMatchedFast(0);
    hfCuts->setCutVzVpdVzMax(10.);
    hfCuts->setCutVzMax(50.);
    hfCuts->setCutVrMax(0.25);

    //track cuts
    hfCuts->setCutNHitsFitMin(17);
    hfCuts->setCutNHitsFitnHitsMax(0.52);
    hfCuts->setCutPrimaryDCAtoVtxMax(1.5); //was 2.0 in DTlusty thesis
    hfCuts->setCutPtMin(0.20);
    hfCuts->setCutEtaMax(1.0);
    hfCuts->setCutEtaMin(-1.0);
    
    /////hfCuts->setCutVzVpdVzMax(100.);
    
    
    
    /*
    hfCuts->setCutRequireHFT(false); //// Usable only for d+Au analysis
    hfCuts->setHybridTof(false); //// Does nothing (Proton PID)
    */
    hfCuts->setHybridTofKaon(false); //// This cut and the one below work for the analysis without BEMC
    hfCuts->setHybridTofPion(false); ////
    ///hfCuts->setCheckHotSpot(false);
    

    hfCuts->setCutTPCNSigmaPionMax(2.4);
    hfCuts->setCutTPCNSigmaPionMin(-2.4);
    hfCuts->setCutTPCNSigmaKaonMax(2.4);
    hfCuts->setCutTPCNSigmaKaonMin(-2.4);

    //hfCuts->setCutDcaMin(0.002,StHFCuts::kPion);//was not mentioned in DTlusty thesis
    //hfCuts->setCutDcaMin(0.002,StHFCuts::kKaon);//was not mentioned in DTlusty thesis
    
    
    hfCuts->setCutTOFNSigmaPionMax(2.4);
    hfCuts->setCutTOFNSigmaPionMin(-1.6);
    /*hfCuts->setCutTOFNSigmaKaon(3.0);
    hfCuts->setCutTOFDeltaOneOverBetaKaon(0.03);
    hfCuts->setCutTOFDeltaOneOverBetaPion(0.03);
    

    

    hfCuts->setHybridTofBetterBetaCuts(false); // Does nothing
    */
    hfCuts->setHybridTofBetterBetaCutsKaon(true); //// This cut and the one below work for the analysis without BEMC, it turns on cuts of TOF 1/beta in a shape of a function
    hfCuts->setHybridTofBetterBetaCutsPion(false); ////
    
    hfCuts->setHybridTofWithBEMC(false);

   /* float dcaDaughtersMax = 10.;  // maximum toto ide
    float decayLengthMin  = 0.000000000; // minimum
    float decayLengthMax  = 99999999.;  //std::numeric_limits<float>::max(); toto ide (cutuje)
    float cosThetaMin     = -20.;   // minimum
    float minMass         = 0.1;
    float maxMass         = 3.5;
    float pairDcaMax      = 99.9;


  hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass, pairDcaMax);
  */

  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, sInputFile, "picoDstMaker"); //for local testing only (akorát že vůbec)
//  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
  StPicoMixedEventMaker* picoMixedEventMaker = new StPicoMixedEventMaker("picoMixedEventMaker", picoDstMaker, hfCuts, outputFile);
  picoMixedEventMaker->setBufferSize(10);

//  clock_t start = clock(); // getting starting time
  chain->Init();
  Int_t nEvents = picoDstMaker->chain()->GetEntries();
  cout << "Total entries = " << nEvents << endl;

  for (Int_t i=0; i<nEvents; i++) {
    if(i%10==0)  cout << "Working on eventNumber " << i << endl;

    chain->Clear();
    int iret = chain->Make(i);
    if (iret) { cout << "Bad return code!" << iret << endl; break;}
  }

  chain->Finish();
//  double duration = (double) (clock() - start) / (double) CLOCKS_PER_SEC;
  cout << "****************************************** " << endl;
  cout << "Work done, total number of events  " << nEvents << endl;
//  cout << "Time needed " << duration << " s" << endl;
  delete chain;
}

