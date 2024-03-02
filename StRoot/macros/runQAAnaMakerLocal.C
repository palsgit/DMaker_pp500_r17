#ifndef __CINT__

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "StChain/StMaker.h"
#include "StChain/StChain.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoHFMaker/StHFMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "macros/loadSharedHFLibraries.C"
#include <iostream>
#include <ctime>
#include <cstdio>
#include "StPicoQAMaker/StPicoQAMaker.h"

using namespace std;

#else
class StChain;
#endif

class StChain;
class StPicoDstMaker;
class StPicoQAMaker;
class StMaker;

void runQAAnaMakerLocal(
			const Char_t *inputFile="./picoLists/runs_local_test.list",
			const Char_t *outputFile="outputLocal",
			const Char_t *badRunListFileName = "./picoLists/picoList_bad.list") {
  string SL_version = "SL22b";
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) {
      cout<<"Environment Star Library does not match the requested library. Exiting..."<<endl;
      exit(1);
  }

#ifdef __CINT__
    gROOT->LoadMacro("loadSharedHFLibraries.C");
  loadSharedHFLibraries();
#endif

  Int_t nEvents = 2000000;
  StChain *chain = new StChain();
  TString sInputFile(inputFile);

  if (!sInputFile.Contains(".list") && !sInputFile.Contains("picoDst.root")) {
    cout << "No input list or picoDst root file provided! Exiting..." << endl;
    exit(1);
  }
  cout<<"event stuff set"<<endl;


      StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");

    hfCuts->setBadRunListFileName(badRunListFileName);
    
    //event cuts
    hfCuts->addTriggerId(570001); //VPDMB-30
    hfCuts->showTriggers();
    hfCuts->setnMatchedFast(0);
    hfCuts->setCutVzVpdVzMax(10.);
    hfCuts->setCutVzMax(50.);
    ////hfCuts->setCutVrMax(0.25);
    hfCuts->setCutVxMax(1.);
    hfCuts->setCutVyMax(1.);
    hfCuts->setCutVxMin(-1.);
    hfCuts->setCutVyMin(-1.);

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


//
    /* float dcaDaughtersMax = 10.;  // maximum toto ide
    float decayLengthMin  = 0.00000000; // minimum
    float decayLengthMax  = 9999999.;  //std::numeric_limits<float>::max(); toto ide (cutuje)
    float cosThetaMin     = -20.;   // minimum
    float minMass         = 0.1;
    float maxMass         = 3.5;
    float pairDcaMax      = 99.9;

    hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass, pairDcaMax);
    */

/*    hfCuts->setCutSecondaryPairPtBin(1,      2,              0.007,          0.012,         0.5,      0.005,    0.009, 0.007);
    hfCuts->setCutSecondaryPairPtBin(2,      3,              0.016,          0.003,         0.5,      0.0065,   0.009, 0.01);
    hfCuts->setCutSecondaryPairPtBin(3,      5,              0.015,          0.009,         0.6,      0.0064,   0.0064, 0.0076);*/

    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, sInputFile, "picoDstMaker"); //for local testing only (akorát že vůbec)
  //  StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
  StPicoQAMaker* PicoQAAnaMaker = new StPicoQAMaker("picoQAAnaMaker", picoDstMaker, outputFile);
  PicoQAAnaMaker->setHFBaseCuts(hfCuts);

//  clock_t start = clock(); // getting starting time
  chain->Init();
  
  int total = picoDstMaker->chain()->GetEntries();
  cout << " Total entries = " << total << endl;
  if(nEvents>total) nEvents = total;

  for (Int_t i=0; i<nEvents; i++) {
    if(i%10==0)       cout << "Working on eventNumber " << i << endl;

    chain->Clear();
    int iret = chain->Make(i);

    if (iret) { cout << "Bad return code!" << iret << endl; break;}
    }
  
  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  chain->Finish();
//  double duration = (double) (clock() - start) / (double) CLOCKS_PER_SEC;
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
//  cout << "Time needed " << duration << " s" << endl;

  delete chain;

}

