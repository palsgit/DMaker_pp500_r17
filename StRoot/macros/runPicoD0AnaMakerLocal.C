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
#include "macros/loadSharedHFLibraries.C"
//#include "StPicoMixedEventMaker/StPicoMixedEventMaker.h"
#include <iostream>
#include <ctime>
#include <cstdio>
#include "StPicoD0AnaMaker/StPicoD0AnaMaker.h"
//#include "StPicoQAMaker/StPicoQAMaker.h"

using namespace std;

#else
class StChain;
#endif

void runPicoD0AnaMakerLocal(
			const Char_t *inputFile="./picoLists/runs_local_test.list",
			const Char_t *outputFile="outputLocal",
			const Char_t *badRunListFileName = "./picoLists/picoList_bad.list") {
    string SL_version = "SL19c";
    string env_SL = getenv ("STAR");
    if (env_SL.find(SL_version)==string::npos) {
        cout<<"Environment Star Library does not match the requested library in run**.C. Exiting..."<<endl;
        exit(1);
    }



#ifdef __CINT__
    gROOT->LoadMacro("loadSharedHFLibraries.C");
  loadSharedHFLibraries();
#endif

    StChain *chain = new StChain();
    TString sInputFile(inputFile);

    if (!sInputFile.Contains(".list") && !sInputFile.Contains("picoDst.root")) {
        cout << "No input list or picoDst root file provided! Exiting..." << endl;
        exit(1);
    }

    StHFCuts* hfCuts = new StHFCuts("hfBaseCuts");

    hfCuts->setBadRunListFileName(badRunListFileName);
    
    //event cuts
    hfCuts->addTriggerId(480001); //VPDMB-5-ssd
    hfCuts->showTriggers();
    hfCuts->setCutVzVpdVzMax(6.);
    hfCuts->setCutVzMax(50.);

    //track cuts
    hfCuts->setCutNHitsFitMin(20);
    hfCuts->setCutNHitsFitnHitsMax(0.52);
    hfCuts->setCutPrimaryDCAtoVtxMax(1.5); //was 2.0 in DTlusty thesis
    hfCuts->setCutPtMin(0.16);
    hfCuts->setCutEtaMax(1.0);
    
    /////hfCuts->setCutVzVpdVzMax(100.);
    
    
    
    /*
    hfCuts->setCutRequireHFT(false); //// Usable only for d+Au analysis
    hfCuts->setHybridTof(false); //// Does nothing (Proton PID)
    */
    hfCuts->setHybridTofKaon(false); //// This cut and the one below work for the analysis without BEMC
    hfCuts->setHybridTofPion(false); ////
    ///hfCuts->setCheckHotSpot(false);
    

    hfCuts->setCutTPCNSigmaPion(2.0);
    hfCuts->setCutTPCNSigmaKaon(2.0);

    //hfCuts->setCutDcaMin(0.002,StHFCuts::kPion);//was not mentioned in DTlusty thesis
    //hfCuts->setCutDcaMin(0.002,StHFCuts::kKaon);//was not mentioned in DTlusty thesis
    
    
    hfCuts->setCutTOFNSigmaPion(2.0);
    hfCuts->setCutTOFNSigmaKaon(2.0);
    /*hfCuts->setCutTOFDeltaOneOverBetaKaon(0.03);
    hfCuts->setCutTOFDeltaOneOverBetaPion(0.03);
    

    

    hfCuts->setHybridTofBetterBetaCuts(false); // Does nothing
    */
    hfCuts->setHybridTofBetterBetaCutsKaon(false); //// This cut and the one below work for the analysis without BEMC, it turns on cuts of TOF 1/beta in a shape of a function
    hfCuts->setHybridTofBetterBetaCutsPion(false); ////
    
    hfCuts->setHybridTofWithBEMC(false);
    

//
    float dcaDaughtersMax = 10.;  // maximum toto ide
    float decayLengthMin  = 0.00000000; // minimum
    float decayLengthMax  = 9999999.;  //std::numeric_limits<float>::max(); toto ide (cutuje)
    float cosThetaMin     = -20.;   // minimum
    float minMass         = 0.1;
    float maxMass         = 3.5;
    float pairDcaMax      = 99.9;

    hfCuts->setCutSecondaryPair(dcaDaughtersMax, decayLengthMin, decayLengthMax, cosThetaMin, minMass, maxMass, pairDcaMax);

/*    hfCuts->setCutSecondaryPairPtBin(1,      2,              0.007,          0.012,         0.5,      0.005,    0.009, 0.007);
    hfCuts->setCutSecondaryPairPtBin(2,      3,              0.016,          0.003,         0.5,      0.0065,   0.009, 0.01);
    hfCuts->setCutSecondaryPairPtBin(3,      5,              0.015,          0.009,         0.6,      0.0064,   0.0064, 0.0076);*/

    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, sInputFile, "picoDstMaker"); //for local testing only (akorát že vůbec)
//    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
    StPicoD0AnaMaker* PicoD0AnaMaker = new StPicoD0AnaMaker("picoD0AnaMaker", picoDstMaker, outputFile);
    PicoD0AnaMaker->workWithRefit(false);
    PicoD0AnaMaker->setHFBaseCuts(hfCuts);



//    StPicoMixedEventMaker* picoMixedEventMaker = new StPicoMixedEventMaker("picoMixedEventMaker", picoDstMaker, hfCuts, outputFile, inputFile);
//    picoMixedEventMaker->setBufferSize(7);
    
//    clock_t start = clock(); // getting starting time
    chain->Init();
    Int_t nEvents = picoDstMaker->chain()->GetEntries();
    cout << " Total entries = " << nEvents << endl;
    for (Int_t i=0; i<nEvents; ++i) {
        if(i%10==0)       cout << "Working on eventNumber " << i << endl;
        chain->Clear();
        int iret = chain->Make(i);
        if (iret) { cout << "Bad return code!" << iret << endl; break;}
    }
    
    chain->Finish();
//    double duration = (double) (clock() - start) / (double) CLOCKS_PER_SEC;
    cout << "****************************************** " << endl;
    cout << "Work done, total number of events  " << nEvents << endl;
//    cout << "Time needed " << duration << " s" << endl;
    delete chain;
}
