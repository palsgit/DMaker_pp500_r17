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

void runQAAnaMaker(
    const char*  inputFile,
    const Char_t *outputFile,  
    const Char_t *badRunListFileName) {
    string SL_version = "SL22b";
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
    hfCuts->addTriggerId(570001); //VPDMB-30


    hfCuts->setCutPrimaryDCAtoVtxMax(1.5);
    hfCuts->setCutVzMax(30.);
    hfCuts->setCutVzVpdVzMax(100.);
    hfCuts->setCutNHitsFitMin(20);
    hfCuts->setCutNHitsFitnHitsMax(0.52);
    hfCuts->setCutRequireHFT(false);
    hfCuts->setHybridTof(false); // Does nothing (Proton PID)
    hfCuts->setHybridTofKaon(true);
    hfCuts->setHybridTofPion(true);
    hfCuts->setCheckHotSpot(false);

    hfCuts->setCutTPCNSigmaPion(3.0);
    hfCuts->setCutTPCNSigmaKaon(2.0);
//    hfCuts->setCutTOFDeltaOneOverBetaKaon(0.03);
//    hfCuts->setCutTOFDeltaOneOverBetaPion(0.03);
//    hfCuts->setCutPtMin(0.15);

    hfCuts->setCutDcaMin(0.002,StHFCuts::kPion);
    hfCuts->setCutDcaMin(0.002,StHFCuts::kKaon);

    hfCuts->setHybridTofBetterBetaCuts(false); // Does nothing
    hfCuts->setHybridTofBetterBetaCutsKaon(true);
    hfCuts->setHybridTofBetterBetaCutsPion(false);

    //Single track pt
    hfCuts->setCutPtRange(0.15,50.0,StHFCuts::kPion); //0.2 , 50.0
    hfCuts->setCutPtRange(0.15,50.0,StHFCuts::kKaon); //0.2, 50.0
    //TPC setters
    hfCuts->setCutTPCNSigmaPion(10.); //3
    hfCuts->setCutTPCNSigmaKaon(10.); //3
    //TOF setters, need to set pt range as well
    hfCuts->setCutTOFDeltaOneOverBeta(0.1, StHFCuts::kKaon); // v podstate 5 sigma; nastavene = f * (sigmaTOF), sigma TOF je 0.013
    hfCuts->setCutPtotRangeHybridTOF(0.2,50.0,StHFCuts::kKaon);
    hfCuts->setCutTOFDeltaOneOverBeta(0.1, StHFCuts::kPion); // v podstate 6 sigma
    hfCuts->setCutPtotRangeHybridTOF(0.2,50.0,StHFCuts::kPion);


    float dcaDaughtersMax = 0.2;  // maximum
    float decayLengthMin  = 0.000; // minimum
    float decayLengthMax  = 999999; //std::numeric_limits<float>::max();
    float cosThetaMin     = -20.;   // minimum
    float minMass         = 0.6;
    float maxMass         = 2.6;
    float pairDcaMax      = 99.9;


    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, sInputFile, "picoDstMaker"); //for local testing only (akorát že vůbec)
//    StPicoDstMaker* picoDstMaker = new StPicoDstMaker(static_cast<StPicoDstMaker::PicoIoMode>(StPicoDstMaker::IoRead), inputFile, "picoDstMaker");
    StPicoQAMaker* PicoQAMaker = new StPicoQAMaker("picoQAMaker", picoDstMaker, outputFile);
    PicoQAMaker->setHFBaseCuts(hfCuts);

//    clock_t start = clock(); // getting starting time
    chain->Init();
    Int_t nEvents = picoDstMaker->chain()->GetEntries();
    cout << " Total entries = " << nEvents << endl;

    for (Int_t i=0; i<nEvents; ++i) {
//        if(i%10==0)       cout << "Working on eventNumber " << i << endl;
        chain->Clear();
        int iret = chain->Make(i);
        if (iret) { cout << "Bad return code!" << iret << endl; break;}
    }

    chain->Finish();
  //  double duration = (double) (clock() - start) / (double) CLOCKS_PER_SEC;
    cout << "****************************************** " << endl;
    cout << "Work done, total number of events  " << nEvents << endl;
//    cout << "Time needed " << duration << " s" << endl;
    delete chain;
}
