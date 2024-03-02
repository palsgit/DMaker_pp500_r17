#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBEmcPidTraits.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "phys_constants.h"
#include "StPicoD0AnaMaker.h"
#include "TComplex.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include <iostream>
#include <fstream>
#include <vector>

#include "Math/Vector4D.h"

#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"

/////#include "StUpsCandidate.h"

#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"

#include "StPicoHFMaker/StHFTriplet.h"
#include "StPicoHFMaker/StHFRotPair.h"

#include "StPicoCutsBase/StPicoCutsBase.h"

//#include "StPicoKFVertexFitter/StPicoKFVertexFitter.h

ClassImp(StPicoD0AnaMaker)

using namespace std;

const int nptBins=3;
//float const bdtCuts[nptBins] = {0.36, 0.3, 0.29}; //original CM
float const bdtCuts[nptBins] = {0.46, 0.4, 0.39};
//float const bdtCuts[nptBins] = {-1, -1, -1};
const float momBins[nptBins+1] = {1,2,3,5};
TString ptbin[nptBins] = {"12", "23", "35"};

// _________________________________________________________
StPicoD0AnaMaker::StPicoD0AnaMaker(char const* name, StPicoDstMaker* picoMaker, char const* outputBaseFileName) :
        StPicoHFMaker(name, picoMaker, outputBaseFileName),
        mOutFileBaseName(outputBaseFileName), mSwitchRefit(false){
    // constructor
}

// _________________________________________________________
StPicoD0AnaMaker::~StPicoD0AnaMaker() {
    // destructor
}

// _________________________________________________________
int StPicoD0AnaMaker::InitHF() {
    // -- INITIALIZE USER HISTOGRAMS ETC HERE -------------------
    //    add them to the output list mOutList which is automatically written
    //
    // EXAMPLE //  mOutList->Add(new TH1F(...));
    // EXAMPLE //  TH1F* hist = static_cast<TH1F*>(mOutList->Last());
//    mOutList->Add(new TH2F("h_piTOF","h_piTOF",100,0,10, 250, -1, 1.5));
//    mOutList->Add(new TH2F("h_kTOF","h_kTOF",100,0,10, 250, -1, 1.5));
//    mOutList->Add(new TH2F("h_pTOF","h_pTOF",100,0,10, 250, -1, 1.5));
//
//    mOutList->Add(new TH2F("h_piTOF_20","h_piTOF_20",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_kTOF_20","h_kTOF_20",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_pTOF_20","h_pTOF_20",100,0,10, 300, 0, 1));
//
//    mOutList->Add(new TH2F("h_piTOF_HFT","h_piTOF_HFT",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_kTOF_HFT","h_kTOF_HFT",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_pTOF_HFT","h_pTOF_HFT",100,0,10, 300, 0, 1));
//
//    mOutList->Add(new TH2F("h_piTOF_HFT_20","h_piTOF_HFT_20",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_kTOF_HFT_20","h_kTOF_HFT_20",100,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_pTOF_HFT_20","h_pTOF_HFT_20",100,0,10, 300, 0, 1));
//
//    mOutList->Add(new TH2F("h_piTOFbeta","h_piTOFbeta",500,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_kTOFbeta","h_kTOFbeta",500,0,10, 300, 0, 1));
//    mOutList->Add(new TH2F("h_pTOFbeta","h_pTOFbeta",500,0,10, 300, 0, 1));
//
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Event Selection Cuts %%%%%%%%%%%%%%%%%%%%%%//

    /*
    mOutList->Add(new TH2F("hD0VsRemoved","hD0VsRemoved", 100, -0.5, 99.5, 100, -0.5, 99.5));

    mOutList->Add(new TH1F("hNTracksRemoved","hNTracksRemoved", 2000, -0.5, 1999.5));
    mOutList->Add(new TH1F("hNTracksPrimary","hNTracksPrimary", 2000, -0.5, 1999.5));
    mOutList->Add(new TH1F("hNTracksDiffRemovedPrimary","hNTracksDiffRemovedPrimary", 2000, -0.5, 1999.5));
    mOutList->Add(new TH1F("hNTracksDiffRemovedGlobal","hNTracksDiffRemovedGlobal", 5000, -0.5, 4999.5));
    mOutList->Add(new TH1F("hHotSpotDiffRemovedPrimary","hHotSpotDiffRemovedPrimary", 2000, -0.5, 1999.5));
    mOutList->Add(new TH1F("hNTracksGoodToFit","hNTracksGoodToFit", 2000, -0.5, 1999.5));
    */

    //mOutList->Add(new TH1F("hNTracksRemoved","hNTracksRemoved", 5000, 1000.5, 6000.5));
    //mOutList->Add(new TH1F("hNTracksPrimary","hNTracksPrimary", 5000, 1000.5, 6000.5));
    //mOutList->Add(new TH1F("hNTracksDiffRemovedPrimary","hNTracksDiffRemovedPrimary", 2000, -0.5, 1999.5));
    //mOutList->Add(new TH1F("hNTracksDiffRemovedGlobal","hNTracksDiffRemovedGlobal", 5000, -0.5, 4999.5));
    //mOutList->Add(new TH1F("hHotSpotDiffRemovedPrimary","hHotSpotDiffRemovedPrimary", 2000, -0.5, 1999.5));
    //mOutList->Add(new TH1F("hNTracksGoodToFit","hNTracksGoodToFit", 5000, 1000.5, 6000.5));
    
    mOutList->Add(new TH1F("hPrimVtZ", "hPrimVtZ before event selection; #it{v}_{z} [cm]; Entries", 500, -250, 250));
    mOutList->Add(new TH2F("hPrimVtXY", "hPrimVtXY before event selection; #it{v}_{x} [cm]; #it{v}_{y} [cm]", 250, -2.5, 2.5, 250, -2.5, 2.5));
    mOutList->Add(new TH2F("hVzVPDvsVzTPC", "Vz and VzVPD before event selection; #it{v}^{TPC}_{z} [cm]; #it{v}^{VPD}_{z} [cm]; Entries", 1500, -750, 750, 1500, -750, 750));
    mOutList->Add(new TH1F("hVzVPDvsVzTPCDiff", "Vz - VzVPD before event selection; #it{v}^{TPC}_{z} - #it{v}^{VPD}_{z} [cm]; Entries", 5000, -500, 500));
    /*//////////mOutList->Add(new TH1F("hPrimVtZ_beforematch", "hPrimVtZ before 2x matching in Fast", 500, -250, 250));
    mOutList->Add(new TH2F("hPrimVtXY_beforematch", "hPrimVtXY before 2x matching in Fast; #it{v}_{x}; #it{v}_{y}", 250, -2.5, 2.5, 250, -2.5, 2.5));
    mOutList->Add(new TH2F("hVzVPDvsVzTPC_beforematch", "Vz and VzVPD before 2x matching in Fast; #it{v}_{z}; #it{v}^{VPD}_{z}", 1500, -750, 750, 1500, -750, 750));
    mOutList->Add(new TH1F("hVzVPDvsVzTPCDiff_beforematch", "Vz - VzVPD before 2x matching in Fast", 5000, -500, 500));
    mOutList->Add(new TH1F("hPrimVtZ_checkvzvpdtpc", "hPrimVtZ to check vzVPD vzTPC correlation", 500, -250, 250));
    mOutList->Add(new TH2F("hPrimVtXY_checkvzvpdtpc", "hPrimVtXY to check vzVPD vzTPC correlation; #it{v}_{x}; #it{v}_{y}", 250, -2.5, 2.5, 250, -2.5, 2.5));
    mOutList->Add(new TH2F("hVzVPDvsVzTPC_checkvzvpdtpc", "Vz and VzVPD to check vzVPD vzTPC correlation; #it{v}_{z}; #it{v}^{VPD}_{z}", 1500, -750, 750, 1500, -750, 750));
    mOutList->Add(new TH1F("hVzVPDvsVzTPCDiff_checkvzvpdtpc", "Vz - VzVPD to check vzVPD vzTPC correlation", 5000, -500, 500));
    mOutList->Add(new TH1F("hPrimVtZ_wmatch", "hPrimVtZ after 2x matching in Fast", 500, -250, 250));
    mOutList->Add(new TH2F("hPrimVtXY_wmatch", "hPrimVtXY after 2x matching in Fast; #it{v}_{x}; #it{v}_{y}", 250, -2.5, 2.5, 250, -2.5, 2.5));
    mOutList->Add(new TH2F("hVzVPDvsVzTPC_wmatch", "Vz and VzVPD after 2x matching in Fast; #it{v}_{z}; #it{v}^{VPD}_{z}", 1500, -750, 750, 1500, -750, 750));
    mOutList->Add(new TH1F("hVzVPDvsVzTPCDiff_wmatch", "Vz - VzVPD after 2x matching in Fast", 5000, -500, 500));
    *///////////
    mOutList->Add(new TH1F("hPrimVtZ_evcut", "hPrimVtZ after event selection; #it{v}_{z} [cm]; Entries", 500, -250, 250));
    mOutList->Add(new TH2F("hPrimVtXY_evcut", "hPrimVtXY after event selection; #it{v}_{x} [cm]; #it{v}_{y} [cm]", 250, -2.5, 2.5, 250, -2.5, 2.5));
    mOutList->Add(new TH2F("hVzVPDvsVzTPC_evcut", "Vz and VzVPD after event selection; #it{v}^{TPC}_{z} [cm]; #it{v}^{VPD}_{z} [cm]", 1500, -750, 750, 1500, -750, 750));
    mOutList->Add(new TH1F("hVzVPDvsVzTPCDiff_evcut", "Vz - VzVPD before event selection; #it{v}^{TPC}_{z} - #it{v}^{VPD}_{z} [cm]; Entries", 5000, -500, 500));
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QA histograms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    //////////mOutList->Add(new TH1F("h_RefMult", "RefMult distribution after event cuts; RefMult; Events", 1000, 0, 1000));
    //////////mOutList->Add(new TH1F("h_gRefMult", "gRrefMult distribution after event cuts; gRefMult; Events", 1000, 0, 1000));
    mOutList->Add(new TH2F("h_gRefMultVsRefMult","gRefMult vs. RefMult before event selection;RefMult;gRefMult",100,0,100,100,0,100));
    mOutList->Add(new TH2F("h_gRefMultVsRefMult_evcut","gRefMult vs. RefMult after event selection;RefMult;gRefMult",100,0,100,100,0,100));
    mOutList->Add(new TH2F("hTofMultVsRefMult","TofMult vs. RefMult before event selection;RefMult;TofMult",100,0,100,100,0,100));
    mOutList->Add(new TH2F("hTofMultVsRefMult_evcut","TofMult vs. RefMult after event selection;RefMult;TofMult",100,0,100,100,0,100));
   
   mOutList->Add(new TProfile("prof_nPrimaryTracks_vs_BBCx", "nPrimaryTracks vs BBCx;BBC Coincidence Rate / 1000; primary track multiplicity", 5000, 1000.5, 6000.5));
   mOutList->Add(new TProfile("prof_nGoodTracks_vs_BBCx", "nGoodTracks vs BBCx;BBC Coincidence Rate / 1000; good track multiplicity", 5000, 1000.5, 6000.5));
   mOutList->Add(new TProfile("prof_ntracksPionTPC_vs_BBCx", "nTPCpionTracks vs BBCx;BBC Coincidence Rate / 1000; TPCpion track multiplicity", 5000, 1000.5, 6000.5));
   mOutList->Add(new TProfile("prof_ntracksKaonTPC_vs_BBCx", "nTPCkaonTracks vs BBCx;BBC Coincidence Rate / 1000; TPCkaon track multiplicity", 5000, 1000.5, 6000.5));
   mOutList->Add(new TProfile("prof_ntracksPion_vs_BBCx", "npionTracks vs BBCx;BBC Coincidence Rate / 1000; pion track multiplicity", 5000, 1000.5, 6000.5));
   mOutList->Add(new TProfile("prof_ntracksKaon_vs_BBCx", "nkaonTracks vs BBCx;BBC Coincidence Rate / 1000; kaon track multiplicity", 5000, 1000.5, 6000.5));

   mOutList->Add(new TProfile("prof_nPrimaryTracks_vs_runID", "nPrimaryTracks vs runID (converted);runID (converted); primary track multiplicity", 4000, -499.5, 3500.5));
   mOutList->Add(new TProfile("prof_nGoodTracks_vs_runID", "nGoodTracks vs runID (converted);runID (converted); good track multiplicity", 4000, -499.5, 3500.5));
   mOutList->Add(new TProfile("prof_ntracksKaonTPC_vs_runID", "nTPCkaonTracks vs runID (converted); runID (converted); TPCkaon track multiplicity", 4000, -499.5, 3500.5));
   mOutList->Add(new TProfile("prof_ntracksPionTPC_vs_runID", "nTPCpionTracks vs runID (converted); runID (converted); TPCpion track multiplicity", 4000, -499.5, 3500.5));
   mOutList->Add(new TProfile("prof_ntracksKaon_vs_runID", "nkaonTracks vs runID (converted);runID (converted); kaon track multiplicity", 4000, -499.5, 3500.5));
   mOutList->Add(new TProfile("prof_ntracksPion_vs_runID", "npionTracks vs runID (converted);runID (converted); pion track multiplicity", 4000, -499.5, 3500.5));
    
    //////////mOutList->Add(new TH2F("h_nPrimaryTracks_vs_BBCx", "nPrimaryTracks vs BBCx;BBC Coincidence Rate / 1000; primary track multiplicity", 5000, 1000.5, 6000.5, 200, 0, 200));
    //////////mOutList->Add(new TH2F("h_nGoodTracks_vs_BBCx", "nGoodTracks vs BBCx;BBC Coincidence Rate / 1000; good track multiplicity", 5000, 1000.5, 6000.5, 200, 0, 200));
     //////////mOutList->Add(new TH2F("h_ntracksKaonTPC_vs_BBCx", "nTPCkaonTracks vs BBCx;BBC Coincidence Rate / 1000; TPCkaon track multiplicity", 5000, 1000.5, 6000.5, 200, 0, 200));
     //////////mOutList->Add(new TH2F("h_ntracksPionTPC_vs_BBCx", "nTPCpionTracks vs BBCx;BBC Coincidence Rate / 1000; TPCpion track multiplicity", 5000, 1000.5, 6000.5, 200, 0, 200));
    //////////mOutList->Add(new TH2F("h_ntracksKaon_vs_BBCx", "nkaonTracks vs BBCx;BBC Coincidence Rate / 1000; kaon track multiplicity", 5000, 1000.5, 6000.5, 200, 0, 200));
    //////////mOutList->Add(new TH2F("h_ntracksPion_vs_BBCx", "npionTracks vs BBCx;BBC Coincidence Rate / 1000; pion track multiplicity", 5000, 1000.5, 6000.5, 200, 0, 200));

    //////////mOutList->Add(new TH2F("h_ntracksKaonTPC_vs_runID", "nTPCkaonTracks vs runID (converted); runID (converted); TPCkaon track multiplicity", 4000, -499.5, 3500.5, 200, 0, 200));
    //////////mOutList->Add(new TH2F("h_ntracksPionTPC_vs_runID", "nTPCpionTracks vs runID (converted); runID (converted); TPCpion track multiplicity", 4000, -499.5, 3500.5, 200, 0, 200));


    
    /*mOutList->Add(new TH2F("h_nPrimaryTracks_vs_runID", "nPrimaryTracks vs runID (converted);runID (converted); primary track multiplicity", 4000, -499.5, 3500.5, 200, 0, 200));
    mOutList->Add(new TH2F("h_nGoodTracks_vs_runID", "nGoodTracks vs runID (converted);runID (converted); good track multiplicity", 4000, -499.5, 3500.5, 200, 0, 200));
    mOutList->Add(new TH2F("h_ntracksKaonTPC_vs_runID", "nTPCkaonTracks vs runID (converted); runID (converted); TPCkaon track multiplicity", 4000, -499.5, 3500.5, 200, 0, 200));
    mOutList->Add(new TH2F("h_ntracksPionTPC_vs_runID", "nTPCpionTracks vs runID (converted); runID (converted); TPCpion track multiplicity", 4000, -499.5, 3500.5, 200, 0, 200));
    mOutList->Add(new TH2F("h_ntracksKaon_vs_runID", "nkaonTracks vs runID (converted);runID (converted); kaon track multiplicity", 4000, -499.5, 3500.5, 200, 0, 200));
    mOutList->Add(new TH2F("h_ntracksPion_vs_runID", "npionTracks vs runID (converted);runID (converted); pion track multiplicity", 4000, -499.5, 3500.5, 200, 0, 200));

    mOutList->Add(new TH2F("h_nPrimaryTrackswMatching_vs_BBCx", "nPrimaryTracks with matching vs BBCx;BBC Coincidence Rate / 1000; primary track multiplicity", 5000, 1000.5, 6000.5, 200, 0, 200));
    mOutList->Add(new TH2F("h_nGoodTrackswMatching_vs_BBCx", "nGoodTracks with matching vs BBCx;BBC Coincidence Rate / 1000; good track multiplicity", 5000, 1000.5, 6000.5, 200, 0, 200));
    mOutList->Add(new TH2F("h_ntracksKaonTPCwMatching_vs_BBCx", "nTPCkaonTracks with matching vs BBCx;BBC Coincidence Rate / 1000; TPCkaon track multiplicity", 5000, 1000.5, 6000.5, 200, 0, 200));
    mOutList->Add(new TH2F("h_ntracksPionTPCwMatching_vs_BBCx", "nTPCpionTracks with matching vs BBCx;BBC Coincidence Rate / 1000; TPCpion track multiplicity", 5000, 1000.5, 6000.5, 200, 0, 200));
    
    
    mOutList->Add(new TH2F("h_nPrimaryTrackswMatching_vs_runID", "nPrimaryTracks with matching vs runID (converted);runID (converted); primary track multiplicity", 4000, -499.5, 3500.5, 200, 0, 200));
    mOutList->Add(new TH2F("h_nGoodTrackswMatching_vs_runID", "nGoodTracks with matching vs runID (converted);runID (converted); good track multiplicity", 4000, -499.5, 3500.5, 200, 0, 200));
    mOutList->Add(new TH2F("h_ntracksKaonTPCwMatching_vs_runID", "nTPCkaonTracks with matching vs runID (converted); runID (converted); TPCkaon track multiplicity", 4000, -499.5, 3500.5, 200, 0, 200));
    mOutList->Add(new TH2F("h_ntracksPionTPCwMatching_vs_runID", "nTPCpionTracks with matching vs runID (converted); runID (converted); TPCpion track multiplicity", 4000, -499.5, 3500.5, 200, 0, 200));*/
    
    
 
    /*mOutList->Add(new TH2F("h_nPrimaryTracks_vs_runID", "nPrimaryTracks vs runID;runID; primary track multiplicity", 100001, 49999.5, 150000.5, 200, 0, 200));
    
    //////////mOutList->Add(new TProfile("prof_gRefmult_vs_BBCx", "gRefmult vs BBCx;BBC Coincidence Rate / 1000; gRefmult", 5000, 1000.5, 6000.5));
    */
    
    
    
    
   
    //////////mOutList->Add(new TH2F("h_gRefmult_vs_BBCx", "gRefmult vs BBCx;BBC Coincidence Rate / 1000; gRefmult", 5000, 1000.5, 6000.5, 200, 0, 200));
    
    
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Track Quality Cuts %%%%%%%%%%%%%%%%%%%%%%%%%//

    // ---------------------------Track cuts check--------------------------------------------------------
    mOutList->Add(new TH1F("hpT_tr", "track pT before cut; p_{T} [GeV/c]; Entries", 500, 0, 50));
    mOutList->Add(new TH1F("hpT_tr_cut", "track pT after cut; p_{T} [GeV/c]; Entries", 500, 0, 50));
    mOutList->Add(new TH1F("hEta_tr", "track #eta before cut;#eta; Entries", 100, -5, 5));
    mOutList->Add(new TH1F("hEta_tr_cut", "track #eta after cut;#eta; Entries", 100, -5, 5));
    //////////mOutList->Add(new TH1F("hPhi_tr", "track #varphi before cut;#phi", 120, -TMath::TwoPi(), TMath::TwoPi()));
    //////////mOutList->Add(new TH1F("hPhi_tr_cut", "track #varphi after cut;#phi", 120, -TMath::TwoPi(), TMath::TwoPi()));
    mOutList->Add(new TH2F("hEtaPhi_tr", "#varphi vs #eta distribution before cut; #varphi; #eta", 120, -TMath::TwoPi(), TMath::TwoPi(), 100, -5, 5));
    mOutList->Add(new TH2F("hEtaPhi_tr_cut", "#varphi vs #eta distribution after cut; #varphi; #eta", 120, -TMath::TwoPi(), TMath::TwoPi(), 100, -5, 5));
    mOutList->Add(new TH1F("hNHitsFit", "nHitsFit before cut;nHitsFit; Entries", 60, 0, 60));
    mOutList->Add(new TH1F("hNHitsFit_cut", "nHitsFit after cut;nHitsFit; Entries", 60, 0, 60));
    //////////mOutList->Add(new TH1F("hNHitsMax", "nHitsMax before cut;nHitsMax", 50, 0, 50));
    mOutList->Add(new TH1F("hNHitsFitnHitsMax", "hNHitsFitnHitsMax before cut;NHitsFit/NHitsMax; Entries", 30, 0, 1.5));
    mOutList->Add(new TH1F("hNHitsFitnHitsMax_cut", "hNHitsFitnHitsMax after cut;NHitsFit/NHitsMax; Entries", 30, 0, 1.5));
    mOutList->Add(new TH1F("hDca", "Dca before cut;DCA [cm]; Entries", 1000, -1, 50));
    mOutList->Add(new TH1F("hDca_cut", "Dca after cut;DCA [cm]; Entries", 1000, -1, 50));
    mOutList->Add(new TH2F("hDcaXY", "DcaXY before cut;DCAx [cm]; DCAy [cm]", 150, -7.5, 7.5, 150, -7.5, 7.5));
    mOutList->Add(new TH2F("hDcaXY_cut", "DcaXY after cut;DCAx [cm]; DCAy [cm]", 150, -7.5, 7.5, 150, -7.5, 7.5));
    mOutList->Add(new TH1F("hDcaZ", "DcaZ before cut; DCAz [cm]; Entries", 1000, -50, 50));
    mOutList->Add(new TH1F("hDcaZ_cut", "DcaZ after cut;DCAz [cm]; Entries", 1000, -50, 50));
    
    //------------------------ Energy Loss Check----------------------------------------//

    mOutList->Add(new TH2F("hDedx_tr", "dE/dx after track cut;p*charge[GeV/c];dE/dx [keV/cm]", 1000, -10, 10, 5000, 0, 75));
    /////mOutList->Add(new TH2D("h_nSigmadEdXPion", "n#sigma_{#pi} vs. p;p;n#sigma_{#pi}", 1000, -10, 10, 5000, -50, 50));
	/////mOutList->Add(new TH2D("h_nSigmadEdXKaon", "n#sigma_{K} vs. p;p;n#sigma_{K}", 1000, -10, 10, 5000, -50, 50));
    mOutList->Add(new TH2F("hPiDedx_TOF", "Pion dE/dx after TOF cut;p*charge[GeV/c];dE/dx [keV/cm]", 500, -5, 5, 3000, 0, 30));
    mOutList->Add(new TH2F("hPiDedx_TOF_TPC", "Pion dE/dx after TOF and TPC cut;p*charge[GeV/c];dE/dx [keV/cm]", 500, -5, 5, 3000, 0, 30));
    mOutList->Add(new TH2F("hKDedx_TOF", "Kaon dE/dx after TOF cut;p*charge[GeV/c];dE/dx [keV/cm]", 500, -5, 5, 3000, 0, 30));
    mOutList->Add(new TH2F("hKDedx_TOF_TPC", "Kaon dE/dx after TOF and TPC cut;p*charge[GeV/c];dE/dx [keV/cm]", 500, -5, 5, 3000, 0, 30));
    
    //----------------------- 1/Beta Check ------------------------------------------- //

    mOutList->Add(new TH2F("hBetavsP_tr", "1/#beta vs p after track cuts;p[GeV/c];1/#beta", 250, 0, 10, 1000, 0.0, 10.0));
    mOutList->Add(new TH2F("hPiBetavsP_TPC", "1/#beta vs p for pions after TPC cut;p[GeV/c];1/#beta",  250, 0, 10, 1000, 0.0, 10.0));
    mOutList->Add(new TH2F("hPiBetavsP_TOF_TPC", "1/#beta vs p for pions after TOF&TPC cut;p[GeV/c];1/#beta", 250, 0, 10, 1000, 0.0, 10.0));
    mOutList->Add(new TH2F("hKBetavsP_TPC", "1/#beta vs p for kaons after TPC cut;p[GeV/c];1/#beta",  250, 0, 10, 1000, 0.0, 10.0));
    mOutList->Add(new TH2F("hKBetavsP_TOF_TPC", "1/#beta vs p for kions after TOF&TPC cut;p[GeV/c];1/#beta", 250, 0, 10, 1000, 0.0, 10.0));
    

   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PID cuts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

    mOutList->Add(new TH2F("h_nSigmaOneOverBetaPi_tr", "n#sigma^{#pi}_{1/#beta} after track cut; p[GeV/c]; n#sigma^{#pi}_{TOF}", 250, 0,10,400,-40,40));
    mOutList->Add(new TH2F("h_nSigmaOneOverBetaPi_TPC", "n#sigma^{#pi}_{1/#beta} after TPC cut; p[GeV/c]; n#sigma^{#pi}_{TOF}", 250, 0,10,400,-40,40));
    mOutList->Add(new TH2F("h_nSigmaOneOverBetaPi_TOF_TPC", "n#sigma^{#pi}_{1/#beta} after TOF and TPC cut;p[GeV/c] ;n#sigma^{#pi}_{TOF}", 250, 0,10,400,-40,40));
    mOutList->Add(new TH2F("h_nSigmaOneOverBetaK_tr", "n#sigma^{K}_{1/#beta} after track cut; p[GeV/c]; n#sigma^{K}_{TOF}", 250, 0,10,400,-40,40));
    mOutList->Add(new TH2F("h_nSigmaOneOverBetaK_TPC", "n#sigma^{K}_{1/#beta} after TPC cut; p[GeV/c]; n#sigma^{K}_{TOF}", 250, 0,10,400,-40,40));
    mOutList->Add(new TH2F("h_nSigmaOneOverBetaK_TOF_TPC", "n#sigma^{K}_{1/#beta} after TOF and TPC cut; p[GeV/c]; n#sigma^{K}_{TOF}", 250, 0,10,400,-40,40));

    mOutList->Add(new TH2F("h_DeltaOneOverBetaPi_tr", "#Delta^{#pi}_{1/#beta} after track cut; p[GeV/c]; #Delta^{#pi}_{1/#beta}", 250, 0,10,6000,-5,55));
    mOutList->Add(new TH2F("h_DeltaOneOverBetaPi_TPC", "#Delta^{#pi}_{1/#beta} after TPC cut; p[GeV/c]; #Delta^{#pi}_{1/#beta}", 250, 0,10,6000,-5,55));
    mOutList->Add(new TH2F("h_DeltaOneOverBetaPi_TOF_TPC", "#Delta^{#pi}_{1/#beta} after TOF and TPC cut;p[GeV/c] ;#Delta^{#pi}_{1/#beta}", 250, 0,10,6000,-5,55));
    mOutList->Add(new TH2F("h_DeltaOneOverBetaK_tr", "#Delta^{K}_{1/#beta} after track cut; p[GeV/c]; #Delta^{K}_{1/#beta}", 250, 0,10,6000,-5,55));
    mOutList->Add(new TH2F("h_DeltaOneOverBetaK_TPC", "#Delta^{K}_{1/#beta} after TPC cut; p[GeV/c]; #Delta^{K}_{1/#beta}", 250, 0,10,6000,-5,55));
    mOutList->Add(new TH2F("h_DeltaOneOverBetaK_TOF_TPC", "#Delta^{K}_{1/#beta} after TOF and TPC cut; p[GeV/c]; #Delta^{K}_{1/#beta}", 250, 0,10,6000,-5,55));
    
    
    mOutList->Add(new TH2F("h_nSigmadEdxPi_tr", "n#sigma^{#pi}_{dE/dx} after track cut; p [GeV/c]; n#sigma^{#pi}_{TPC}", 250, 0,10,400,-40,40));
    mOutList->Add(new TH2F("h_nSigmadEdxPi_TOF", "n#sigma^{#pi}_{dE/dx} after TOF cut; p [GeV/c]; n#sigma^{#pi}_{TPC}", 250, 0,10,400,-40,40));
    mOutList->Add(new TH2F("h_nSigmadEdxPi_TOF_TPC", "n#sigma^{#pi}_{dE/dx} after TOF and TPC cut; p [GeV/c]; n#sigma^{#pi}_{TPC}", 250, 0,10,400,-40,40));
    mOutList->Add(new TH2F("h_nSigmadEdxK_tr", "n#sigma^{K}_{dE/dx} after track cut; p [GeV/c]; n#sigma^{K}_{TPC}", 250, 0,10,400,-40,40));
    mOutList->Add(new TH2F("h_nSigmadEdxK_TOF", "n#sigma^{K}_{dE/dx} after TOF cut; p [GeV/c]; n#sigma^{K}_{TPC}", 250, 0,10,400,-40,40));
    mOutList->Add(new TH2F("h_nSigmadEdxK_TOF_TPC", "n#sigma^{K}_{dE/dx} after TOF and TPC cut; p [GeV/c]; n#sigma^{K}_{TPC}", 250, 0,10,400,-40,40));
   
    mOutList->Add(new TH1F("hPionPt","hPionPt after all cuts; p_{T} [GeV]; Entries", 200, 0, 20));
    mOutList->Add(new TH1F("hKaonPt","hKaonPt after all cuts; p_{T} [GeV]; Entries", 200, 0, 20));

//    mOutList->Add(new TH2F("h_pnsigma","h_pnsigma",1000,0,10, 99, -5, 5));

//    mOutList->Add(new TH2F("h_dedx","h_dedx", 1000, 0, 10, 1000, 0, 10));
//h_tracktest
   //////mOutList->Add(new TH1D("h_tracktest","h_tracktest", 6, 0.5, 6.5));
    //////mOutList->Add(new TH1D("h_tracktest_TOF","h_tracktest_TOF", 6, 0.5, 6.5));


    mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
    TString ntpVars = /*"grefMult:eventId:ZDC:BBC:hotSpot:primary:diffRemovedPrimary:*/"runId:pi1_pt:pi1_p:pi1_nSigmaTPC:pi1_nSigmaTOF:pi1_TOFDeltainvbeta:pi1_charge:pi1_isTOFmatched:pi1_isBEMCmatched:k_pt:k_p:k_nSigmaTPC:k_nSigmaTOF:k_TOFDeltainvbeta:k_charge:k_isTOFmatched:k_isBEMCmatched:dcaDaughters:D_theta:D_decayL:dcaD0ToPv:D_cosThetaStar:D_pt:D_mass:D_rapidity:D_phi:D_eta:refMult:k_nHitsFit:pi1_nHitsFit:k_NHitsDedx:pi1_NHitsDedx:k_btofYLocal:pi1_btofYLocal:VzTPC:VzVPD:k_gDCA:pi1_gDCA";
    TString ntpVars2 = /*"grefMult:eventId:ZDC:BBC:hotSpot:primary:diffRemovedPrimary*/"runId:pi1_pt:pi1_p:pi1_nSigmaTPC:pi1_nSigmaTOF:pi1_TOFDeltainvbeta:pi1_charge:pi1_isTOFmatched:pi1_isBEMCmatched:k_pt:k_p:k_nSigmaTPC:k_nSigmaTOF:k_TOFDeltainvbeta:k_charge:k_isTOFmatched:k_isBEMCmatched:dcaDaughters:D_theta:D_decayL:dcaD0ToPv:D_cosThetaStar:D_pt:D_mass:D_rapidity:D_phi:D_eta:refMult:k_nHitsFit:pi1_nHitsFit:k_NHitsDedx:pi1_NHitsDedx:k_btofYLocal:pi1_btofYLocal:VzTPC:VzVPD:k_gDCA:pi1_gDCA";
    ///////////TString ntpVars2 = "grefMult:refMult:runId:eventId:ZDC:BBC:hotSpot:primary:diffRemovedPrimary:pi1_pt:pi1_p:pi1_dca:pi1_nSigma:pi1_nHits:pi1_nHitFit:pi1_TOFinvbeta:pi1_betaBase:pi1_charge:k_pt:k_p:k_dca:k_nSigma:k_nHits:k_nHitFit:k_TOFinvbeta:k_betaBase:k_charge:dcaDaughters:primVz:primVzVpd:D_theta:cosTheta:D_decayL:dcaD0ToPv:D_cosThetaStar:D_pt:D_mass:D_rapidity:VzVPD_VzTPC";
    TString ntpVars3 = /*"grefMult:refMult:eventId:ZDC:BBC:hotSpot:primary:diffRemovedPrimary:*/"runId:pi1_pt:pi1_p:pi1_nSigmaTPC:pi1_nSigmaTOF:pi1_TOFDeltainvbeta:pi1_charge:pi1_isTOFmatched:pi1_isBEMCmatched:k_pt:k_p:k_nSigmaTPC:k_nSigmaTOF:k_TOFDeltainvbeta:k_charge:k_isTOFmatched:k_isBEMCmatched:pi2_pt:pi2_p:pi2_nSigmaTPC:pi2_nSigmaTOF:pi2_TOFDeltainvbeta:pi2_charge:pi2_isTOFmatched:pi2_isBEMCmatched:D_mass:D_pt:D_rapidity:D_cosThetaStar:D_theta:Dstar_pt:Dstar_mass:Dstar_rapidity:refMult:k_nHitsFit:pi1_nHitsFit:pi2_nHitsFit:k_NHitsDedx:pi1_NHitsDedx:pi2_NHitsDedx:k_btofYLocal:pi1_btofYLocal:pi2_btofYLocal:VzTPC:VzVPD:k_gDCA:pi1_gDCA:pi2_gDCA";

//    ntp_kaon = new TNtuple("ntp_kaon", "kaon tree","k_pt:k_phi:k_eta:k_nSigma:k_nHitFit:k_TOFinvbeta:pi_eventId:pi_runId");
//    ntp = new TNtuple("ntp", "pion tree","pi_pt:pi_phi:pi_eta:pi_nSigma:pi_nHitFit:pi_TOFinvbeta:k_eventId:k_runId");
    ntp_DMeson_UnlikeSign = new TNtuple("ntp_UnlikeSign","DMeson TreeUnlikeSign", ntpVars);
    ntp_DMeson_LikeSign = new TNtuple("ntp_LikeSign","DMeson TreeLikeSign", ntpVars);
    ntp_DMeson_Rotated = new TNtuple("ntp_Rotated","DMeson TreeRotated", ntpVars2);
    
    
    ntp_DstarMeson_RightSign = new TNtuple("Dntp_RightSign","DstarMeson TreeRightSign", ntpVars3);
    ntp_DstarMeson_WrongSign = new TNtuple("Dntp_Wrongsign","DstarMeson TreeWrongsign", ntpVars3);
    ntp_DstarMeson_SideBand = new TNtuple("Dntp_SideBand","DstarMeson TreeSideband", ntpVars3);
    

    if (mSwitchRefit) {
        TString dir = "./StRoot/weights/";
        TString prefix = "TMVAClassification";

        for (int pT = 0; pT < nptBins; pT++) {
            reader[pT] = new TMVA::Reader("!Color:!Silent");
            reader[pT]->AddVariable("k_dca", &k_dca[pT]);
            reader[pT]->AddVariable("pi1_dca", &pi1_dca[pT]);
            reader[pT]->AddVariable("dcaDaughters", &dcaDaughters[pT]);
            reader[pT]->AddVariable("cosTheta", &cosTheta[pT]);
            reader[pT]->AddVariable("D_decayL", &D_decayL[pT]);
            reader[pT]->AddVariable("dcaD0ToPv", &dcaD0ToPv[pT]);
//            reader[pT]->AddVariable("D_cosThetaStar", &thetaStar[pT]);

            TString methodName = "BDT method";
            TString weightfile = dir + prefix + TString("_BDT.weights.pt") + ptbin[pT] + TString(".xml");
            reader[pT]->BookMVA(methodName, weightfile);
        }
    }

    return kStOK;
}

// _________________________________________________________
void StPicoD0AnaMaker::ClearHF(Option_t *opt="") {
    return;
}

// _________________________________________________________
int StPicoD0AnaMaker::FinishHF() {
    ntp_DMeson_UnlikeSign -> Write(ntp_DMeson_UnlikeSign->GetName(), TObject::kOverwrite);
    ntp_DMeson_Rotated -> Write(ntp_DMeson_Rotated->GetName(), TObject::kOverwrite);
    ntp_DMeson_LikeSign -> Write(ntp_DMeson_LikeSign->GetName(), TObject::kOverwrite);
    ntp_DstarMeson_RightSign -> Write(ntp_DstarMeson_RightSign->GetName(), TObject::kOverwrite);
    ntp_DstarMeson_WrongSign -> Write(ntp_DstarMeson_WrongSign->GetName(), TObject::kOverwrite);
    ntp_DstarMeson_SideBand -> Write(ntp_DstarMeson_SideBand->GetName(), TObject::kOverwrite);
    
//    ntp -> Write(ntp->GetName(), TObject::kOverwrite);
//    ntp_kaon -> Write(ntp_kaon->GetName(), TObject::kOverwrite);
    return kStOK;
}
// _________________________________________________________
int StPicoD0AnaMaker::MakeHF() {
    createCandidates();
//    analyzeCandidates();

//    TH2F *h_pTOF = static_cast<TH2F*>(mOutList->FindObject("h_pTOF"));
//
//    TH2F *h_piTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_20"));
//    TH2F *h_kTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_20"));
//    TH2F *h_pTOF_20 = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_20"));
//
//    TH2F *h_piTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_HFT"));
//    TH2F *h_kTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_HFT"));
//    TH2F *h_pTOF_HFT = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_HFT"));
//
//    TH2F *h_piTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_piTOF_HFT_20"));
//    TH2F *h_kTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_kTOF_HFT_20"));
//    TH2F *h_pTOF_HFT_20 = static_cast<TH2F*>(mOutList->FindObject("h_pTOF_HFT_20"));
//
//    TH2F *h_piTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_piTOFbeta"));
//    TH2F *h_kTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_kTOFbeta"));
//    TH2F *h_pTOFbeta = static_cast<TH2F*>(mOutList->FindObject("h_pTOFbeta"));
//
//    TH2F *h_pinsigma = static_cast<TH2F*>(mOutList->FindObject("h_pinsigma"));
//    TH2F *h_knsigma = static_cast<TH2F*>(mOutList->FindObject("h_knsigma"));
//    TH2F *h_pnsigma = static_cast<TH2F*>(mOutList->FindObject("h_pnsigma"));
//
//    TH2F *h_dedx = static_cast<TH2F*>(mOutList->FindObject("h_dedx"));
//
//    TVector3 pVtx = mPicoDst->event()->primaryVertex();
//
//    UInt_t nTracks = mPicoDst->numberOfTracks();
//    for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack){
//        StPicoTrack const* trk = mPicoDst->track(iTrack);
//        if (!trk) continue;
//        StPhysicalHelixD helix = trk->helix(mBField);
//        TVector3 momentum = trk->pMom(pVtx, mPicoDst->event()->bField());
//
//        if (!(trk->nHitsFit()>=15)) continue;
//        if (!(fabs(momentum.pseudoRapidity()) <= 1.0)) continue;
//
//        if (fabs(trk->nSigmaPion())<3.0){
//            if (mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion)) {
//                    h_piTOF->Fill(trk->pPt(),mHFCuts->getTofBetaBase(trk));
//                    float oneOverBeta = getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kPion);
//                    h_piTOFbeta->Fill(trk->pPt(),oneOverBeta);
//                    if (trk->nHitsFit()>=20) h_piTOF_20->Fill(trk->pPt(),oneOverBeta);
//                    if (trk->isHFTTrack()) h_piTOF_HFT->Fill(trk->pPt(),oneOverBeta);
//                    if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_piTOF_HFT_20->Fill(trk->pPt(),oneOverBeta);
//            }
//        }
//
//        if (fabs(trk->nSigmaKaon())<3.0){
//            if (mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon)) {
//                h_kTOF->Fill(trk->pPt(),mHFCuts->getTofBetaBase(trk));
//                float oneOverBeta = getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kKaon);
//                h_kTOFbeta->Fill(trk->pPt(),oneOverBeta);
//                if (trk->nHitsFit()>=20) h_kTOF_20->Fill(trk->pPt(),oneOverBeta);
//                if (trk->isHFTTrack()) h_kTOF_HFT->Fill(trk->pPt(),oneOverBeta);
//                if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_kTOF_HFT_20->Fill(trk->pPt(),oneOverBeta);
//            }
//        }
//
//        if (fabs(trk->nSigmaProton())<3.0){
//            if (mHFCuts->isTOFHadron(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kProton)) {
//                h_pTOF->Fill(trk->pPt(),mHFCuts->getTofBetaBase(trk));
//                float oneOverBeta = getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StHFCuts::kProton);
//                h_pTOFbeta->Fill(trk->pPt(),oneOverBeta);
//                if (trk->nHitsFit()>=20) h_pTOF_20->Fill(trk->pPt(),oneOverBeta);
//                if (trk->isHFTTrack()) h_pTOF_HFT->Fill(trk->pPt(),oneOverBeta);
//                if ((trk->isHFTTrack()) && (trk->nHitsFit()>=20)) h_pTOF_HFT_20->Fill(trk->pPt(),oneOverBeta);
//            }
//        }
//
//        h_pinsigma->Fill(momentum.Mag(),trk->nSigmaPion());
//        h_knsigma->Fill(momentum.Mag(),trk->nSigmaKaon());
//        h_pnsigma->Fill(momentum.Mag(),trk->nSigmaProton());
//        h_dedx->Fill(momentum.Mag(),trk->dEdx());
//
//    } // .. end tracks loop

    return kStOK;
}

// _________________________________________________________
int StPicoD0AnaMaker::createCandidates() {
    //make it run faster:
 //   if (!(mPicoEvent->BBCx()<950000)) return kStOK;
    nGoodTracks=0;
//    if (!(abs(mPrimVtx.x())<0.6)) return kStOK;
//    if (!(abs(mPrimVtx.y())<0.6)) return kStOK;

 /*/////
    TH2F *hD0VsRemoved = static_cast<TH2F*>(mOutList->FindObject("hD0VsRemoved"));
    TH1F *hNTracksDiffRemovedPrimary = static_cast<TH1F*>(mOutList->FindObject("hNTracksDiffRemovedPrimary"));
    TH1F *hNTracksDiffRemovedGlobal = static_cast<TH1F*>(mOutList->FindObject("hNTracksDiffRemovedGlobal"));
    TH1F *hHotSpotDiffRemovedPrimary = static_cast<TH1F*>(mOutList->FindObject("hHotSpotDiffRemovedPrimary"));
    TH1F *hNTracksGoodToFit = static_cast<TH1F*>(mOutList->FindObject("hNTracksGoodToFit"));
    *//////
    
    //TH1F *hNTracksRemoved = static_cast<TH1F*>(mOutList->FindObject("hNTracksRemoved"));
    //TH1F *hNTracksPrimary = static_cast<TH1F*>(mOutList->FindObject("hNTracksPrimary"));
    //TH1F *hNTracksGoodToFit = static_cast<TH1F*>(mOutList->FindObject("hNTracksGoodToFit"));
    
    TH1F *hPrimVtZ = static_cast<TH1F*>(mOutList->FindObject("hPrimVtZ"));
    TH2F *hPrimVtXY = static_cast<TH2F*>(mOutList->FindObject("hPrimVtXY"));
    TH2F *hVzVPDvsVzTPC = static_cast<TH2F*>(mOutList->FindObject("hVzVPDvsVzTPC"));
    TH1F *hVzVPDvsVzTPCDiff = static_cast<TH1F*>(mOutList->FindObject("hVzVPDvsVzTPCDiff"));
    /*//////////TH1F *hPrimVtZ_beforematch = static_cast<TH1F*>(mOutList->FindObject("hPrimVtZ_beforematch"));
    TH2F *hPrimVtXY_beforematch = static_cast<TH2F*>(mOutList->FindObject("hPrimVtXY_beforematch"));
    TH2F *hVzVPDvsVzTPC_beforematch = static_cast<TH2F*>(mOutList->FindObject("hVzVPDvsVzTPC_beforematch"));
    TH1F *hVzVPDvsVzTPCDiff_beforematch = static_cast<TH1F*>(mOutList->FindObject("hVzVPDvsVzTPCDiff_beforematch"));
    TH1F *hPrimVtZ_checkvzvpdtpc = static_cast<TH1F*>(mOutList->FindObject("hPrimVtZ_checkvzvpdtpc"));
    TH2F *hPrimVtXY_checkvzvpdtpc = static_cast<TH2F*>(mOutList->FindObject("hPrimVtXY_checkvzvpdtpc"));
    TH2F *hVzVPDvsVzTPC_checkvzvpdtpc = static_cast<TH2F*>(mOutList->FindObject("hVzVPDvsVzTPC_checkvzvpdtpc"));
    TH1F *hVzVPDvsVzTPCDiff_checkvzvpdtpc = static_cast<TH1F*>(mOutList->FindObject("hVzVPDvsVzTPCDiff_checkvzvpdtpc"));
    TH1F *hPrimVtZ_wmatch = static_cast<TH1F*>(mOutList->FindObject("hPrimVtZ_wmatch"));
    TH2F *hPrimVtXY_wmatch = static_cast<TH2F*>(mOutList->FindObject("hPrimVtXY_wmatch"));
    TH2F *hVzVPDvsVzTPC_wmatch = static_cast<TH2F*>(mOutList->FindObject("hVzVPDvsVzTPC_wmatch"));
    TH1F *hVzVPDvsVzTPCDiff_wmatch = static_cast<TH1F*>(mOutList->FindObject("hVzVPDvsVzTPCDiff_wmatch"));
    *///////////
    TH1F *hPrimVtZ_evcut = static_cast<TH1F*>(mOutList->FindObject("hPrimVtZ_evcut"));
    TH2F *hPrimVtXY_evcut = static_cast<TH2F*>(mOutList->FindObject("hPrimVtXY_evcut"));
    TH2F *hVzVPDvsVzTPC_evcut = static_cast<TH2F*>(mOutList->FindObject("hVzVPDvsVzTPC_evcut"));
    TH1F *hVzVPDvsVzTPCDiff_evcut = static_cast<TH1F*>(mOutList->FindObject("hVzVPDvsVzTPCDiff_evcut"));

    TH1F *hpT_tr = static_cast<TH1F*>(mOutList->FindObject("hpT_tr"));
    TH1F *hpT_tr_cut = static_cast<TH1F*>(mOutList->FindObject("hpT_tr_cut"));
    TH1F *hEta_tr = static_cast<TH1F*>(mOutList->FindObject("hEta_tr"));
    TH1F *hEta_tr_cut = static_cast<TH1F*>(mOutList->FindObject("hEta_tr_cut"));
    //////////TH1F *hPhi_tr = static_cast<TH1F*>(mOutList->FindObject("hPhi_tr"));
    //////////TH1F *hPhi_tr_cut = static_cast<TH1F*>(mOutList->FindObject("hPhi_tr_cut"));
    TH2F *hEtaPhi_tr = static_cast<TH2F*>(mOutList->FindObject("hEtaPhi_tr"));
    TH2F *hEtaPhi_tr_cut = static_cast<TH2F*>(mOutList->FindObject("hEtaPhi_tr_cut"));
    TH1F *hNHitsFit = static_cast<TH1F*>(mOutList->FindObject("hNHitsFit"));
    TH1F *hNHitsFit_cut = static_cast<TH1F*>(mOutList->FindObject("hNHitsFit_cut"));
    //////////TH1F *hNHitsMax = static_cast<TH1F*>(mOutList->FindObject("hNHitsMax"));
    TH1F *hNHitsFitnHitsMax = static_cast<TH1F*>(mOutList->FindObject("hNHitsFitnHitsMax"));
    TH1F *hNHitsFitnHitsMax_cut = static_cast<TH1F*>(mOutList->FindObject("hNHitsFitnHitsMax_cut"));
    TH1F *hDca = static_cast<TH1F*>(mOutList->FindObject("hDca"));
    TH1F *hDca_cut = static_cast<TH1F*>(mOutList->FindObject("hDca_cut"));
    TH2F *hDcaXY = static_cast<TH2F*>(mOutList->FindObject("hDcaXY"));
    TH2F *hDcaXY_cut = static_cast<TH2F*>(mOutList->FindObject("hDcaXY_cut"));
    TH1F *hDcaZ = static_cast<TH1F*>(mOutList->FindObject("hDcaZ"));
    TH1F *hDcaZ_cut = static_cast<TH1F*>(mOutList->FindObject("hDcaZ_cut"));

    TH1F *hKaonPt = static_cast<TH1F*>(mOutList->FindObject("hKaonPt"));
    TH1F *hPionPt = static_cast<TH1F*>(mOutList->FindObject("hPionPt"));

    //////////TH1F *h_RefMult = static_cast<TH1F*>(mOutList->FindObject("h_RefMult"));
    //////////TH1F *h_gRefMult = static_cast<TH1F*>(mOutList->FindObject("h_gRefMult"));
    TH2F *h_gRefMultVsRefMult = static_cast<TH2F*>(mOutList->FindObject("h_gRefMultVsRefMult"));
    TH2F *h_gRefMultVsRefMult_evcut = static_cast<TH2F*>(mOutList->FindObject("h_gRefMultVsRefMult_evcut"));
    //////////TH1F *h_TofMult_matched = static_cast<TH1F*>(mOutList->FindObject("h_TofMult_matched"));
    TH2F *hTofMultVsRefMult = static_cast<TH2F*>(mOutList->FindObject("hTofMultVsRefMult"));
    TH2F *hTofMultVsRefMult_evcut = static_cast<TH2F*>(mOutList->FindObject("hTofMultVsRefMult_evcut"));

    TProfile *prof_nPrimaryTracks_vs_BBCx = static_cast<TProfile*>(mOutList->FindObject("prof_nPrimaryTracks_vs_BBCx"));
    TProfile *prof_nGoodTracks_vs_BBCx = static_cast<TProfile*>(mOutList->FindObject("prof_nGoodTracks_vs_BBCx"));
    TProfile *prof_ntracksKaonTPC_vs_BBCx = static_cast<TProfile*>(mOutList->FindObject("prof_ntracksKaonTPC_vs_BBCx"));
    TProfile *prof_ntracksPionTPC_vs_BBCx = static_cast<TProfile*>(mOutList->FindObject("prof_ntracksPionTPC_vs_BBCx"));
    TProfile *prof_ntracksKaon_vs_BBCx = static_cast<TProfile*>(mOutList->FindObject("prof_ntracksKaon_vs_BBCx"));
    TProfile *prof_ntracksPion_vs_BBCx = static_cast<TProfile*>(mOutList->FindObject("prof_ntracksPion_vs_BBCx"));

    TProfile *prof_nPrimaryTracks_vs_runID = static_cast<TProfile*>(mOutList->FindObject("prof_nPrimaryTracks_vs_runID"));
    TProfile *prof_nGoodTracks_vs_runID = static_cast<TProfile*>(mOutList->FindObject("prof_nGoodTracks_vs_runID"));
    TProfile *prof_ntracksKaonTPC_vs_runID = static_cast<TProfile*>(mOutList->FindObject("prof_ntracksKaonTPC_vs_runID"));
    TProfile *prof_ntracksPionTPC_vs_runID = static_cast<TProfile*>(mOutList->FindObject("prof_ntracksPionTPC_vs_runID"));
    TProfile *prof_ntracksKaon_vs_runID = static_cast<TProfile*>(mOutList->FindObject("prof_ntracksKaon_vs_runID"));
    TProfile *prof_ntracksPion_vs_runID = static_cast<TProfile*>(mOutList->FindObject("prof_ntracksPion_vs_runID"));

    /*TH2F *h_nPrimaryTracks_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_nPrimaryTracks_vs_BBCx"));
    TH2F *h_nGoodTracks_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_nGoodTracks_vs_BBCx"));
    TH2F *h_ntracksKaonTPC_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_ntracksKaonTPC_vs_BBCx"));
    TH2F *h_ntracksPionTPC_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_ntracksPionTPC_vs_BBCx"));
    TH2F *h_ntracksKaon_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_ntracksKaon_vs_BBCx"));
    TH2F *h_ntracksPion_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_ntracksPion_vs_BBCx"));

    TH2F *h_nPrimaryTrackswMatching_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_nPrimaryTrackswMatching_vs_BBCx"));
    TH2F *h_nGoodTrackswMatching_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_nGoodTrackswMatching_vs_BBCx"));
    TH2F *h_ntracksKaonTPCwMatching_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_ntracksKaonTPCwMatching_vs_BBCx"));
    TH2F *h_ntracksPionTPCwMatching_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_ntracksPionTPCwMatching_vs_BBCx"));
    
    //////////TProfile *prof_gRefmult_vs_BBCx = static_cast<TProfile*>(mOutList->FindObject("prof_gRefmult_vs_BBCx"));
    //////////TProfile *prof_gRefmult_vs_runID = static_cast<TProfile*>(mOutList->FindObject("prof_gRefmult_vs_runID"));

    TH2F *h_nPrimaryTracks_vs_runID = static_cast<TH2F*>(mOutList->FindObject("h_nPrimaryTracks_vs_runID"));
    TH2F *h_nGoodTracks_vs_runID = static_cast<TH2F*>(mOutList->FindObject("h_nGoodTracks_vs_runID"));
    TH2F *h_ntracksKaonTPC_vs_runID = static_cast<TH2F*>(mOutList->FindObject("h_ntracksKaonTPC_vs_runID"));
    TH2F *h_ntracksPionTPC_vs_runID = static_cast<TH2F*>(mOutList->FindObject("h_ntracksPionTPC_vs_runID"));
    TH2F *h_ntracksKaon_vs_runID = static_cast<TH2F*>(mOutList->FindObject("h_ntracksKaon_vs_runID"));
    TH2F *h_ntracksPion_vs_runID = static_cast<TH2F*>(mOutList->FindObject("h_ntracksPion_vs_runID"));

    TH2F *h_nPrimaryTrackswMatching_vs_runID = static_cast<TH2F*>(mOutList->FindObject("h_nPrimaryTrackswMatching_vs_runID"));
    TH2F *h_nGoodTrackswMatching_vs_runID = static_cast<TH2F*>(mOutList->FindObject("h_nGoodTrackswMatching_vs_runID"));
    TH2F *h_ntracksKaonTPCwMatching_vs_runID = static_cast<TH2F*>(mOutList->FindObject("h_ntracksKaonTPCwMatching_vs_runID"));
    TH2F *h_ntracksPionTPCwMatching_vs_runID = static_cast<TH2F*>(mOutList->FindObject("h_ntracksPionTPCwMatching_vs_runID"));

    
    //////////TH2F *h_ntracksKaonTPC_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_ntracksKaonTPC_vs_BBCx"));
    //////////TH2F *h_ntracksPionTPC_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_ntracksPionTPC_vs_BBCx"));

    

          
    //////////TH2F *h_nGoodTracks_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_nGoodTracks_vs_BBCx"));    
    TH2F *h_nPrimaryTracks_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_nPrimaryTracks_vs_BBCx"));
    TH2F *h_ntracksKaon_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_ntracksKaon_vs_BBCx"));
    TH2F *h_ntracksPion_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_ntracksPion_vs_BBCx"));
    TH2F *h_gRefmult_vs_BBCx = static_cast<TH2F*>(mOutList->FindObject("h_gRefmult_vs_BBCx"));
    //////////*/
    TH2F *hDedx_tr = static_cast<TH2F*>(mOutList->FindObject("hDedx_tr"));
    /////TH2D *h_nSigmadEdXPion = static_cast<TH2D*>(mOutList->FindObject("h_nSigmadEdXPion"));
    /////TH2D *h_nSigmadEdXKaon = static_cast<TH2D*>(mOutList->FindObject("h_nSigmadEdXKaon"));
    TH2F *hPiDedx_TOF = static_cast<TH2F*>(mOutList->FindObject("hPiDedx_TOF"));
    TH2F *hKDedx_TOF = static_cast<TH2F*>(mOutList->FindObject("hKDedx_TOF"));
    TH2F *hPiDedx_TOF_TPC = static_cast<TH2F*>(mOutList->FindObject("hPiDedx_TOF_TPC"));
    TH2F *hKDedx_TOF_TPC = static_cast<TH2F*>(mOutList->FindObject("hKDedx_TOF_TPC"));
    
    TH2F *hBetavsP_tr = static_cast<TH2F*>(mOutList->FindObject("hBetavsP_tr"));
    TH2F *hPiBetavsP_TPC = static_cast<TH2F*>(mOutList->FindObject("hPiBetavsP_TPC"));
    TH2F *hKBetavsP_TPC = static_cast<TH2F*>(mOutList->FindObject("hKBetavsP_TPC"));
    TH2F *hPiBetavsP_TOF_TPC = static_cast<TH2F*>(mOutList->FindObject("hPiBetavsP_TOF_TPC"));
    TH2F *hKBetavsP_TOF_TPC = static_cast<TH2F*>(mOutList->FindObject("hKBetavsP_TOF_TPC"));
    
    TH2F *h_nSigmaOneOverBetaPi_tr = static_cast<TH2F*>(mOutList->FindObject("h_nSigmaOneOverBetaPi_tr"));
    TH2F *h_nSigmaOneOverBetaK_tr = static_cast<TH2F*>(mOutList->FindObject("h_nSigmaOneOverBetaK_tr"));
    TH2F *h_nSigmaOneOverBetaPi_TPC = static_cast<TH2F*>(mOutList->FindObject("h_nSigmaOneOverBetaPi_TPC"));
    TH2F *h_nSigmaOneOverBetaK_TPC = static_cast<TH2F*>(mOutList->FindObject("h_nSigmaOneOverBetaK_TPC"));
    TH2F *h_nSigmaOneOverBetaPi_TOF_TPC = static_cast<TH2F*>(mOutList->FindObject("h_nSigmaOneOverBetaPi_TOF_TPC"));
    TH2F *h_nSigmaOneOverBetaK_TOF_TPC = static_cast<TH2F*>(mOutList->FindObject("h_nSigmaOneOverBetaK_TOF_TPC"));

    TH2F *h_DeltaOneOverBetaPi_tr = static_cast<TH2F*>(mOutList->FindObject("h_DeltaOneOverBetaPi_tr"));
    TH2F *h_DeltaOneOverBetaK_tr = static_cast<TH2F*>(mOutList->FindObject("h_DeltaOneOverBetaK_tr"));
    TH2F *h_DeltaOneOverBetaPi_TPC = static_cast<TH2F*>(mOutList->FindObject("h_DeltaOneOverBetaPi_TPC"));
    TH2F *h_DeltaOneOverBetaK_TPC = static_cast<TH2F*>(mOutList->FindObject("h_DeltaOneOverBetaK_TPC"));
    TH2F *h_DeltaOneOverBetaPi_TOF_TPC = static_cast<TH2F*>(mOutList->FindObject("h_DeltaOneOverBetaPi_TOF_TPC"));
    TH2F *h_DeltaOneOverBetaK_TOF_TPC = static_cast<TH2F*>(mOutList->FindObject("h_DeltaOneOverBetaK_TOF_TPC"));
    
    TH2F *h_nSigmadEdxPi_tr = static_cast<TH2F*>(mOutList->FindObject("h_nSigmadEdxPi_tr"));
    TH2F *h_nSigmadEdxK_tr = static_cast<TH2F*>(mOutList->FindObject("h_nSigmadEdxK_tr"));
    TH2F *h_nSigmadEdxPi_TOF = static_cast<TH2F*>(mOutList->FindObject("h_nSigmadEdxPi_TOF"));
    TH2F *h_nSigmadEdxK_TOF = static_cast<TH2F*>(mOutList->FindObject("h_nSigmadEdxK_TOF"));
    TH2F *h_nSigmadEdxPi_TOF_TPC = static_cast<TH2F*>(mOutList->FindObject("h_nSigmadEdxPi_TOF_TPC"));
    TH2F *h_nSigmadEdxK_TOF_TPC = static_cast<TH2F*>(mOutList->FindObject("h_nSigmadEdxK_TOF_TPC"));

   //====================================================================================================// 

    UInt_t nTracks = mPicoDst->numberOfTracks();
    ////UInt_t nTracks_Event = mPicoEvent->grefMult();

    ////cout << "nTracks" << "   " << nTracks << "   " << "nTracks" << "  " << nTracks_Event  << endl;
    nGoodTracks=0;
    
    Int_t nD0 = 0;
    Int_t nDstar = 0;
    Int_t nPrimary = 0;
    Int_t ntrackPionTPC = 0;
    Int_t ntrackKaonTPC = 0;
    Int_t ntrackPionTOF = 0;
    Int_t ntrackKaonTOF = 0;
    Int_t ntrackPionTOF_TPC = 0;
    Int_t ntrackKaonTOF_TPC = 0;
    
    Int_t nGoodTrackswMatching = 0;
    Int_t nPrimarywMatching = 0;
    Int_t ntrackPionTPCwMatching = 0;
    Int_t ntrackKaonTPCwMatching = 0;

    
 
    RunId = ConvertRunID(mPicoEvent->runId());
    float BBCx = mPicoEvent->BBCx() / 1000.;
    float refMult = mPicoEvent->refMult();
    float gRefMult = mPicoEvent->grefMult();
    float tofMult = mPicoEvent->nBTOFMatch();

   hPrimVtZ->Fill(mPrimVtx.z());
   hPrimVtXY->Fill(mPrimVtx.x(), mPrimVtx.y());
   hVzVPDvsVzTPC->Fill(mPrimVtx.z(), mPicoEvent->vzVpd());
   hVzVPDvsVzTPCDiff->Fill(mPrimVtx.z() - mPicoEvent->vzVpd());
   hTofMultVsRefMult->Fill(tofMult, refMult);
   h_gRefMultVsRefMult->Fill(gRefMult, refMult);

    if (fabs(mPrimVtx.z()) < mHFCuts->getCutVzMax()) {
    
    /*//////////hPrimVtZ_beforematch->Fill(mPrimVtx.z());
    hPrimVtXY_beforematch->Fill(mPrimVtx.x(), mPrimVtx.y());
    hVzVPDvsVzTPC_beforematch->Fill(mPrimVtx.z(), mPicoEvent->vzVpd());
    hVzVPDvsVzTPCDiff_beforematch->Fill(mPrimVtx.z() - mPicoEvent->vzVpd());
    

    if (fabs(mPrimVtx.z() - mPicoEvent->vzVpd()) < mHFCuts->getCutVzVpdVzMax()) {
         hPrimVtZ_checkvzvpdtpc->Fill(mPrimVtx.z());
         hPrimVtXY_checkvzvpdtpc->Fill(mPrimVtx.x(), mPrimVtx.y());
         hVzVPDvsVzTPC_checkvzvpdtpc->Fill(mPrimVtx.z(), mPicoEvent->vzVpd());
         hVzVPDvsVzTPCDiff_checkvzvpdtpc->Fill(mPrimVtx.z() - mPicoEvent->vzVpd());
    }
    *//////////

    if (mHFCuts->isMatchedFast(mPicoDst)) {

        //-------------- Filling Histogram before event cuts --------------------//
        
        /*//////////
        hPrimVtZ_wmatch->Fill(mPrimVtx.z());
        hPrimVtXY_wmatch->Fill(mPrimVtx.x(), mPrimVtx.y());
        hVzVPDvsVzTPC_wmatch->Fill(mPrimVtx.z(), mPicoEvent->vzVpd());
        hVzVPDvsVzTPCDiff_wmatch->Fill(mPrimVtx.z() - mPicoEvent->vzVpd());
        *//////////

        //////
        
        if (mHFCuts->isBetterEvent(mPicoDst)){

            //-------------- Filling Histogram after event cuts --------------------//


            hPrimVtZ_evcut->Fill(mPrimVtx.z());
            hPrimVtXY_evcut->Fill(mPrimVtx.x(), mPrimVtx.y());
            hVzVPDvsVzTPC_evcut->Fill(mPrimVtx.z(), mPicoEvent->vzVpd());
            hVzVPDvsVzTPCDiff_evcut->Fill(mPrimVtx.z() - mPicoEvent->vzVpd());

            //////////h_RefMult->Fill(refMult);
            //////////h_gRefMult->Fill(gRefMult);
            
            h_gRefMultVsRefMult_evcut->Fill(gRefMult, refMult);
            hTofMultVsRefMult_evcut->Fill(tofMult, refMult);

            //////////h_gRefmult_vs_BBCx->Fill(BBCx, gRefMult);
            //////////prof_gRefmult_vs_BBCx->Fill(BBCx, gRefMult);


            for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
            StPicoTrack* trk = mPicoDst->track(iTrack);
            //////int bTOFPidTraitsIndex = trk->bTofPidTraitsIndex();
            //////StPicoBTofPidTraits* tofPidTraits = mPicoDst->btofPidTraits(bTOFPidTraitsIndex);

            float nsigmaTOFpion = mHFCuts->getnSigmaTOF(trk, StPicoCutsBase::kPion);
            float nsigmaTOFkaon = mHFCuts->getnSigmaTOF(trk, StPicoCutsBase::kKaon);

                TVector3 pMom = trk->pMom(); //primary momentum
                float invBeta = mHFCuts->getTofBetaBase(trk);
                float BetaPion = mHFCuts->getTofBetaBase(trk);
                float DeltaOneoverBeta_pi = mHFCuts->getOneOverBeta(trk,BetaPion,StPicoCutsBase::kPion) ; // /0.011;
                float BetaKaon = mHFCuts->getTofBetaBase(trk);
                float DeltaOneoverBeta_K = mHFCuts->getOneOverBeta(trk,BetaKaon,StPicoCutsBase::kKaon); // /0.012;
                

       

                if (trk->isPrimary()) { 
              ///if (mHFCuts->isTOFmatched(trk) || mHFCuts->isBEMCmatched(trk))
              //{
              nPrimary++;
              primaryTracks.push_back(iTrack);
           
                    //----------------- Filling histogram before track cut --------------------//

                    hpT_tr->Fill(trk->pPt());
                    hEta_tr->Fill(pMom.PseudoRapidity());
                    //////////hPhi_tr->Fill(pMom.Phi());
                    hEtaPhi_tr->Fill(pMom.Phi(), pMom.PseudoRapidity());
                    hNHitsFit->Fill(trk->nHitsFit());
                    //////////hNHitsMax->Fill(trk->nHitsMax());
                    hNHitsFitnHitsMax->Fill((double) trk->nHitsFit() / (double) trk->nHitsMax());
                    hDca->Fill((mPrimVtx - trk->origin()).Mag());
                    hDcaXY->Fill((mPrimVtx - trk->origin()).x(), (mPrimVtx - trk->origin()).y());
                    hDcaZ->Fill((mPrimVtx - trk->origin()).z());

                    if (mHFCuts->isTOFmatched(trk) || mHFCuts->isBEMCmatched(trk)) nPrimarywMatching++;

            if (mHFCuts->isGoodTrack(trk)){
              
                        nGoodTracks++;
            
               
                        //----------------------- Filling Histogram after track cut --------------//
              
                        hpT_tr_cut->Fill(trk->pPt());
                        hEta_tr_cut->Fill(pMom.PseudoRapidity());
                        //////////hPhi_tr_cut->Fill(pMom.Phi());
                        hEtaPhi_tr_cut->Fill(pMom.Phi(), pMom.PseudoRapidity());
                        hNHitsFit_cut->Fill(trk->nHitsFit());
                        hNHitsFitnHitsMax_cut->Fill((double) trk->nHitsFit() / (double) trk->nHitsMax());
                        hDca_cut->Fill((mPrimVtx - trk->origin()).Mag());
                        hDcaXY_cut->Fill((mPrimVtx - trk->origin()).x(), (mPrimVtx - trk->origin()).y());
                        hDcaZ_cut->Fill((mPrimVtx - trk->origin()).z());
              
                        hDedx_tr->Fill(pMom.Mag() * trk->charge(), trk->dEdx());
                        hBetavsP_tr->Fill(pMom.Mag(), 1.0 / invBeta);

                        h_nSigmadEdxPi_tr->Fill(pMom.Mag(), trk->nSigmaPion());
                        h_nSigmadEdxK_tr->Fill(pMom.Mag(), trk->nSigmaKaon());

                        h_nSigmaOneOverBetaPi_tr->Fill(pMom.Mag(), nsigmaTOFpion);
                        h_nSigmaOneOverBetaK_tr->Fill(pMom.Mag(), nsigmaTOFkaon);

                        h_DeltaOneOverBetaPi_tr->Fill(pMom.Mag(), DeltaOneoverBeta_pi);
                        h_DeltaOneOverBetaK_tr->Fill(pMom.Mag(), DeltaOneoverBeta_K);
                        
                        if (mHFCuts->isTOFmatched(trk) || mHFCuts->isBEMCmatched(trk)) nGoodTrackswMatching++;
                        

                        //-------------------- Filling TOF histogram after TPC cuts---------------//
               
              if (mHFCuts->isPionTPC(trk)){
                /////if (mHFCuts->isTPCPion(trk)){
                /////if (mHFCuts->isTOFmatched(trk) || mHFCuts->isBEMCmatched(trk))
                /////if (mHFCuts->isTOFmatched(trk)) 
                ///{
                            ntrackPionTPC++;
                
                ///float BBCPion = mPicoEvent->BBCx() / 1000.;
                ///float grefMultPion = mPicoEvent->grefMult();
                ////h_QA_OneOverBetaDiffPion->Fill(trk->pPt(), npi1_TOFinvbeta);
                            hPiBetavsP_TPC->Fill(pMom.Mag(), 1.0 / BetaPion);
                            h_nSigmaOneOverBetaPi_TPC->Fill(pMom.Mag(), nsigmaTOFpion);
                            h_DeltaOneOverBetaPi_TPC->Fill(pMom.Mag(), DeltaOneoverBeta_pi);

                            if (mHFCuts->isTOFmatched(trk) || mHFCuts->isBEMCmatched(trk)) ntrackPionTPCwMatching++;
                        
                ////}
              }
            
            
             if (mHFCuts->isKaonTPC(trk)){
             /////if (mHFCuts->isTPCKaon(trk)){
                ////if (mHFCuts->isTOFmatched(trk) || mHFCuts->isBEMCmatched(trk))
                /////if (mHFCuts->isTOFmatched(trk)) 
                ////{
                            ntrackKaonTPC++;

                ////float BBCKaon = mPicoEvent->BBCx() / 1000.;
                ////float grefMultKaon = mPicoEvent->grefMult();
                /////h_ntracks_vs_BBCx_Kaon->Fill(BBCKaon, grefMultKaon);
                /////prof_ntracks_vs_BBCx_Kaon->Fill(BBCKaon, grefMultKaon);
                ////h_QA_OneOverBetaDiffKaon->Fill(trk->pPt(), nk_TOFinvbeta);
                            hKBetavsP_TPC->Fill(pMom.Mag(), 1.0 / BetaKaon);
                            h_nSigmaOneOverBetaK_TPC->Fill(pMom.Mag(), nsigmaTOFkaon);
                            h_DeltaOneOverBetaK_TPC->Fill(pMom.Mag(), DeltaOneoverBeta_K);

                            if (mHFCuts->isTOFmatched(trk) || mHFCuts->isBEMCmatched(trk)) ntrackKaonTPCwMatching++;
                /////}
             }

                        //-------------------------- Filling TPC histograms after TOF cuts ----------------//

             if (mHFCuts->isPionTOF(trk)){

                            ntrackPionTOF++;

                            hPiDedx_TOF->Fill(pMom.Mag() * trk->charge(), trk->dEdx());
                            h_nSigmadEdxPi_TOF->Fill(pMom.Mag(), trk->nSigmaPion());
              }

              if (mHFCuts->isKaonTOF(trk)){

                            ntrackKaonTOF++;

                            hKDedx_TOF->Fill(pMom.Mag() * trk->charge(), trk->dEdx());
                            h_nSigmadEdxK_TOF->Fill(pMom.Mag(), trk->nSigmaKaon());

              }

                        //----------------------- Filling all histograms after TOF&TPC cuts --------------//

             if (mHFCuts->isGoodPion(trk)) {
                /////if (mHFCuts->isTOFmatched(trk) || mHFCuts->isBEMCmatched(trk))
                /////if (mHFCuts->isTOFmatched(trk)) 
                ///{
                    if (trk->charge() != 0) {
                    ntrackPionTOF_TPC++;
                mIdxPicoPions.push_back(iTrack);
                
                            hPiDedx_TOF_TPC->Fill(pMom.Mag() * trk->charge(), trk->dEdx());
                            hPiBetavsP_TOF_TPC->Fill(pMom.Mag(), 1.0 / BetaPion);
                            hPionPt->Fill(trk->pPt());
                            h_nSigmaOneOverBetaPi_TOF_TPC->Fill(pMom.Mag(), nsigmaTOFpion);
                            h_nSigmadEdxPi_TOF_TPC->Fill(pMom.Mag(),trk->nSigmaPion());
                            h_DeltaOneOverBetaPi_TOF_TPC->Fill(pMom.Mag(), DeltaOneoverBeta_pi);
                        
                }
             }

             if (mHFCuts->isGoodKaon(trk)){
                /////if (mHFCuts->isTOFmatched(trk) || mHFCuts->isBEMCmatched(trk))
                /////if (mHFCuts->isTOFmatched(trk)) 
                ///{
                    if (trk->charge() != 0) {
                    ntrackKaonTOF_TPC++;
                mIdxPicoKaons.push_back(iTrack);

                            hKDedx_TOF_TPC->Fill(pMom.Mag() * trk->charge(), trk->dEdx());
                            hKBetavsP_TOF_TPC->Fill(pMom.Mag(), 1.0 / BetaKaon);
                            hKaonPt->Fill(trk->pPt());
                            h_nSigmaOneOverBetaK_TOF_TPC->Fill(pMom.Mag(), nsigmaTOFkaon);
                            h_nSigmadEdxK_TOF_TPC->Fill(pMom.Mag(),trk->nSigmaKaon());
                            h_DeltaOneOverBetaK_TOF_TPC->Fill(pMom.Mag(), DeltaOneoverBeta_K);
                }
              }
            }
            else tracksToRemove.push_back(iTrack);
        }

        
    }
            

            //------------------------- Histograms for Pile-up check -------------------// 
            
            //////////h_nPrimaryTracks_vs_BBCx->Fill(BBCx, nPrimary);
            prof_nPrimaryTracks_vs_BBCx->Fill(BBCx, nPrimary);
            //////////h_nPrimaryTracks_vs_runID->Fill(RunId, nPrimary);
            prof_nPrimaryTracks_vs_runID->Fill(RunId, nPrimary);
    
            //////////h_nGoodTracks_vs_BBCx->Fill(BBCx, nGoodTracks);
            prof_nGoodTracks_vs_BBCx->Fill(BBCx, nGoodTracks);
            //////////h_nGoodTracks_vs_runID->Fill(RunId, nGoodTracks);
            prof_nGoodTracks_vs_runID->Fill(RunId, nGoodTracks);

            //////////h_ntracksKaonTPC_vs_BBCx->Fill(BBCx, ntrackKaonTPC);
            prof_ntracksKaonTPC_vs_BBCx->Fill(BBCx, ntrackKaonTPC);
            //////////h_ntracksKaonTPC_vs_runID->Fill(RunId, ntrackKaonTPC);
            prof_ntracksKaonTPC_vs_runID->Fill(RunId, ntrackKaonTPC);

            //////////h_ntracksPionTPC_vs_BBCx->Fill(BBCx, ntrackPionTPC);
            prof_ntracksPionTPC_vs_BBCx->Fill(BBCx, ntrackPionTPC);
            //////////h_ntracksPionTPC_vs_runID->Fill(RunId, ntrackPionTPC);
            prof_ntracksPionTPC_vs_runID->Fill(RunId, ntrackPionTPC);

            //////////h_ntracksKaon_vs_BBCx->Fill(BBCx, ntrackKaonTOF_TPC);
            prof_ntracksKaon_vs_BBCx->Fill(BBCx, ntrackKaonTOF_TPC);
            //////////h_ntracksKaon_vs_runID->Fill(RunId, ntrackKaonTOF_TPC);
            prof_ntracksKaon_vs_runID->Fill(RunId, ntrackKaonTOF_TPC);

            

            //////////h_ntracksPion_vs_BBCx->Fill(BBCx, ntrackPionTOF_TPC);
            prof_ntracksPion_vs_BBCx->Fill(BBCx, ntrackPionTOF_TPC);
            //////////h_ntracksPion_vs_runID->Fill(RunId, ntrackPionTOF_TPC);
            prof_ntracksPion_vs_runID->Fill(RunId, ntrackPionTOF_TPC);

            /*///////// h_nPrimaryTrackswMatching_vs_BBCx->Fill(BBCx, nPrimarywMatching);
            //////////prof_nPrimaryTracks_vs_BBCx->Fill(BBCx, nPrimary);
            //////////h_nPrimaryTrackswMatching_vs_runID->Fill(RunId, nPrimarywMatching);
    
            h_nGoodTrackswMatching_vs_BBCx->Fill(BBCx, nGoodTrackswMatching);
            //////////prof_nGoodTracks_vs_BBCx->Fill(BBCx, nGoodTracks);
            h_nGoodTrackswMatching_vs_runID->Fill(RunId, nGoodTrackswMatching);

            h_ntracksKaonTPCwMatching_vs_BBCx->Fill(BBCx, ntrackKaonTPCwMatching);
            //////////prof_ntracksKaonTPC_vs_BBCx->Fill(BBCx, ntrackKaonTPC);
            h_ntracksKaonTPCwMatching_vs_runID->Fill(RunId, ntrackKaonTPCwMatching);

            h_ntracksPionTPCwMatching_vs_BBCx->Fill(BBCx, ntrackPionTPCwMatching);
            //////////prof_ntracksPionTPC_vs_BBCx->Fill(BBCx, ntrackPionTPC);
            h_ntracksPionTPCwMatching_vs_runID->Fill(RunId, ntrackPionTPCwMatching);*/

            

            


    TVector3 useVertex(mPrimVtx.x(), mPrimVtx.y(), mPrimVtx.z());
    if (mSwitchRefit) useVertex=refitVertex(true);

//        if (tracksToRemove.size()==nPrimary) {
//            cout<<useVertex.x()<<" "<<useVertex.y()<<" "<<useVertex.z()<<endl;
//            cout<<mPrimVtx.x()<<" "<<mPrimVtx.y()<<" "<<mPrimVtx.z()<<endl;
//        }

    //////////vect_daughter_id = {{0.5, 0.5}};
    //////////daughter_id;

    if (mIdxPicoPions.size() > 0 && mIdxPicoKaons.size() > 0) {

    

    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
        StPicoTrack *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
        for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
            if ( mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion1] ) continue;
            //////////daughter_id = {mIdxPicoPions[idxPion1], mIdxPicoKaons[idxKaon]};
            //////////if (isVectorInVectorOfVectors(vect_daughter_id, daughter_id)) continue;
            StPicoTrack *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);
            StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion),mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), mIdxPicoPions[idxPion1],mIdxPicoKaons[idxKaon], useVertex, mBField, kTRUE);
            StHFRotPair *rotpair = new StHFRotPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion),mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), mIdxPicoPions[idxPion1],mIdxPicoKaons[idxKaon], useVertex, mBField, kTRUE);

            /////if (!mHFCuts->isGoodSecondaryVertexPair(pair)) continue;
            /////if((pair->pt())<=0) continue;
			/////if((TMath::Abs(pair->rapidity()))>1) continue;
            ////if(pair->cosThetaStar() == 1) continue;

            /////bool isD0 = false;
            ////bool isBgD0 = false;
            
            int charge = kaon->charge() + pion1->charge();
            /////if(charge == 0) isD0=true;

            /*if(kaon->charge() == 1 && pion1->charge() == -1 ) isD0=true;
            if(kaon->charge() == -1 && pion1->charge() == 1 ) isD0=true;

            if(kaon->charge() == 1 && pion1->charge() == 1 ) isBgD0=true;
            if(kaon->charge() == -1 && pion1->charge() == -1 ) isBgD0=true;
            */

            Float_t ispion1TOFmatched = 0;
            Float_t ispion1BEMCmatched = 0;
            Float_t iskaonTOFmatched = 0;
            Float_t iskaonBEMCmatched = 0;

            if (mHFCuts->isTOFmatched(pion1)) ispion1TOFmatched = 1;
            if (mHFCuts->isBEMCmatched(pion1)) ispion1BEMCmatched = 1;
            if (mHFCuts->isTOFmatched(kaon)) iskaonTOFmatched = 1;
            if (mHFCuts->isBEMCmatched(kaon)) iskaonBEMCmatched = 1;

            //////////Float_t hotSpot=0;
            //////////if (mHFCuts->checkHotSpot(&mPrimVtx)) hotSpot=1;

            float pairM = pair->m();

            const int nNtVars = ntp_DMeson_UnlikeSign->GetNvar();
            float ntVar[nNtVars];
            int ii=0;

            /*float VzTPC = mPrimVtx.z();
            float VzVPD = mPicoEvent->vzVpd();
            float VzVPD_VzTPC = VzVPD - VzTPC;*/

            /*//////////ntVar[ii++] = mPicoEvent->grefMult();
            
            ntVar[ii++] = mPicoEvent->runId();
            ntVar[ii++] = mPicoEvent->eventId();
            ntVar[ii++] = mPicoEvent->ZDCx();
            ntVar[ii++] = mPicoEvent->BBCx();
            ntVar[ii++] = hotSpot;

            ntVar[ii++] = primary;
            ntVar[ii++] = nGoodTracks;
            *///////////
            ntVar[ii++] = mPicoEvent->runId();

            ntVar[ii++] = pion1->pPt();
            ntVar[ii++] = pion1->pPtot();
            //////////ntVar[ii++] = pair->particle1Dca();
            ntVar[ii++] = pion1->nSigmaPion();
            ntVar[ii++] = mHFCuts->getnSigmaTOF(pion1, StPicoCutsBase::kPion);
            //////////ntVar[ii++] = pion1->nHits();
            //////////ntVar[ii++] = pion1->nHitsFit();
            ntVar[ii++] = mHFCuts->getOneOverBeta(pion1,mHFCuts->getTofBetaBase(pion1),StPicoCutsBase::kPion);
            //////////ntVar[ii++] = mHFCuts->getTofBetaBase(pion1);
            ntVar[ii++] = pion1->charge();
            ntVar[ii++] = ispion1TOFmatched;
            ntVar[ii++] = ispion1BEMCmatched;

            ntVar[ii++] = kaon->pPt();
            ntVar[ii++] = kaon->pPtot();
            //////////
            ntVar[ii++] = kaon->nSigmaKaon();
            ntVar[ii++] = mHFCuts->getnSigmaTOF(kaon, StPicoCutsBase::kKaon);
            //////////ntVar[ii++] = kaon->nHits();
            //////////ntVar[ii++] = kaon->nHitsFit();
            ntVar[ii++] = mHFCuts->getOneOverBeta(kaon,mHFCuts->getTofBetaBase(kaon),StPicoCutsBase::kKaon);
            //////////ntVar[ii++] = mHFCuts->getTofBetaBase(kaon);
            ntVar[ii++] = kaon->charge();
            ntVar[ii++] = iskaonTOFmatched;
            ntVar[ii++] = iskaonBEMCmatched;

            ntVar[ii++] = pair->dcaDaughters();
            

            ntVar[ii++] = pair->pointingAngle();
            //////////ntVar[ii++] = cos(pair->pointingAngle());
            ntVar[ii++] = pair->decayLength();
            ntVar[ii++] = pair->DcaToPrimaryVertex(); //(pair->decayLength())*sin(pair->pointingAngle());
            ntVar[ii++] = pair->cosThetaStar();

            ntVar[ii++] = pair->pt();
            ntVar[ii++] = pair->m();
            ntVar[ii++] = pair->rapidity();
            ntVar[ii++] = pair->phi();
            ntVar[ii++] = pair->eta();

            ntVar[ii++] = mPicoEvent->refMult();
            ntVar[ii++] = kaon->nHitsFit();
            ntVar[ii++] = pion1->nHitsFit();
            ntVar[ii++] = kaon->nHitsDedx();
            ntVar[ii++] = pion1->nHitsDedx();
            ntVar[ii++] = mHFCuts->getbtofYLocal(kaon);
            ntVar[ii++] = mHFCuts->getbtofYLocal(pion1);

            ntVar[ii++] = mPrimVtx.z();
            ntVar[ii++] = mPicoEvent->vzVpd();

            ntVar[ii++] = pair->particle2Dca();
            ntVar[ii++] = pair->particle1Dca();






            //////////ntVar[ii++] = VzVPD_VzTPC;

            const int nNtVars2 = ntp_DMeson_Rotated->GetNvar();
            float ntVar2[nNtVars2];
            int iiii=0;

            /*//////////ntVar2[iiii++] = mPicoEvent->grefMult();
            ntVar2[iiii++] = mPicoEvent->refMult();
            ntVar2[iiii++] = mPicoEvent->runId();
            ntVar2[iiii++] = mPicoEvent->eventId();
            ntVar2[iiii++] = mPicoEvent->ZDCx();
            ntVar2[iiii++] = mPicoEvent->BBCx();
            ntVar2[iiii++] = hotSpot;

            ntVar2[iiii++] = primary;
            ntVar2[iiii++] = nGoodTracks;
            *///////////

            ntVar2[iiii++] = mPicoEvent->runId();

            ntVar2[iiii++] = pion1->pPt();
            ntVar2[iiii++] = pion1->pPtot();
            //////////ntVar2[iiii++] = rotpair->particle1Dca();
            ntVar2[iiii++] = pion1->nSigmaPion();
            ntVar2[iiii++] = mHFCuts->getnSigmaTOF(pion1, StPicoCutsBase::kPion);
            //////////ntVar2[iiii++] = pion1->nHitsFit();
            ntVar2[iiii++] = mHFCuts->getOneOverBeta(pion1,mHFCuts->getTofBetaBase(pion1),StPicoCutsBase::kPion);
            //////////ntVar2[iiii++] = mHFCuts->getTofBetaBase(pion1);
            ntVar2[iiii++] = pion1->charge();
            ntVar2[iiii++] = ispion1TOFmatched;
            ntVar2[iiii++] = ispion1BEMCmatched;

            ntVar2[iiii++] = kaon->pPt();
            ntVar2[iiii++] = kaon->pPtot();
            //////////ntVar2[iiii++] = rotpair->particle2Dca();
            ntVar2[iiii++] = kaon->nSigmaKaon();
            ntVar2[iiii++] = mHFCuts->getnSigmaTOF(kaon, StPicoCutsBase::kKaon);
            //////////ntVar2[iiii++] = kaon->nHits();
            //////////ntVar2[iiii++] = kaon->nHitsFit();

            ntVar2[iiii++] = mHFCuts->getOneOverBeta(kaon,mHFCuts->getTofBetaBase(kaon),StPicoCutsBase::kKaon);
            //////////ntVar2[iiii++] = mHFCuts->getTofBetaBase(kaon);
            ntVar2[iiii++] = kaon->charge();
            ntVar2[iiii++] = iskaonTOFmatched;
            ntVar2[iiii++] = iskaonBEMCmatched;

            ntVar2[iiii++] = rotpair->dcaDaughters();
            //////////ntVar2[iiii++] = mPrimVtx.z();
            //////////ntVar2[iiii++] = mPicoEvent->vzVpd();

            ntVar2[iiii++] = rotpair->pointingAngle();
            //////////ntVar2[iiii++] = cos(rotpair->pointingAngle());
            ntVar2[iiii++] = rotpair->decayLength();
            ntVar2[iiii++] = rotpair->DcaToPrimaryVertex(); //(rotpair->decayLength())*sin(rotpair->pointingAngle());
            ntVar2[iiii++] = rotpair->cosThetaStar();

            ntVar2[iiii++] = rotpair->pt();
            ntVar2[iiii++] = rotpair->m();
            ntVar2[iiii++] = rotpair->rapidity();
            ntVar2[iiii++] = rotpair->phi();
            ntVar2[iiii++] = rotpair->eta();

            //////////ntVar2[iiii++] = VzVPD_VzTPC;

            ntVar2[iiii++] = mPicoEvent->refMult();
            ntVar2[iiii++] = kaon->nHitsFit();
            ntVar2[iiii++] = pion1->nHitsFit();
            ntVar2[iiii++] = kaon->nHitsDedx();
            ntVar2[iiii++] = pion1->nHitsDedx();
            ntVar2[iiii++] = mHFCuts->getbtofYLocal(kaon);
            ntVar2[iiii++] = mHFCuts->getbtofYLocal(pion1);

            ntVar2[iiii++] = mPrimVtx.z();
            ntVar2[iiii++] = mPicoEvent->vzVpd();

            ntVar2[iiii++] = rotpair->particle2Dca();
            ntVar2[iiii++] = rotpair->particle1Dca();





            
            if (charge == 0) {
                ntp_DMeson_UnlikeSign->Fill(ntVar);
                ntp_DMeson_Rotated->Fill(ntVar2);
                //////////vect_daughter_id.push_back(daughter_id);
                nD0++;
            
            if((pairM > 1.84) && (pairM < 1.89)) {

                for (unsigned short idxPion2 = 0; idxPion2 < mIdxPicoPions.size(); ++idxPion2)
                {

                    StPicoTrack *pion2 = mPicoDst->track(mIdxPicoPions[idxPion2]);
                    if (mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion1] || mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion2] || mIdxPicoPions[idxPion1] == mIdxPicoPions[idxPion2]) continue;

                    StHFTriplet *triplet = new StHFTriplet(pion1, kaon, pion2, mHFCuts->getHypotheticalMass(StHFCuts::kPion), mHFCuts->getHypotheticalMass(StHFCuts::kKaon), mHFCuts->getHypotheticalMass(StHFCuts::kPion), mIdxPicoPions[idxPion1], mIdxPicoKaons[idxKaon], mIdxPicoPions[idxPion2], mPrimVtx, mBField);


                    ///bool isDstar = false;
                    // isDstar = true;

                    Float_t ispion2TOFmatched = 0;
                    Float_t ispion2BEMCmatched = 0;

                    if (mHFCuts->isTOFmatched(pion2)) ispion2TOFmatched = 1;
                    if (mHFCuts->isBEMCmatched(pion2)) ispion2BEMCmatched = 1;

                    /*Float_t hotSpot2 = 0;
                    if (mHFCuts->checkHotSpot(&mPrimVtx)) hotSpot2 = 1;
                    
                    float tripletM = triplet->m();
                    float triplet_pair = tripletM - pairM;
                    */

                    const int nNtVars3 = ntp_DstarMeson_RightSign->GetNvar();
                    float ntVar3[nNtVars3];
                    int iii = 0;

                    /*ntVar3[iii++] = mPicoEvent->grefMult();
                    ntVar3[iii++] = mPicoEvent->refMult();
                    ntVar3[iii++] = mPicoEvent->runId();
                    ntVar3[iii++] = mPicoEvent->eventId();
                    ntVar3[iii++] = mPicoEvent->ZDCx();
                    ntVar3[iii++] = mPicoEvent->BBCx();
                    ntVar3[iii++] = hotSpot2;

                    ntVar3[iii++] = primary2;
                    ntVar3[iii++] = nGoodTracks;
                    */
                    ntVar3[iii++] = mPicoEvent->runId();
                    ntVar3[iii++] = pion1->pPt();
                    ntVar3[iii++] = pion1->pPtot();
                    //////////ntVar3[iii++] = triplet->particle1Dca();
                    ntVar3[iii++] = pion1->nSigmaPion();
                    //////////ntVar3[iii++] = (float)((float)pion1->nHitsFit()/(float)pion1->nHitsMax());
                    //////////ntVar3[iii++] = pion1->nHitsFit();
                    ntVar3[iii++] = mHFCuts->getnSigmaTOF(pion1, StPicoCutsBase::kPion);
                    ntVar3[iii++] = mHFCuts->getOneOverBeta(pion1,mHFCuts->getTofBetaBase(pion1),StPicoCutsBase::kPion);
                    ntVar3[iii++] = pion1->charge();
                    ntVar3[iii++] = ispion1TOFmatched;
                    ntVar3[iii++] = ispion1BEMCmatched;
                    //////////ntVar3[iii++] = pion1->pMom().PseudoRapidity();

                    ntVar3[iii++] = kaon->pPt();
                    ntVar3[iii++] = kaon->pPtot();
                    //////////ntVar3[iii++] = triplet->particle2Dca();
                    ntVar3[iii++] = kaon->nSigmaKaon();
                    //////////ntVar3[iii++] = (float)((float)kaon->nHitsFit()/(float)kaon->nHitsMax());
                    //////////ntVar3[iii++] = kaon->nHitsFit();
                    ntVar3[iii++] = mHFCuts->getnSigmaTOF(kaon, StPicoCutsBase::kKaon);
                    ntVar3[iii++] = mHFCuts->getOneOverBeta(kaon,mHFCuts->getTofBetaBase(kaon),StPicoCutsBase::kKaon);
                    ntVar3[iii++] = kaon->charge();
                    ntVar3[iii++] = iskaonTOFmatched;
                    ntVar3[iii++] = iskaonBEMCmatched;
                    //////////ntVar3[iii++] = kaon->pMom().PseudoRapidity();


                    ntVar3[iii++] = pion2->pPt();
                    ntVar3[iii++] = pion2->pPtot();
                    //////////ntVar3[iii++] = triplet->particle3Dca();
                    ntVar3[iii++] = pion2->nSigmaPion();
                    //////////ntVar3[iii++] = (float)((float)pion2->nHitsFit()/(float)pion2->nHitsMax());
                    //////////ntVar3[iii++] = pion2->nHitsFit();
                    ntVar3[iii++] = mHFCuts->getnSigmaTOF(pion2, StPicoCutsBase::kPion);
                    ntVar3[iii++] = mHFCuts->getOneOverBeta(pion2,mHFCuts->getTofBetaBase(pion2),StPicoCutsBase::kPion);
                    ntVar3[iii++] = pion2->charge();
                    ntVar3[iii++] = ispion2TOFmatched;
                    ntVar3[iii++] = ispion2BEMCmatched;
                    //////////ntVar3[iii++] = pion2->pMom().PseudoRapidity();


                    //////////ntVar3[iii++] = mPrimVtx.z();
                    //////////ntVar3[iii++] = mPicoEvent->vzVpd();
                    //////////ntVar3[iii++] = sqrt(pow(mPrimVtx.x(),2) + pow(mPrimVtx.y(),2));

                    ntVar3[iii++] = pair->m();
                    ntVar3[iii++] = pair->pt();
                    ntVar3[iii++] = pair->rapidity();
                    ntVar3[iii++] = pair->cosThetaStar();
                    ntVar3[iii++] = pair->pointingAngle();

                    
                    ntVar3[iii++] = triplet->pt();
                    ntVar3[iii++] = triplet->m();
                    ntVar3[iii++] = triplet->rapidity();
                    //////////ntVar3[iii++] = triplet_pair;

                    ntVar3[iii++] = mPicoEvent->refMult();
                    ntVar3[iii++] = kaon->nHitsFit();
                    ntVar3[iii++] = pion1->nHitsFit();
                    ntVar3[iii++] = pion2->nHitsFit();
                    ntVar3[iii++] = kaon->nHitsDedx();
                    ntVar3[iii++] = pion1->nHitsDedx();
                    ntVar3[iii++] = pion2->nHitsDedx();
                    ntVar3[iii++] = mHFCuts->getbtofYLocal(kaon);
                    ntVar3[iii++] = mHFCuts->getbtofYLocal(pion1);
                    ntVar3[iii++] = mHFCuts->getbtofYLocal(pion2);

                    ntVar3[iii++] = mPrimVtx.z();
                    ntVar3[iii++] = mPicoEvent->vzVpd();

                    ntVar3[iii++] = triplet->particle2Dca();
                    ntVar3[iii++] = triplet->particle1Dca();
                    ntVar3[iii++] = triplet->particle3Dca();



                    
                    
                    if ((pion1->charge() == pion2->charge())){
                        ntp_DstarMeson_RightSign->Fill(ntVar3);
                        nDstar++;
                    } 
                    
                    if ((pion1->charge() == -1*(pion2->charge()))) {
                      ntp_DstarMeson_WrongSign->Fill(ntVar3);
                    }
                    


                }  //  for (unsigned short idxPion2 = idxPion1+1; idxPion2 < mIdxPicoPions.size(); ++idxPion2)

            }

             if((pairM>1.7 && pairM<1.8)||(pairM>1.92 && pairM<2.02)) {

                for (unsigned short idxPion2 = 0; idxPion2 < mIdxPicoPions.size(); ++idxPion2)
                {

                    StPicoTrack *pion2 = mPicoDst->track(mIdxPicoPions[idxPion2]);
                    if (mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion1] || mIdxPicoKaons[idxKaon] == mIdxPicoPions[idxPion2] || mIdxPicoPions[idxPion1] == mIdxPicoPions[idxPion2]) continue;

                    StHFTriplet *triplet = new StHFTriplet(pion1, kaon, pion2, mHFCuts->getHypotheticalMass(StHFCuts::kPion), mHFCuts->getHypotheticalMass(StHFCuts::kKaon), mHFCuts->getHypotheticalMass(StHFCuts::kPion), mIdxPicoPions[idxPion1], mIdxPicoKaons[idxKaon], mIdxPicoPions[idxPion2], mPrimVtx, mBField);

                    Float_t ispion2TOFmatched = 0;
                    Float_t ispion2BEMCmatched = 0;

                    if (mHFCuts->isTOFmatched(pion2)) ispion2TOFmatched = 1;
                    if (mHFCuts->isBEMCmatched(pion2)) ispion2BEMCmatched = 1;

                    /*Float_t hotSpot2 = 0;
                    if (mHFCuts->checkHotSpot(&mPrimVtx)) hotSpot2 = 1;
                    float tripletM = triplet->m();
                    float triplet_pair = tripletM - pairM;
                    */

                    const int nNtVars3 = ntp_DstarMeson_RightSign->GetNvar();
                    float ntVar3[nNtVars3];
                    int iii = 0;

                    /*ntVar3[iii++] = mPicoEvent->grefMult();
                    ntVar3[iii++] = mPicoEvent->refMult();
                    ntVar3[iii++] = mPicoEvent->runId();
                    ntVar3[iii++] = mPicoEvent->eventId();
                    ntVar3[iii++] = mPicoEvent->ZDCx();
                    ntVar3[iii++] = mPicoEvent->BBCx();
                    ntVar3[iii++] = hotSpot2;

                    ntVar3[iii++] = primary2;
                    ntVar3[iii++] = nGoodTracks;
                    */

                    ntVar3[iii++] = mPicoEvent->runId();
                    ntVar3[iii++] = pion1->pPt();
                    ntVar3[iii++] = pion1->pPtot();
                    //////////ntVar3[iii++] = triplet->particle1Dca();
                    ntVar3[iii++] = pion1->nSigmaPion();
                    //////////ntVar3[iii++] = (float)((float)pion1->nHitsFit()/(float)pion1->nHitsMax());
                    //////////ntVar3[iii++] = pion1->nHitsFit();
                    ntVar3[iii++] = mHFCuts->getnSigmaTOF(pion1, StPicoCutsBase::kPion);
                    ntVar3[iii++] = mHFCuts->getOneOverBeta(pion1,mHFCuts->getTofBetaBase(pion1),StPicoCutsBase::kPion);
                    ntVar3[iii++] = pion1->charge();
                    ntVar3[iii++] = ispion1TOFmatched;
                    ntVar3[iii++] = ispion1BEMCmatched;
                    //////////ntVar3[iii++] = pion1->pMom().PseudoRapidity();

                    ntVar3[iii++] = kaon->pPt();
                    ntVar3[iii++] = kaon->pPtot();
                    //////////ntVar3[iii++] = triplet->particle2Dca();
                    ntVar3[iii++] = kaon->nSigmaKaon();
                    //////////ntVar3[iii++] = (float)((float)kaon->nHitsFit()/(float)kaon->nHitsMax());
                    //////////ntVar3[iii++] = kaon->nHitsFit();
                    ntVar3[iii++] = mHFCuts->getnSigmaTOF(kaon, StPicoCutsBase::kKaon);
                    ntVar3[iii++] = mHFCuts->getOneOverBeta(kaon,mHFCuts->getTofBetaBase(kaon),StPicoCutsBase::kKaon);
                    ntVar3[iii++] = kaon->charge();
                    ntVar3[iii++] = iskaonTOFmatched;
                    ntVar3[iii++] = iskaonBEMCmatched;
                    //////////ntVar3[iii++] = kaon->pMom().PseudoRapidity();


                    ntVar3[iii++] = pion2->pPt();
                    ntVar3[iii++] = pion2->pPtot();
                    //////////ntVar3[iii++] = triplet->particle3Dca();
                    ntVar3[iii++] = pion2->nSigmaPion();
                    //////////ntVar3[iii++] = (float)((float)pion2->nHitsFit()/(float)pion2->nHitsMax());
                    //////////ntVar3[iii++] = pion2->nHitsFit();
                    ntVar3[iii++] = mHFCuts->getnSigmaTOF(pion2, StPicoCutsBase::kPion);
                    ntVar3[iii++] = mHFCuts->getOneOverBeta(pion2,mHFCuts->getTofBetaBase(pion2),StPicoCutsBase::kPion);
                    ntVar3[iii++] = pion2->charge();
                    ntVar3[iii++] = ispion2TOFmatched;
                    ntVar3[iii++] = ispion2BEMCmatched;
                    //////////ntVar3[iii++] = pion2->pMom().PseudoRapidity();

                    //////////ntVar3[iii++] = mPrimVtx.z();
                    //////////ntVar3[iii++] = mPicoEvent->vzVpd();
                    //////////ntVar3[iii++] = sqrt(pow(mPrimVtx.x(),2) + pow(mPrimVtx.y(),2));

                    ntVar3[iii++] = pair->m();
                    ntVar3[iii++] = pair->pt();
                    ntVar3[iii++] = pair->rapidity();
                    ntVar3[iii++] = pair->cosThetaStar();
                    ntVar3[iii++] = pair->pointingAngle();

                    ntVar3[iii++] = triplet->pt();
                    ntVar3[iii++] = triplet->m();
                    ntVar3[iii++] = triplet->rapidity();
                    //////////ntVar3[iii++] = triplet_pair;

                    ntVar3[iii++] = mPicoEvent->refMult();
                    ntVar3[iii++] = kaon->nHitsFit();
                    ntVar3[iii++] = pion1->nHitsFit();
                    ntVar3[iii++] = pion2->nHitsFit();
                    ntVar3[iii++] = kaon->nHitsDedx();
                    ntVar3[iii++] = pion1->nHitsDedx();
                    ntVar3[iii++] = pion2->nHitsDedx();
                    ntVar3[iii++] = mHFCuts->getbtofYLocal(kaon);
                    ntVar3[iii++] = mHFCuts->getbtofYLocal(pion1);
                    ntVar3[iii++] = mHFCuts->getbtofYLocal(pion2);

                    ntVar3[iii++] = mPrimVtx.z();
                    ntVar3[iii++] = mPicoEvent->vzVpd();

                    ntVar3[iii++] = triplet->particle2Dca();
                    ntVar3[iii++] = triplet->particle1Dca();
                    ntVar3[iii++] = triplet->particle3Dca();

                    
                    if ((pion1->charge() == pion2->charge())){
                        ntp_DstarMeson_SideBand->Fill(ntVar3);
                    } 

                }
                   ////else{continue;}
            }
            } 
              
        if (charge == -2 || charge == 2) {
                ntp_DMeson_LikeSign->Fill(ntVar);
                //////////vect_daughter_id.push_back(daughter_id);
            }
        
             } // for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon)
    } // for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1)
    
    } // if (mIdxPicoPions.size() > 0 && mIdxPicoKaons.size() > 0)

    /*/////
    if (nD0 > 0) {
        hNTracksRemoved->Fill(tracksToRemove.size());
        hNTracksGoodToFit->Fill(nGoodTracks);
        hNTracksPrimary->Fill(nPrimary);
        hNTracksDiffRemovedPrimary->Fill(nPrimary-tracksToRemove.size());

       ///// cout << nTracks << "   " << tracksToRemove.size() << endl;
        hNTracksDiffRemovedGlobal->Fill(nTracks-tracksToRemove.size());
        if(mPrimVtx.x() > -0.25 && mPrimVtx.x() < -0.16 && mPrimVtx.y() > -0.25 && mPrimVtx.y() < -0.16) hHotSpotDiffRemovedPrimary->Fill(nPrimary-tracksToRemove.size());
    }

   hD0VsRemoved->Fill(tracksToRemove.size(), nD0);
   *///////

    mIdxPicoPions.clear();
    mIdxPicoPions.shrink_to_fit();

    primaryTracks.clear();
    primaryTracks.shrink_to_fit();

    mIdxPicoKaons.clear();
    mIdxPicoKaons.shrink_to_fit();

    tracksToRemove.clear();
    tracksToRemove.shrink_to_fit();

        } //----end loop for betterevents

    
    } // end loop for 2X matching


    } // end loop for vzmax cut

   
   
    return kStOK;
}

//____________________________________________________________________________________
/*TVector3 StPicoD0AnaMaker::refitVertex(bool always){
    bool pairRem=true;
    bool singleTrack=true;
//    bool singleTrack=!pairRem;
    std::vector<int> goodTracksToFit;

    TH1F *hPVDiffX = static_cast<TH1F*>(mOutList->FindObject("hPVDiffX"));
    TH1F *hPVDiffY = static_cast<TH1F*>(mOutList->FindObject("hPVDiffY"));
    TH1F *hPVDiffZ = static_cast<TH1F*>(mOutList->FindObject("hPVDiffZ"));
    TH1F *hPVDiffXRemoved = static_cast<TH1F*>(mOutList->FindObject("hPVDiffXRemoved"));
    TH1F *hPVDiffYRemoved = static_cast<TH1F*>(mOutList->FindObject("hPVDiffYRemoved"));
    TH1F *hPVDiffZRemoved = static_cast<TH1F*>(mOutList->FindObject("hPVDiffZRemoved"));

    TH1F *hRemovedPairMass = static_cast<TH1F*>(mOutList->FindObject("hRemovedPairMass"));
    TH1F *hInvMass[3] = {static_cast<TH1F*>(mOutList->FindObject("hInvMassBDT12")), static_cast<TH1F*>(mOutList->FindObject("hInvMassBDT23")), static_cast<TH1F*>(mOutList->FindObject("hInvMassBDT35"))};

    bool isRemovedtrack=false;

    //removing according single track cuts

    if (singleTrack) {
        float dca;
        for (unsigned short iTrack = 0; iTrack < primaryTracks.size(); ++iTrack) {
            StPicoTrack* trk = mPicoDst->track(primaryTracks[iTrack]);
                dca = (mPrimVtx - trk->origin()).Mag();
                if (dca > 0.009) tracksToRemove.push_back(primaryTracks[iTrack]);
        }
    }


    //removing with pair cuts
    if (pairRem) {
        for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
            StPicoTrack const *pion1 = mPicoDst->track(mIdxPicoPions[idxPion1]);
            for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
                StPicoTrack const *kaon = mPicoDst->track(mIdxPicoKaons[idxKaon]);

                if(kaon->charge()+pion1->charge() != 0) continue;
                StHFPair *pair = new StHFPair(pion1, kaon, mHFCuts->getHypotheticalMass(StPicoCutsBase::kPion), mHFCuts->getHypotheticalMass(StPicoCutsBase::kKaon), mIdxPicoPions[idxPion1], mIdxPicoKaons[idxKaon], mPrimVtx, mBField, kTRUE);
                if(!(pair->m()>1.81 && pair->m()<1.92)) continue;

                //CUTS

                if (cos(pair->pointingAngle()) > 0.9 && pair->dcaDaughters() < 0.007 && pair->m()>1.81 && pair->m()<1.92) {
                    tracksToRemove.push_back(mIdxPicoPions[idxPion1]);
                    tracksToRemove.push_back(mIdxPicoKaons[idxKaon]);
                    hRemovedPairMass->Fill(pair->m());
                    isRemovedtrack = true;
                }


                //BDT
                //find the correct pT bin

                int pTbin = 0;
                for (int pT = 0; pT < nptBins; pT++) {
                    if(pair->pt() >= momBins[pT] && pair->pt() < momBins[pT+1]) pTbin = pT;
                }

                k_dca[pTbin] = pair->particle2Dca();
                pi1_dca[pTbin] = pair->particle1Dca();
                D_decayL[pTbin] = pair->decayLength();
                cosTheta[pTbin] = cos(pair->pointingAngle());
                dcaD0ToPv[pTbin] = pair->DcaToPrimaryVertex();
                dcaDaughters[pTbin] = pair->dcaDaughters();
//                thetaStar[pTbin] = pair->cosThetaStar();

                //evaluate BDT, continue just pairs that have passed BDT cut
                float valueMVA = reader[pTbin]->EvaluateMVA("BDT method");
                if(valueMVA > bdtCuts[pTbin]) {
                    //filling plots of invariant mass for unlike and like sign pairs
                    hInvMass[pTbin]->Fill(pair->m());
                    //excluding daughter tracks
                    tracksToRemove.push_back(mIdxPicoPions[idxPion1]);
                    tracksToRemove.push_back(mIdxPicoKaons[idxKaon]);

                }
            }
        }
    }

    for (int iTrk = 0; iTrk < primaryTracks.size(); ++iTrk) {
        if(std::binary_search(tracksToRemove.begin(), tracksToRemove.end(), primaryTracks[iTrk])) continue;
        goodTracksToFit.push_back(primaryTracks[iTrk]);
    }

    TVector3 newKFVertex(mPrimVtx.x(), mPrimVtx.y(), mPrimVtx.z());
    if (goodTracksToFit.size()>2) {
//    if (isRemovedtrack || always) {
        //Make new vertex and evaluate stuff:
        StPicoKFVertexFitter kfVertexFitter;
//        KFVertex kfVertex = kfVertexFitter.primaryVertexRefit(mPicoDst, tracksToRemove); //when removed tracks need to be checked
        KFVertex kfVertex = kfVertexFitter.primaryVertexRefitUsingTracks(mPicoDst, goodTracksToFit); //when you have array with tracks for refit
        if (kfVertex.GetX()) {
            newKFVertex.SetXYZ(kfVertex.GetX(), kfVertex.GetY(), kfVertex.GetZ());
            hPVDiffX->Fill(mPrimVtx.x() - newKFVertex.x());
            hPVDiffY->Fill(mPrimVtx.y() - newKFVertex.y());
            hPVDiffZ->Fill(mPrimVtx.z() - newKFVertex.z());

            if (tracksToRemove.size() > 0) {
                hPVDiffXRemoved->Fill(mPrimVtx.x() - newKFVertex.x());
                hPVDiffYRemoved->Fill(mPrimVtx.y() - newKFVertex.y());
                hPVDiffZRemoved->Fill(mPrimVtx.z() - newKFVertex.z());
            }
        }
    }
    nGoodTracks=goodTracksToFit.size();

    goodTracksToFit.clear();
    goodTracksToFit.shrink_to_fit();

    return newKFVertex;
}*/

//____________________________________________________________________________________
int StPicoD0AnaMaker::analyzeCandidates() {
    
    /*/////////
    TH1D *h_tracktest = static_cast<TH1D*>(mOutList->FindObject("h_tracktest"));
    TH1D *h_tracktest_TOF = static_cast<TH1D*>(mOutList->FindObject("h_tracktest_TOF"));

    /////const char *aNames[]   = {"all", "NHitFit>20", "pT>0.6", "dca<1.5","TPC pion", "TPC kaon"};
    /////const char *aNamesTOF[]   = {"all TOF match", "NHitFit>20 & TOF", "pT>0.6 & TOF","dca<1.5 & TOF","TPC & TOF pion", "TPC & TOF kaon"};
    
    const char *aNames[]   = {"all", "NHitFit>20", "pT>0.16", "dca<1.5","TPC pion", "TPC kaon"};
    const char *aNamesFast[]   = {"all TOF/BEMC match", "NHitFit>20 & TOF/BEMC", "pT>0.16 & TOF/BEMC","dca<1.5 & TOF/BEMC","TPC & TOF/BEMC pion", "TPC & TOF/BEMC kaon"};
    for (unsigned int ii = 0; ii < 6; ii++) {
        h_tracktest->GetXaxis()->SetBinLabel(ii+1, aNames[ii]);
        h_tracktest_Fast->GetXaxis()->SetBinLabel(ii+1, aNamesTOF[ii]);
    }

    for(unsigned int i=0;i<mPicoDst->numberOfTracks();i++) {
        StPicoTrack const *t = mPicoDst->track(i);
        if (!t) continue;
        if (!t->isHFTTrack()) continue;
        h_tracktest->Fill(1);
        if (mHFCuts->isGoodTrack(t)) h_tracktest->Fill(2); // NhitsFit

        float pt=t->pPt();
        if ((pt>0.16) && mHFCuts->isGoodTrack(t)) h_tracktest->Fill(3);
        float dca = (mPrimVtx - t->origin()).Mag();
        if (dca<1.5 && (pt>0.16) && mHFCuts->isGoodTrack(t)) h_tracktest->Fill(4);
        bool tpcPion = mHFCuts->isTPCHadron(t, StPicoCutsBase::kPion);
        bool tpcKaon = mHFCuts->isTPCHadron(t, StPicoCutsBase::kKaon);
        if ((pt>0.16) && (dca<1.5 ) && (mHFCuts->isGoodTrack(t))) {
            if (tpcPion) h_tracktest->Fill(5);
           if (tpcKaon) h_tracktest->Fill(6);
        }
        if (mHFCuts->isTOFmatched(t)||mHFCuts->isBEMCmatched(t)) {
            h_tracktest_Fast->Fill(1);
            if (mHFCuts->isGoodTrack(t)) h_tracktest_TOF->Fill(2); // NhitsFit
            if (pt>0.16 && mHFCuts->isGoodTrack(t)) h_tracktest_TOF->Fill(3);
            if (dca<1.5 && pt>0.16 && mHFCuts->isGoodTrack(t)) h_tracktest_TOF->Fill(4);
            if ((pt>0.16) && (dca<1.5 ) && (mHFCuts->isGoodTrack(t))) {
//                if (tpcPion && mHFCuts->isTOFHadronPID(t, mHFCuts->getTofBetaBase(t), StPicoCutsBase::kPion))  h_tracktest_TOF->Fill(5);
                if (tpcPion)  h_tracktest_TOF->Fill(5);
                if (tpcKaon && mHFCuts->isTOFHadronPID(t, mHFCuts->getTofBetaBase(t), StPicoCutsBase::kKaon))  h_tracktest_TOF->Fill(6);
            }
        }


    }
//    for (unsigned short idxPion1 = 0; idxPion1 < mIdxPicoPions.size(); ++idxPion1) {
//        StPicoTrack const *t = mPicoDst->track(mIdxPicoPions[idxPion1]);
//        ntp->Fill(t->pPt(), t->pMom().phi(), t->pMom().pseudoRapidity(), t->nSigmaPion(), t->nHitsFit(), getOneOverBeta(t, mHFCuts->getTofBetaBase(t), StPicoCutsBase::kPion), mPicoEvent->eventId(), mPicoEvent->runId());
//    }
//
//    for (unsigned short idxKaon = 0; idxKaon < mIdxPicoKaons.size(); ++idxKaon) {
//        StPicoTrack const *t = mPicoDst->track(mIdxPicoKaons[idxKaon]);
//        ntp_kaon->Fill(t->pPt(), t->pMom().phi(), t->pMom().pseudoRapidity(), t->nSigmaKaon(), t->nHitsFit(), getOneOverBeta(t, mHFCuts->getTofBetaBase(t), StPicoCutsBase::kKaon), mPicoEvent->eventId(), mPicoEvent->runId());
//    }
     *//////////
    return kStOK;
}
    

bool StPicoD0AnaMaker::isVectorInVectorOfVectors(const std::vector<std::vector<float>> vectorOfVectors, const std::vector<float> targetVector) {
    return false;
    for (const auto vector : vectorOfVectors) {
        if (vector == targetVector) {
            return true;
    }
    }
}

int StPicoD0AnaMaker::ConvertRunID(int number) {
    std::ifstream file("/gpfs/mnt/gpfs01/star/pwg/palsp/DMaker_pp500_r17/picoLists/runs_numbers.list");
    if (!file.is_open()) {
        std::cerr << "Failed to open file:    " << "runs_numbers.list" << std::endl;
        return -1;
    }

    std::string line;
    int lineNumber = 1;
    while (std::getline(file, line)) {
        int value = std::stoi(line);
        if (value == number) {
            file.close();
            return lineNumber;
        }
        lineNumber++;
    }

    file.close();
    return -1; // Number not found
}