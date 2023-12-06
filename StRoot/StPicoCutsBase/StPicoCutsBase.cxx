#include <limits>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>

#include "StPicoCutsBase.h"

///#include "StPicoEvent/StPicoPhysicalHelix.h"*.root
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "phys_constants.h"
#include "SystemOfUnits.h"

#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"

ClassImp(StPicoCutsBase)

// _________________________________________________________
StPicoCutsBase::StPicoCutsBase() : TNamed("PicoCutsBase", "PicoCutsBase"),
                                   mTOFCorr(new StV0TofCorrection), mPicoDst(NULL), mEventStatMax(7), mTOFResolution(0.013),
                                   mBadRunListFileName("picoList_bad.list"), mnMatchedFast(2),mVzMax(50.), mVzVpdVzMax(6.),
                                   mNHitsFitMin(20), mRequireHFT(false), mNHitsFitnHitsMax(0.52), mPrimaryDCAtoVtxMax(6.0), mPtMin(0.2), mEtaMax(1.0), mHybridTof(false), mHybridTofKaon(false), mHybridTofPion(false), mHybridTofBetterBetaCuts(false), mHybridTofBetterBetaCutsKaon(false), mHybridTofBetterBetaCutsPion(false), mTPCBetterCutsPion(false), mOnlyHotSpot(false) {

    for (Int_t idx = 0; idx < kPicoPIDMax; ++idx) {
        mPtRange[idx][0] = std::numeric_limits<float>::lowest();
        mPtRange[idx][1] = std::numeric_limits<float>::max();
        mDcaMin[idx] = std::numeric_limits<float>::lowest();
        mDcaMinTertiary[idx] = std::numeric_limits<float>::lowest();
        mPtotRangeTOF[idx][0] = std::numeric_limits<float>::lowest();
        mPtotRangeTOF[idx][1] = std::numeric_limits<float>::max();
        mPtotRangeHybridTOF[idx][0] = std::numeric_limits<float>::lowest();
        mPtotRangeHybridTOF[idx][1] = std::numeric_limits<float>::max();
        mTPCNSigmaMax[idx] = std::numeric_limits<float>::max();
        mTPCNSigmaMin[idx] = std::numeric_limits<float>::lowest();
        mTOFNSigmaMax[idx] = std::numeric_limits<float>::max();
        mTOFNSigmaMin[idx] = std::numeric_limits<float>::lowest();
        mTOFDeltaOneOverBetaMax[idx] = std::numeric_limits<float>::max();
        mTOFDeltaOneOverBetaMin[idx] = std::numeric_limits<float>::lowest();
    }

    mHypotheticalMass[kPion]      = M_PION_PLUS;
    mHypotheticalMass2[kPion]     = M_PION_PLUS*M_PION_PLUS;
    mHypotheticalMass[kKaon]      = M_KAON_PLUS;
    mHypotheticalMass2[kKaon]     = M_KAON_PLUS*M_KAON_PLUS;
    mHypotheticalMass[kProton]    = M_PROTON;
    mHypotheticalMass2[kProton]   = M_PROTON*M_PROTON;
    mHypotheticalMass[kElectron]  = M_ELECTRON;
    mHypotheticalMass2[kElectron] = M_ELECTRON*M_ELECTRON;
    mHypotheticalMass[kMuon]      = M_MUON_PLUS;
    mHypotheticalMass2[kMuon]     = M_MUON_PLUS*M_MUON_PLUS;
    mHypotheticalMass[kK0Short]   = M_KAON_0_SHORT;
    mHypotheticalMass2[kK0Short]  = M_KAON_0_SHORT*M_KAON_0_SHORT;
    mHypotheticalMass[kLambda]    = M_LAMBDA;
    mHypotheticalMass2[kLambda]   = M_LAMBDA*M_LAMBDA;
}

// _________________________________________________________
StPicoCutsBase::StPicoCutsBase(const Char_t *name) : TNamed(name, name),
                                                     mTOFCorr(new StV0TofCorrection), mPicoDst(NULL), mEventStatMax(7), mTOFResolution(0.013),
                                                     mBadRunListFileName("picoList_bad_MB.list"), mVzMax(50.), mVzVpdVzMax(6.),
                                                     mNHitsFitMin(20), mRequireHFT(true), mNHitsFitnHitsMax(0.52), mPrimaryDCAtoVtxMax(6.0), mPtMin(0.2), mEtaMax(1.0), mHybridTof(false),mHybridTofKaon(false), mHybridTofPion(false), mHybridTofBetterBetaCuts(false), mHybridTofBetterBetaCutsKaon(false), mHybridTofBetterBetaCutsPion(false), mTPCBetterCutsPion(false), mOnlyHotSpot(false) {
    // -- constructor

    for (Int_t idx = 0; idx < kPicoPIDMax; ++idx) {
        mPtRange[idx][0] = std::numeric_limits<float>::lowest();
        mPtRange[idx][1] = std::numeric_limits<float>::max();
        mDcaMin[idx] = std::numeric_limits<float>::lowest();
        mDcaMinTertiary[idx] = std::numeric_limits<float>::lowest();
        mPtotRangeTOF[idx][0] = std::numeric_limits<float>::lowest();
        mPtotRangeTOF[idx][1] = std::numeric_limits<float>::max();
        mPtotRangeHybridTOF[idx][0] = std::numeric_limits<float>::lowest();
        mPtotRangeHybridTOF[idx][1] = std::numeric_limits<float>::max();
        mTPCNSigmaMax[idx] = std::numeric_limits<float>::max();
        mTPCNSigmaMin[idx] = std::numeric_limits<float>::lowest();
        mTOFNSigmaMax[idx] = std::numeric_limits<float>::max();
        mTOFNSigmaMin[idx] = std::numeric_limits<float>::lowest();
        mTOFDeltaOneOverBetaMax[idx] = std::numeric_limits<float>::max();
        mTOFDeltaOneOverBetaMin[idx] = std::numeric_limits<float>::lowest();
    }

    mHypotheticalMass[kPion]      = M_PION_PLUS;
    mHypotheticalMass2[kPion]     = M_PION_PLUS*M_PION_PLUS;
    mHypotheticalMass[kKaon]      = M_KAON_PLUS;
    mHypotheticalMass2[kKaon]     = M_KAON_PLUS*M_KAON_PLUS;
    mHypotheticalMass[kProton]    = M_PROTON;
    mHypotheticalMass2[kProton]   = M_PROTON*M_PROTON;
    mHypotheticalMass[kElectron]  = M_ELECTRON;
    mHypotheticalMass2[kElectron] = M_ELECTRON*M_ELECTRON;
    mHypotheticalMass[kMuon]      = M_MUON_PLUS;
    mHypotheticalMass2[kMuon]     = M_MUON_PLUS*M_MUON_PLUS;
    mHypotheticalMass[kK0Short]   = M_KAON_0_SHORT;
    mHypotheticalMass2[kK0Short]  = M_KAON_0_SHORT*M_KAON_0_SHORT;
    mHypotheticalMass[kLambda]    = M_LAMBDA;
    mHypotheticalMass2[kLambda]   = M_LAMBDA*M_LAMBDA;

}
// _________________________________________________________
StPicoCutsBase::~StPicoCutsBase() {
    // destructor

    if (mTOFCorr)
        delete mTOFCorr;
    mTOFCorr = NULL;
}

// _________________________________________________________
void StPicoCutsBase::initBase() {
    // -- init cuts class

    // -- Read in bad run list and fill vector
    // -----------------------------------------

    // -- open list
    ifstream runs;

    // -- open in working dir
    runs.open(mBadRunListFileName.Data());
    if (!runs.is_open()) {
        runs.open(Form("picoLists/%s", mBadRunListFileName.Data()));
        if (!runs.is_open()) {
            cout << "StPicoCutsBase::initBase -- Bad run list NOT found :" << mBadRunListFileName << endl;
            cout << "StPicoCutsBase::initBase -- continue without bad run selection! " << endl;
            //exit(EXIT_FAILURE);
        }
    }

    if (runs.is_open()) {
        Int_t runId = 0;
        while( runs >> runId )
            mVecBadRunList.push_back(runId);

        runs.close();

        // -- sort bad runs vector
        std::sort(mVecBadRunList.begin(), mVecBadRunList.end());
    }
}

// _________________________________________________________
bool StPicoCutsBase::isGoodEvent(StPicoDst const * const picoDst, int *aEventNums) {
    
    // -- set current mPicoDst
    mPicoDst = picoDst;

    // -- get picoDst event
    StPicoEvent* picoEvent = mPicoDst->event();

    // -- set current primary vertex
    mPrimVtx = picoEvent->primaryVertex();

    //    if(mOnlyHotSpot) cout<<"m hot spor true"<<endl;
    if(mOnlyHotSpot && !(checkHotSpot(&mPrimVtx))) return false;

    // -- quick method without providing stats
    if (!aEventNums) {
        return (isGoodRun(picoEvent) && isGoodTrigger(picoEvent));
    }

    // -- reset event cuts
    for (unsigned int ii = 0; ii < mEventStatMax; ++ii)
        aEventNums[ii] = 0;

    unsigned int iCut = 0;
    // -- 0 - before event cuts
    aEventNums[iCut] = 0;

    // -- 1 - is bad run
    ++iCut;
    if (!isGoodRun(picoEvent)) aEventNums[iCut] = 1;

    // -- 2 - No Trigger fired
    ++iCut;
    //  trigger - is ok?
    if (!isGoodTrigger(picoEvent)) aEventNums[iCut] = 1;

    ++iCut;

    //if the event is wrong, the array member is 1 (eccept [0])
    // -- is rejected
    bool isGoodEvent = true;
    for (unsigned int ii = 0; ii < mEventStatMax-1; ++ii) {
        if (aEventNums[ii]) isGoodEvent = false;
    }

    if(!isGoodEvent) aEventNums[mEventStatMax-1]=1;

    return isGoodEvent;
}

//_________________________________________________________________________-

bool StPicoCutsBase::isBetterEvent(StPicoDst const * const picoDst, int *aEventCuts) {
    
    // -- set current mPicoDst
    mPicoDst = picoDst;

    // -- get picoDst event
    StPicoEvent* picoEvent = mPicoDst->event();

    // -- set current primary vertex
    mPrimVtx = picoEvent->primaryVertex();

//    if(mOnlyHotSpot) cout<<"m hot spor true"<<endl;
    if(mOnlyHotSpot && !(checkHotSpot(&mPrimVtx))) return false;

    // -- quick method without providing stats
   /* if (!aEventCuts) {
        return (isGoodRun(picoEvent) && isGoodTrigger(picoEvent) &&
                (fabs(picoEvent->primaryVertex().z()) < mVzMax) &&
                (fabs(picoEvent->primaryVertex().z() - picoEvent->vzVpd()) < mVzVpdVzMax) && isMatchedFast(mPicoDst));
    }*/
    if (!aEventCuts) {
        return (isGoodRun(picoEvent) && isGoodTrigger(picoEvent) &&
                (fabs(picoEvent->primaryVertex().z()) < mVzMax) &&
                (fabs(picoEvent->primaryVertex().z() - picoEvent->vzVpd()) < mVzVpdVzMax) && (sqrt(pow(picoEvent->primaryVertex().x(),2) + pow(picoEvent->primaryVertex().y(),2)) < mVrMax));
    }

    // -- reset event cuts
    for (unsigned int ii = 0; ii < mEventStatMax; ++ii)
        aEventCuts[ii] = 0;

    unsigned int iCut = 0;
    // -- 0 - before event cuts
    aEventCuts[iCut] = 0;

    // -- 1 - is bad run
    ++iCut;
    if (!isGoodRun(picoEvent)) aEventCuts[iCut] = 1;

    // -- 2 - No Trigger fired
    ++iCut;
//  trigger - is ok?
    if (!isGoodTrigger(picoEvent)) aEventCuts[iCut] = 1;

    // -- 3 - Vertex z outside cut window
    ++iCut;
    if (fabs(picoEvent->primaryVertex().z()) >= mVzMax) aEventCuts[iCut] = 1;

    // -- 4 Vertex z - vertex_z(vpd) outside cut window
    ++iCut;
    if (fabs(picoEvent->primaryVertex().z() - picoEvent->vzVpd()) >= mVzVpdVzMax) aEventCuts[iCut] = 1;

    // -- 5 Vertex r - vertex_r(vpd) outside cut window
    ++iCut;
    if (sqrt(pow(picoEvent->primaryVertex().x(),2) + pow(picoEvent->primaryVertex().y(),2)) >= mVrMax) aEventCuts[iCut] = 1;

   /* // -- 5 - No. of matched tracks in Fast detectors
    ++iCut;
    if (!isMatchedFast(mPicoDst)) aEventCuts[iCut] = 1;
   */

    ++iCut;

    //if the event is wrong, the array member is 1 (eccept [0])
    // -- is rejected
    bool isBetterEvent = true;
    for (unsigned int ii = 0; ii < mEventStatMax-1; ++ii) {
        if (aEventCuts[ii]) isBetterEvent = false;
    }

    if(!isBetterEvent) aEventCuts[mEventStatMax-1]=1;

    return isBetterEvent;
}

// _________________________________________________________
bool StPicoCutsBase::isGoodRun(StPicoEvent const * const picoEvent) const {
    // -- is good run (not in bad runlist)
    return (!(std::binary_search(mVecBadRunList.begin(), mVecBadRunList.end(), picoEvent->runId())));
}

// _________________________________________________________
bool StPicoCutsBase::isGoodTrigger(StPicoEvent const * const picoEvent) const {
    // -- is good trigger in list of good triggerIds

    for(std::vector<unsigned int>::const_iterator iter = mVecTriggerIdList.begin(); iter != mVecTriggerIdList.end(); ++iter)
        if(picoEvent->isTrigger(*iter))
            return true;

    return false;
}

//_______________________________________________________________

bool StPicoCutsBase::isMatchedFast(StPicoDst const * const picoDst) {
mPicoDst = picoDst;
UInt_t nTracks = mPicoDst->numberOfTracks();
Int_t nMatchedFast = 0;
   for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
    StPicoTrack* trk = mPicoDst->track(iTrack);
    if (trk->isPrimary()){
        if (isTOFmatched(trk) || isBEMCmatched(trk)) nMatchedFast += 1;
    }
   }

   return (nMatchedFast >= mnMatchedFast);
}
// _________________________________________________________
bool StPicoCutsBase::isGoodTrack(StPicoTrack const * const trk) const {
    //  int tofIndex = trk->bTofPidTraitsIndex();
//  bool TofMatch = kFALSE;
//  StPicoBTofPidTraits* tofPidTraits;
//  if (tofIndex >= 0)  tofPidTraits = mPicoDst->btofPidTraits(tofIndex);
//  if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
    return (trk->nHitsFit() >= mNHitsFitMin && (((float)trk->nHitsFit())/(trk->nHitsMax())) > mNHitsFitnHitsMax && cutMaxDcaToPrimVertex(trk) && trk->pPt() > mPtMin && trk->pMom().PseudoRapidity() < mEtaMax && trk->pMom().PseudoRapidity() > mEtaMin);
//    return ((!mRequireHFT || trk->isHFTTrack()) && trk->nHitsFit() >= mNHitsFitMin && cutMaxDcaToPrimVertex(trk) && trk->gPt() > mPtMin);

}

// _________________________________________________________
bool StPicoCutsBase::isPionTPC(StPicoTrack const *trk) const {
    if (!isGoodTrack(trk)) return false;
    ///////if (!cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
    if (!isTPCPion(trk)) return false;
    return true;
}

// _________________________________________________________
bool StPicoCutsBase::isKaonTPC(StPicoTrack const *trk) const {
    if (!isGoodTrack(trk)) return false;
    ///////if (!cutMinDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
    if (!isTPCKaon(trk)) return false;
    return true;
}

// __________________________________________________________

bool StPicoCutsBase::isPionTOF(StPicoTrack const *trk) const {
    if (!isGoodTrack(trk)) return false;
    ///////if (!cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
    ///////if (!isTPCPion(trk)) return false;
    bool tofWtpc = false;
    double pt = trk->pPt();
    if (mHybridTofWithBEMC) {

         if(isTOFmatched(trk)) tofWtpc = isTOFPionBetterCuts(trk);
         else if (TMath::Abs(trk->nSigmaPion())<4) tofWtpc = true;
         /*else if((pt > 1.6) && isBEMCmatched(trk)) {
            if (TMath::Abs(trk->nSigmaPion())<2) tofWtpc = true;
         }*/
        /*if (ptot < 1.3) {
            if(isTOFmatched(trk)) tofWtpc = isTOFPionBetterCuts(trk);
            if(!isTOFmatched(trk)&&(isBEMCmatched(trk))) tofWtpc = true;
        }

        if (ptot > 1.3 && ptot < 2.07) {
            if(isTOFmatched(trk)) tofWtpc = isTOFPion(trk);
            if(!isTOFmatched(trk)&&(isBEMCmatched(trk))) tofWtpc = true;
        }

        if (ptot > 2.07) {
            if(isTOFmatched(trk)||(isBEMCmatched(trk))) tofWtpc = isTOFPion(trk);

        }*/
    }
    if (!mHybridTofWithBEMC) {
        if (mHybridTofBetterBetaCutsPion) {
           /* if (mHybridTofPion) tofWtpc = isHybridTOFPionBetterCuts(trk);
            if (!mHybridTofPion)*/
        tofWtpc = isTOFPionBetterCuts(trk);
        }


        if (!mHybridTofBetterBetaCutsPion) { // Original constant cut
           /* if (mHybridTofPion) tofWtpc = isHybridTOFPion(trk);
            if (!mHybridTofPion)*/
            tofWtpc = isTOFPion(trk); //actual entrypoint
        }
    }
    return tofWtpc;
}

bool StPicoCutsBase::isKaonTOF(StPicoTrack const *const trk) const {
    if (!isGoodTrack(trk)) return false;
    //////if (!cutMinDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
    //////if (!isTPCKaon(trk)) return false;
    bool tofWtpc = false;
    ///////double ptot = trk->gPtot();
    double pt = trk->pPt();
    if (mHybridTofWithBEMC) {
        if(isTOFmatched(trk)) tofWtpc = isTOFKaonBetterCuts(trk);
        else if (TMath::Abs(trk->nSigmaKaon())<4) tofWtpc = true;
         /*else if((pt > 1.6) && isBEMCmatched(trk)) {
            if (TMath::Abs(trk->nSigmaKaon())<2) tofWtpc = true;
         }*/
        
        /*
        if (ptot < 1.3) {
            tof = isTOFKaonBetterCuts(trk);
        }

        if (ptot > 1.3 && ptot < 2.07) {
            if(isTOFmatched(trk)) tof = isTOFKaonBetterCuts(trk);
            if(!isTOFmatched(trk)&&(isBEMCmatched(trk))) tof = true;
        }

        if (ptot > 2.07) {
            if(isTOFmatched(trk)||(isBEMCmatched(trk))) tof = isTOFKaon(trk);
        }*/



    }
    if (!mHybridTofWithBEMC) {

        if (mHybridTofBetterBetaCutsKaon) {
            /*if (mHybridTofKaon) tofWtpc = isHybridTOFKaonBetterCuts(trk);
            if (!mHybridTofKaon)*/
             tofWtpc = isTOFKaonBetterCuts(trk); //actual entrypoint
        }


        if (!mHybridTofBetterBetaCutsKaon) { // Original constant cut
            /*if (mHybridTofKaon) tofWtpc = isHybridTOFKaon(trk);
            if (!mHybridTofKaon)*/
             tofWtpc = isTOFKaon(trk);
        }
    }
    return tofWtpc;
}

// _________________________________________________________
bool StPicoCutsBase::isGoodPion(StPicoTrack const *const trk) const {
    ///return true;
    if (isGoodTrack(trk) && isPionTOF(trk) && isPionTPC(trk)) return true;
    else return false;
    ///////if (!cutMinDcaToPrimVertex(trk, StPicoCutsBase::kPion)) return false;
    ///if (!isPionTOF(trk)) return false;
    ////if (!isPionTPC(trk)) return false;
    /*bool tofWtpc = false;
    double pt = trk->pPt();
    if (mHybridTofWithBEMC) {

         if(isTOFmatched(trk)) tofWtpc = isTOFPion(trk);
         else if((pt > 1.6) && isBEMCmatched(trk)) {
            if (TMath::Abs(trk->nSigmaPion())<2) tofWtpc = true;
         }
        ////if (ptot < 1.3) {
        ////    if(isTOFmatched(trk)) tofWtpc = isTOFPionBetterCuts(trk);
         ////   if(!isTOFmatched(trk)&&(isBEMCmatched(trk))) tofWtpc = true;
        ////}

        ////if (ptot > 1.3 && ptot < 2.07) {
        ////    if(isTOFmatched(trk)) tofWtpc = isTOFPion(trk);
        ////    if(!isTOFmatched(trk)&&(isBEMCmatched(trk))) tofWtpc = true;
        ////}

        ////if (ptot > 2.07) {
            if(isTOFmatched(trk)||(isBEMCmatched(trk))) tofWtpc = isTOFPion(trk);

        ///}
    }
    if (!mHybridTofWithBEMC) {
        if (mHybridTofBetterBetaCutsPion) {
            /////if (mHybridTofPion) tofWtpc = isHybridTOFPionBetterCuts(trk);
            ////if (!mHybridTofPion)
            tofWtpc = isTOFPionBetterCuts(trk);
        }


        if (!mHybridTofBetterBetaCutsPion) { // Original constant cut
            ///if (mHybridTofPion) tofWtpc = isHybridTOFPion(trk);
            ////if (!mHybridTofPion)
            tofWtpc = isTOFPion(trk); //actual entrypoint
        }
    }
    return tofWtpc;*/
}

// _________________________________________________________
/*old bool StPicoCutsBase::isTOFPionCutOK(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const {
    if (tofBeta <= 0) {return false;}
    //////double ptot    = trk->gPtot();
    bool tofpion = false;
    double ptot    = trk->pPtot();
    float betaInv = ptot / sqrt(ptot*ptot + mHypotheticalMass2[pidFlag]);
  /*  float pion_higher = 6-8/3*ptot;
    float pion_lower = -6+8/3*ptot;


    if(ptot<1.5) {
        return ((1 / tofBeta - 1 / betaInv) / 0.011 < pion_higher && (1 / tofBeta - 1 / betaInv) / 0.011 > pion_lower);
    }
    if(ptot>1.5) {
        return ((1 / tofBeta - 1 / betaInv) / 0.011 < 2 && (1 / tofBeta - 1 / betaInv) / 0.011 > -2);
    }
    float pion_higher = 6-2*ptot;
    float pion_lower = -6+2*ptot;


    if(ptot<1.5) {
        tofpion =  ((1 / tofBeta - 1 / betaInv) / 0.011 < pion_higher && (1 / tofBeta - 1 / betaInv) / 0.011 > pion_lower);
    }
    if(ptot>1.5) {
        tofpion =  ((1 / tofBeta - 1 / betaInv) / 0.011 < 3 && (1 / tofBeta - 1 / betaInv) / 0.011 > -3);
    }

    return tofpion;


    }
    */
    
bool StPicoCutsBase::isTOFPionCutOK(StPicoTrack const *trk) const {
    
    if (!(trk->isTofTrack())) return false;
    float nSigma = std::numeric_limits<float>::quiet_NaN();
    int bTOFPidTraitsIndex = trk->bTofPidTraitsIndex();
    StPicoBTofPidTraits* tofPidTraits = mPicoDst->btofPidTraits(bTOFPidTraitsIndex);  
    
    bool tofpion = false;
    double ptot    = trk->pPtot();
    float pion_higher = 12-(5*ptot);
    float pion_lower = -12+(5*ptot);

    nSigma = tofPidTraits->nSigmaPion();


    if(ptot<1.2) {
        tofpion =  (nSigma < pion_higher && nSigma > pion_lower);
    }
    if(ptot>1.2) {
        tofpion =  (nSigma < 6 && nSigma > -6);
    }

    return tofpion; 

}

bool StPicoCutsBase::isTPCPionBetterCut(StPicoTrack const *trk) const {
    
    /////if (!(trk->isTofTrack())) return false;
    float nSigma = std::numeric_limits<float>::quiet_NaN();
    /////int bTOFPidTraitsIndex = trk->bTofPidTraitsIndex();
    /////StPicoBTofPidTraits* tofPidTraits = mPicoDst->btofPidTraits(bTOFPidTraitsIndex);  
    
    bool tpcpion = false;
    double ptot    = trk->pPtot();

    /////double ptot    = trk->pPtot();
    float f_res =  0.9548 + (0.0002/(pow((ptot + 0.1046), 3.9041))); //sigma
    float f_pos =  0.0988 + (-0.0012/(pow((ptot - 0.0706), 3.1974)));  //mean 

    /////float pion_higher = 2.92*f_res + f_pos;///done
    /////float pion_lower = -2.70*f_res + f_pos;///done

    float pion_higher = 6.0*f_res + f_pos; ///loose
    float pion_lower = -6.0*f_res + f_pos;///loose

    nSigma = trk->nSigmaPion();

    tpcpion = (nSigma < pion_higher && nSigma > pion_lower);
    /*if(ptot<1.5) {
        tofpion =  (nSigma < pion_higher && nSigma > pion_lower);
    }
    if(ptot>1.5) {
        tofpion =  (nSigma < 3 && nSigma > -3);
    }*/

    return tpcpion; 

}
// _________________________________________________________
/* old bool StPicoCutsBase::isTOFBetterPion(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const {
    return isTOFPionCutOK(trk, tofBeta, pidFlag);
}
*/

bool StPicoCutsBase::isTOFBetterPion(StPicoTrack const *trk) const {
    return isTOFPionCutOK(trk);

}

// _________________________________________________________
/*
bool StPicoCutsBase::isHybridTOFBetterPion(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const {
    if (tofBeta <= 0 || tofBeta != tofBeta )
        return true;

    return isTOFPionCutOK(trk, tofBeta, pidFlag);

}
*/

// __________________________________________________________

// _________________________________________________________
bool StPicoCutsBase::isGoodKaon(StPicoTrack const *const trk) const {
    if (isGoodTrack(trk) && isKaonTOF(trk) && isKaonTPC(trk)) return true;
    else return false;
    //////if (!cutMinDcaToPrimVertex(trk, StPicoCutsBase::kKaon)) return false;
    ////if (!isKaonTOF(trk)) return false;
    ////if (!isKaonTPC(trk)) return false;
    /*bool tofWtpc = false;
    ///////double ptot = trk->gPtot();
    double pt = trk->pPt();
    if (mHybridTofWithBEMC) {
        if(isTOFmatched(trk)) tofWtpc = isTOFKaonBetterCuts(trk);
         else if((pt > 1.6) && isBEMCmatched(trk)) {
            if (TMath::Abs(trk->nSigmaKaon())<2) tofWtpc = true;
         }
        
        
        ////if (ptot < 1.3) {
        ////    tof = isTOFKaonBetterCuts(trk);
        ////}

        ////if (ptot > 1.3 && ptot < 2.07) {
        ////   if(isTOFmatched(trk)) tof = isTOFKaonBetterCuts(trk);
         /////   if(!isTOFmatched(trk)&&(isBEMCmatched(trk))) tof = true;
        ////}

        /////if (ptot > 2.07) {
        /////   if(isTOFmatched(trk)||(isBEMCmatched(trk))) tof = isTOFKaon(trk);
        ////}



    }
    if (!mHybridTofWithBEMC) {

        if (mHybridTofBetterBetaCutsKaon) {
            //// if (mHybridTofKaon) tofWtpc = isHybridTOFKaonBetterCuts(trk);
            ////if (!mHybridTofKaon) 
            tofWtpc = isTOFKaonBetterCuts(trk); //actual entrypoint
        }


        if (!mHybridTofBetterBetaCutsKaon) { // Original constant cut
            ////if (mHybridTofKaon) tofWtpc = isHybridTOFKaon(trk);
            ////if (!mHybridTofKaon)
            tofWtpc = isTOFKaon(trk);
        }
    }
    return tofWtpc;*/
}

// _________________________________________________________
/* old bool StPicoCutsBase::isTOFKaonCutOK(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const {
    if (tofBeta <= 0) {return false;}
    /////double ptot    = trk->gPtot();
    double ptot    = trk->pPtot();
    float betaInv = ptot / sqrt(ptot*ptot + mHypotheticalMass2[pidFlag]);

    /////float f_res = 1.12190e+00 + 1.04744e-01/pow((ptot - 1.41968e-01),1.37291e+00);  //sigma
    /////float f_pos = 7.84921e-02 + 4.22205e-04/pow((ptot + 6.17544e-02),6.95657e+00);  //mean

    float f_res =  1.07078 + 0.0194723/pow((ptot + 0.0183138), 3.41348);  //sigma
    float f_pos = 0.0283416 + 0.00230815/pow((ptot + 0.0270787), 4.96298);  //mean


    float kaon_higher = 1.8*f_res + f_pos;
    float kaon_lower = -1.2*f_res + f_pos;


    return ( (1/tofBeta - 1/betaInv)/0.012 < kaon_higher && (1/tofBeta - 1/betaInv)/0.012 > kaon_lower );
}
*/

bool StPicoCutsBase::isTOFKaonCutOK(StPicoTrack const *trk) const {
    if (!(trk->isTofTrack())) return false;
    float nSigma = std::numeric_limits<float>::quiet_NaN();
    int bTOFPidTraitsIndex = trk->bTofPidTraitsIndex();
    StPicoBTofPidTraits* tofPidTraits = mPicoDst->btofPidTraits(bTOFPidTraitsIndex);   
    
    double ptot    = trk->pPtot();
    float f_res =  1.21 + (0.10/(pow((ptot - 0.17), 1.21))); //sigma
    float f_pos =  -0.01 + (0.02/(pow((ptot - 0.17), 1.73)));  //mean

   /* float mSigmalower = -2.31;
    if (ptot > 1.25 && ptot < 1.50) mSigmalower = -1.56;
    if (ptot > 1.50 && ptot < 1.65) mSigmalower = -1.18;
    if (ptot > 1.65) mSigmalower = -0.8;

    float mSigmahigher = 3.18;
    if (ptot > 2.0) mSigmahigher = 2.41;

    float kaon_higher = mSigmahigher*f_res + f_pos;
    float kaon_lower = mSigmalower*f_res + f_pos; *//////done

    float kaon_higher = 6.0*f_res + f_pos;
    float kaon_lower = -6.0*f_res + f_pos;
    /*float kaon_higher = 2.41*f_res + f_pos;
    float kaon_lower = -1.18*f_res + f_pos;*/


    nSigma = tofPidTraits->nSigmaKaon();
    
    return (nSigma < kaon_higher && nSigma > kaon_lower);

}

// _________________________________________________________
 /*old bool StPicoCutsBase::isTOFBetterKaon(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const {
    return isTOFKaonCutOK(trk, tofBeta, pidFlag);

}
*/

bool StPicoCutsBase::isTOFBetterKaon(StPicoTrack const *trk) const {
    return isTOFKaonCutOK(trk);

}

// _________________________________________________________
/*old bool StPicoCutsBase::isHybridTOFBetterKaon(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const {
    if (tofBeta <= 0 || tofBeta != tofBeta )
        return true;

    return isTOFKaonCutOK(trk, tofBeta, pidFlag);

}
*/


// _________________________________________________________
bool StPicoCutsBase::isGoodProton(StPicoTrack const *const trk) const {
    if (!isGoodTrack(trk)) return false;
    if (!cutMinDcaToPrimVertex(trk, StPicoCutsBase::kProton)) return false;
    if (!isTPCProton(trk)) return false;
    bool tof = false;
    if (mHybridTof) tof = isHybridTOFProton(trk);
    if (!mHybridTof) tof = isTOFProton(trk);

    return tof;
}

// _________________________________________________________
bool StPicoCutsBase::cutMinDcaToPrimVertex(StPicoTrack const * const trk, int pidFlag) const {
    // -- check on min dca for identified particle
    float dca = (mPrimVtx - trk->origin()).Mag();
    return (dca >= mDcaMin[pidFlag]);
}

// _________________________________________________________
bool StPicoCutsBase::cutMaxDcaToPrimVertex(StPicoTrack const * const trk) const {
    // -- check on max dca for all particles
    float dca = (mPrimVtx - trk->origin()).Mag();
    return (dca < mPrimaryDCAtoVtxMax);
}

// _________________________________________________________
bool StPicoCutsBase::cutMinDcaToPrimVertexTertiary(StPicoTrack const * const trk, int pidFlag) const {
    // -- check on min dca for identified particle - used for tertiary particles only

    StPicoPhysicalHelix helix = trk->helix(mPicoDst->event()->bField());
    helix.moveOrigin(helix.pathLength(mPrimVtx));
    float dca = (mPrimVtx - helix.origin()).Mag();

    return (dca >= mDcaMinTertiary[pidFlag]);
}

// _________________________________________________________
bool StPicoCutsBase::isTPCHadron(StPicoTrack const * const trk, int pidFlag) const {
    // -- check for good hadron in TPC
    float nSigma = std::numeric_limits<float>::quiet_NaN();

    if (pidFlag == kPion)
        nSigma = trk->nSigmaPion();
    else if (pidFlag == kKaon)
        nSigma = trk->nSigmaKaon();
    else if (pidFlag == kProton)
        nSigma = trk->nSigmaProton();
     
    return (nSigma < mTPCNSigmaMax[pidFlag] && nSigma > mTPCNSigmaMin[pidFlag] && trk->nHitsFit() >= mNHitsFitMin);
}



// _________________________________________________________
bool StPicoCutsBase::isTOFHadronPID(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const {
    if (tofBeta <= 0) {return false;}
    /////double ptot    = trk->gPtot();
    double ptot    = trk->pPtot();
    float betaInv = ptot / sqrt(ptot*ptot + mHypotheticalMass2[pidFlag]);
    double deltabetainv = (1/tofBeta - 1/betaInv);
    return (deltabetainv > mTOFDeltaOneOverBetaMin[pidFlag] && mTOFDeltaOneOverBetaMax[pidFlag] > deltabetainv);
}



// _________________________________________________________
 /*Old bool StPicoCutsBase::isTOFHadron(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const {
    // -- check for good hadron in TOF in ptot range
    //    use for
    //      - primary hadrons
    //      - secondarys from charm decays (as an approximation)
    //    return:
    //      not in ptot range : true

    // -- only apply, if in ptot range


//  float ptot = trk->gPtot();
//  if (ptot < mPtotRangeTOF[pidFlag][0] || ptot >= mPtotRangeTOF[pidFlag][1])
//    return true;

    return isTOFHadronPID(trk, tofBeta, pidFlag);
}*/

// _________________________________________________________

bool StPicoCutsBase::isTOFHadron(StPicoTrack const * const trk, int pidFlag) const {
    // -- check for good hadron in TOF
    if (!(trk->isTofTrack())) return false;
    float nSigma = std::numeric_limits<float>::quiet_NaN();
    int bTOFPidTraitsIndex = trk->bTofPidTraitsIndex();
    StPicoBTofPidTraits* tofPidTraits = mPicoDst->btofPidTraits(bTOFPidTraitsIndex);    

    if (pidFlag == kPion)
        nSigma = tofPidTraits->nSigmaPion();
    else if (pidFlag == kKaon)
        nSigma = tofPidTraits->nSigmaKaon();
    else if (pidFlag == kProton)
        nSigma = tofPidTraits->nSigmaProton();
     
    return (nSigma < mTOFNSigmaMax[pidFlag] && nSigma > mTOFNSigmaMin[pidFlag]);
}


// _________________________________________________________
bool StPicoCutsBase::isTOFmatched(StPicoTrack const *trk) const {
    /*OLD
    int tofIndex = trk->bTofPidTraitsIndex();
    trk->isTofTrack();
    bool TofMatch = kFALSE;
    StPicoBTofPidTraits* tofPidTraits;
    if (tofIndex >= 0)  tofPidTraits = mPicoDst->btofPidTraits(tofIndex);
    if (tofIndex >= 0 && tofPidTraits && tofPidTraits->btofMatchFlag() > 0)  TofMatch = kTRUE;
     */
    return trk->isTofTrack();
}

// _________________________________________________________
bool StPicoCutsBase::isBEMCmatched(StPicoTrack const *trk) const {

    return trk->isBemcTrack();
}

// _________________________________________________________
bool StPicoCutsBase::isHybridTOFHadron(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const {
    // -- check for good hadron in TOF in ptot range
    //    use for
    //      - primary hadrons
    //      - secondarys from charm decays (as an approximation)
    //    return:
    //      not in ptot range : true
    //      no TOF info       : true

    // -- only apply, if in ptot range
//  float ptot = trk->gPtot();
//  if (ptot < mPtotRangeHybridTOF[pidFlag][0] || ptot >= mPtotRangeHybridTOF[pidFlag][1])
//    return true;

    // -- only apply, if has TOF information
    if (tofBeta <= 0 || tofBeta != tofBeta )
        return true;

    return isTOFHadronPID(trk, tofBeta, pidFlag);
}

// _________________________________________________________
StPicoBTofPidTraits* StPicoCutsBase::hasTofPid(StPicoTrack const * const trk) const {
    // -- check if track has TOF pid information
    //    return NULL otherwise

    int index2tof = trk->bTofPidTraitsIndex();
    return (index2tof >= 0) ? mPicoDst->btofPidTraits(index2tof) : NULL;
}

// _________________________________________________________
float StPicoCutsBase::getTofBetaBase(StPicoTrack const * const trk) const {
    int index2tof = trk->bTofPidTraitsIndex(); //if smaller than 0 => not TOF track
    float beta = std::numeric_limits<float>::quiet_NaN();

    if(index2tof >= 0) {
        StPicoBTofPidTraits *tofPid = mPicoDst->btofPidTraits(index2tof);
        if(tofPid)  beta = tofPid->btofBeta();

        if (beta < 1e-4) {
            TVector3 const btofHitPos = tofPid->btofHitPos();
            StPicoPhysicalHelix helix = trk->helix(mPicoDst->event()->bField());

            float L = tofPathLength(&mPrimVtx, &btofHitPos, helix.curvature());
            float tof = tofPid->btof();
            if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
            else beta = std::numeric_limits<float>::quiet_NaN();
        }
    }

    return beta;
}

// _________________________________________________________
float StPicoCutsBase::getOneOverBeta(StPicoTrack const * const trk,  float const & tofBeta, int pidFlag) const {
    if ((tofBeta <= 0) || (tofBeta!=tofBeta))
        return std::numeric_limits<float>::quiet_NaN();
    else {
        float m2 = mHypotheticalMass[pidFlag]*mHypotheticalMass[pidFlag];
    /////float ptot = trk->gPtot();
    float ptot = trk->pPtot();
    float betaInv = ptot / sqrt(ptot*ptot + m2);
    return (1/tofBeta - 1/betaInv);
    }
}

// _________________________________________________________

float StPicoCutsBase::getnSigmaTOF(StPicoTrack const * const trk, int pidFlag) const {
    int index2tof = trk->bTofPidTraitsIndex(); //if smaller than 0 => not TOF track
    float nSigmaTOF = std::numeric_limits<float>::quiet_NaN();

    if (index2tof >= 0) {
        StPicoBTofPidTraits *tofPid = mPicoDst->btofPidTraits(index2tof);
        if(tofPid) {
            if (pidFlag == kPion) nSigmaTOF = tofPid->nSigmaPion();
            if (pidFlag == kKaon) nSigmaTOF = tofPid->nSigmaKaon();
            if (pidFlag == kProton) nSigmaTOF = tofPid->nSigmaProton();
        } 
    }

    return nSigmaTOF;

}

// _________________________________________________________

float StPicoCutsBase::getTofBeta(StPicoTrack const * const trk) const {
    return getTofBetaBase(trk);
}

// _________________________________________________________
float StPicoCutsBase::getTofBeta(StPicoTrack const * const trk,
                                 TVector3 const & secondaryMother, TVector3 const & secondaryVtx) const {
    // -- provide correced beta of TOF for pico track
    //    use for
    //      - secondaries

    float beta = std::numeric_limits<float>::quiet_NaN();

    StPicoBTofPidTraits *tofPid = hasTofPid(trk);
    if (!tofPid)
        return beta;
//
//
//  // -- set waypoints
//  mTOFCorr->setVectors3D(mPrimVtx)(secondaryVtx)(tofHit);
//
//  // -- set mother track
//  mTOFCorr->setMotherTracks(secondaryMother);
//
//  float tof = tofPid->btof();
//  StPicoPhysicalHelix helix = trk->helix(mPicoDst->event()->bField());
//
//  // -- correct beta
//  mTOFCorr->correctBeta(helix, tof, beta);
//
//  // -- clean up
//  mTOFCorr->clearContainers();
//
    return beta;
}

// _________________________________________________________
float StPicoCutsBase::getTofBeta(StPicoTrack const * const trk,
                                 TVector3 const & secondaryMother, TVector3 const & secondaryVtx,
                                 TVector3 const & tertiaryMother,  TVector3 const & tertiaryVtx) const {
    // -- provide correced beta of TOF for pico track
    //    use for
    //      - tertiaries

    float beta = std::numeric_limits<float>::quiet_NaN();

//  StPicoBTofPidTraits *tofPid = hasTofPid(trk);
//  if (!tofPid)
//    return beta;
//
//  StThreeVectorD tofHit = tofPid->btofHitPos();
//
//  // -- set waypoints
//  mTOFCorr->setVectors3D(mPrimVtx)(secondaryVtx)(tertiaryVtx)(tofHit);
//
//  // -- set mother track
//  mTOFCorr->setMotherTracks(secondaryMother)(tertiaryMother);
//
//  float tof = tofPid->btof();
//  StPicoPhysicalHelix helix = trk->helix(mPicoDst->event()->bField());
//
//  // -- correct beta
////  mTOFCorr->correctBeta(helix, tof, beta);
//
//  // -- clean up
//  mTOFCorr->clearContainers();

    return beta;
}

// _________________________________________________________
float StPicoCutsBase::tofPathLength(const TVector3* beginPoint, const TVector3* endPoint, float curvature) const {
    float xdif =  endPoint->x() - beginPoint->x();
    float ydif =  endPoint->y() - beginPoint->y();

    float C = sqrt(xdif*xdif + ydif*ydif);
    float s_perp = C;
    if (curvature){
        float R = 1/curvature;
        s_perp = 2*R * asin(C/(2*R)); //arc length
    }

    float s_z = fabs(endPoint->z() - beginPoint->z());
    float value = sqrt(s_perp*s_perp + s_z*s_z);

    return value;
}

// _________________________________________________________
bool StPicoCutsBase::checkHotSpot(TVector3* vertex) const {
    if(vertex->x()>-0.25 && vertex->x()<-0.16 && vertex->y()>-0.25 && vertex->y()<-0.16) return true;
    else return false;
}
