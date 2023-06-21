#ifndef STPICOCUTSBASE_H
#define STPICOCUTSBASE_H
#include <limits>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>

#include "TNamed.h"
#include "TString.h"
#include "TVector3.h"

#include "StBTofUtil/StV0TofCorrection.h"

class StPicoTrack;
class StPicoEvent;
class StPicoDst;
class StPicoBTofPidTraits;

class StPicoCutsBase : public TNamed
{
public:

    StPicoCutsBase();
    StPicoCutsBase(const Char_t *name);
    ~StPicoCutsBase();

    void initBase();

    virtual void init() { initBase(); }

    bool isGoodEvent(StPicoDst const * const picoDst, int *aEventNums = NULL);
    bool isBetterEvent(StPicoDst const * const picoDst, int *aEventCuts = NULL);
    bool isGoodRun(StPicoEvent const * const picoEvent) const;
    bool isGoodTrigger(StPicoEvent const * const picoEvent) const;
    bool isMatchedFast(StPicoDst const * const picoDst);
    bool isGoodTrack(StPicoTrack const * const trk) const;
    bool isGoodPion(StPicoTrack const * const trk) const;
    bool isGoodKaon(StPicoTrack const * const trk) const;
    bool isGoodProton(StPicoTrack const * const trk) const;
    bool checkHotSpot(TVector3*) const;
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- DCA to Primary vertex
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    bool cutMinDcaToPrimVertex(StPicoTrack const * const trk, int pidFlag) const;
    bool cutMinDcaToPrimVertexTertiary(StPicoTrack const * const trk, int pidFlag) const;
    bool cutMaxDcaToPrimVertex(StPicoTrack const * const trk) const;

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- PID
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    enum ePicoPID {kPion, kKaon, kProton, kElectron, kMuon, kK0Short, kLambda, kPicoPIDMax};

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- TOF PID
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    StPicoBTofPidTraits* hasTofPid(StPicoTrack const * const trk) const;

    bool isTPCHadron(StPicoTrack const * const trk, int pidFlag) const;
    bool isTPCPion(StPicoTrack const *trk) const;
    bool isTPCKaon(StPicoTrack const *trk) const;
    bool isTPCProton(StPicoTrack const *trk) const;

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- TOF PID
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    bool isTOFHadronPID(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const;

    bool isTOFHadron(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const;
    bool isHybridTOFHadron(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const;

    // -- Is TOF particle in ptot range
    //    if track has no TOF information - return false
    //    use for
    //      - primary hadrons
    //      - secondarys from charm decays (as an approximation)


    bool isTOFmatched(StPicoTrack const *trk) const;
    bool isTOFPion(StPicoTrack const *trk) const;
    bool isTOFKaon(StPicoTrack const *trk) const;
    bool isTOFProton(StPicoTrack const *trk) const;

    bool isTOFPion(StPicoTrack const *trk,   float const & tofBeta) const;
    bool isTOFKaon(StPicoTrack const *trk,   float const & tofBeta) const;
    bool isTOFProton(StPicoTrack const *trk, float const & tofBeta) const;

    bool isBEMCmatched(StPicoTrack const *trk) const;

    bool isPionTPC(StPicoTrack const *trk) const;
    bool isKaonTPC(StPicoTrack const *trk) const;
    bool isPionTOF(StPicoTrack const *trk) const;
    bool isKaonTOF(StPicoTrack const *trk) const;

    // -- Is TOF particle in ptot range
    //    if track has no TOF information - return true
    //    use for
    //      - primary hadrons
    //      - secondarys from charm decays (as an approximation)
    bool isHybridTOFPion(StPicoTrack const *trk) const;
    bool isHybridTOFKaon(StPicoTrack const *trk) const;
    bool isHybridTOFProton(StPicoTrack const *trk) const;

    bool isHybridTOFPion(StPicoTrack const *trk,   float const & tofBeta) const;
    bool isHybridTOFKaon(StPicoTrack const *trk,   float const & tofBeta) const;
    bool isHybridTOFProton(StPicoTrack const *trk, float const & tofBeta) const;


    // Comparing kaon parameters with functions and not constants

    bool isTOFKaonCutOK(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const;
    bool isTOFBetterKaon(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const;
    bool isHybridTOFBetterKaon(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const;
    bool isTOFKaonBetterCuts(StPicoTrack const *trk) const;
    bool isTOFKaonBetterCuts(StPicoTrack const *trk, float const & tofBeta) const;
    bool isHybridTOFKaonBetterCuts(StPicoTrack const *trk) const;
    bool isHybridTOFKaonBetterCuts(StPicoTrack const *trk, float const & tofBeta) const;


    bool isTOFPionCutOK(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const;
    bool isTOFBetterPion(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const;
    bool isHybridTOFBetterPion(StPicoTrack const *trk, float const & tofBeta, int pidFlag) const;
    bool isTOFPionBetterCuts(StPicoTrack const *trk) const;
    bool isTOFPionBetterCuts(StPicoTrack const *trk, float const & tofBeta) const;
    bool isHybridTOFPionBetterCuts(StPicoTrack const *trk) const;
    bool isHybridTOFPionBetterCuts(StPicoTrack const *trk, float const & tofBeta) const;


    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    const unsigned int&  eventStatMax()  const { return mEventStatMax; }

    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    // -- SETTER for CUTS
    // -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

    void setBadRunListFileName(const char* fileName);
    
    void addTriggerId(unsigned int triggerId);
    void addTriggerByName(std::string name);
    void showTriggers();
    
    void setnMatchedFast(int i);
    void setCutVzMax(float f);
    void setCutVzVpdVzMax(float f);
    void setCheckHotSpot(bool b);

    void setCutNHitsFitMin(int i);
    void setCutRequireHFT(bool b);
    void setCutNHitsFitnHitsMax(float f);

    void setCutPrimaryDCAtoVtxMax(float f);
    void setCutPtRange(float min, float max, int pidFlag);
    void setCutPtMin(float min);
    void setCutEtaMax(float max);
    void setCutEtaMin(float min);
    void setCutDcaMin(float min, int pidFlag);
    void setCutDcaMinTertiary(float min, int pidFlag);
    void setCutTPCNSigmaMax(float f, int pidFlag);
    void setCutTPCNSigmaMin(float f, int pidFlag);
    void setCutTOFNSigmaMax(float f, int pidFlag);
    void setCutTOFNSigmaMin(float f, int pidFlag);
    void setCutTOFDeltaOneOverBeta(float f, int pidFlag);
    void setCutPtotRangeTOF(float min, float max, int pidFlag);
    void setCutPtotRangeHybridTOF(float min, float max, int pidFlag);

    void setCutPionPtRange(float min, float max);
    void setCutPionDcaMin(float min);
    void setCutPionDcaMinTertiary(float min);
    void setCutTPCNSigmaPionMax(float f);
    void setCutTPCNSigmaPionMin(float f);
    void setCutTOFNSigmaPionMax(float f);
    void setCutTOFNSigmaPionMin(float f);
    void setCutTOFDeltaOneOverBetaPion(float f);
    void setCutPionPtotRangeTOF(float min, float max);
    void setCutPionPtotRangeHybridTOF(float min, float max);

    void setCutKaonPtRange(float min, float max);
    void setCutKaonDcaMin(float min);
    void setCutKaonDcaMinTertiary(float min);
    void setCutTPCNSigmaKaonMax(float f);
    void setCutTPCNSigmaKaonMin(float f);
    void setCutTOFNSigmaKaonMax(float f);
    void setCutTOFNSigmaKaonMin(float f);
    void setCutTOFDeltaOneOverBetaKaon(float f);
    void setCutKaonPtotRangeTOF(float min, float max);
    void setCutKaonPtotRangeHybridTOF(float min, float max);

    void setCutProtonPtRange(float min, float max);
    void setCutProtonDcaMin(float min);
    void setCutProtonDcaMinTertiary(float min);
    void setCutTPCNSigmaProtonMax(float f);
    void setCutTPCNSigmaProtonMin(float f);
    void setCutTOFNSigmaProtonMax(float f);
    void setCutTOFNSigmaProtonMin(float f);
    void setCutTOFDeltaOneOverBetaProton(float f);
    void setCutProtonPtotRangeTOF(float min, float max);
    void setCutProtonPtotRangeHybridTOF(float min, float max);

    void setHybridTof(bool t);
    void setHybridTofKaon(bool t);
    void setHybridTofPion(bool t);

    void setHybridTofBetterBetaCuts(bool t);
    void setHybridTofBetterBetaCutsKaon(bool t);
    void setHybridTofBetterBetaCutsPion(bool t);

    void setHybridTofWithBEMC(bool t);

    float getCutVzMax();
    float getCutVzVpdVzMax();

    float tofPathLength(const TVector3* beginPoint, const TVector3* endPoint, float curvature) const;

    // -- calculate beta of track -- basic calculation
    float getTofBetaBase(StPicoTrack const* const trk) const;
    float getOneOverBeta(StPicoTrack const * const trk,  float const & tofBeta, int pidFlag) const;

    // -- calculate beta of track -- for primary particles
    float getTofBeta(StPicoTrack const* const trk) const;

    // -- calculate corrected beta of track -- for secondary particles
    float getTofBeta(StPicoTrack const * const trk,
                     TVector3 const & secondaryMother, TVector3 const & secondaryVtx) const;

    // -- calculate corrected beta of track -- for tertiary particles
    float getTofBeta(StPicoTrack const * const trk,
                     TVector3 const & secondaryMother, TVector3 const & secondaryVtx,
                     TVector3 const & tertiaryMother,  TVector3 const & tertiaryVtx) const;

    const float& getHypotheticalMass(int pidFlag)           const;

private:

    StPicoCutsBase(StPicoCutsBase const &);
    StPicoCutsBase& operator=(StPicoCutsBase const &);

    StV0TofCorrection* mTOFCorr;  // TOF correction

    TVector3    mPrimVtx;   // primary vertex of current event
    const StPicoDst*  mPicoDst;   //! ptr to picoDst

    unsigned int mEventStatMax;   // number of event cuts

    float        mTOFResolution;  // TOF resolution = 0.013

    // -- bad run list
    TString mBadRunListFileName;
    std::vector<int> mVecBadRunList;

    // -- trigger id list
    std::vector<unsigned int> mVecTriggerIdList;

    // -- event cuts
    int   mnMatchedFast;
    float mVzMax;
    float mVzVpdVzMax;

    // -- tracking
    int   mNHitsFitMin;
    bool  mRequireHFT;
    float mNHitsFitnHitsMax;
    float mPrimaryDCAtoVtxMax;         // used for primary selection for TOF Beta recalculation
    float mPtMin;
    float mEtaMax;
    float mEtaMin;
    bool  mHybridTof;
    bool  mHybridTofKaon;
    bool  mHybridTofPion;
    bool  mHybridTofBetterBetaCuts;
    bool  mHybridTofBetterBetaCutsKaon;
    bool  mHybridTofBetterBetaCutsPion;
    bool  mHybridTofWithBEMC;
    bool  mOnlyHotSpot;


    // -- acceptance - per particle type [ePicoPID]
    float mPtRange[kPicoPIDMax][2];

    // -- dca to primary vertex - per particle type [ePicoPID]
    float mDcaMin[kPicoPIDMax];
    float mDcaMinTertiary[kPicoPIDMax];

    // -- PID cuts - per particle type [ePicoPID]
    float mHypotheticalMass[kPicoPIDMax];        // hypothetical mass
    float mHypotheticalMass2[kPicoPIDMax];       // hypothetical mass squared

    float mTPCNSigmaMax[kPicoPIDMax];
    float mTPCNSigmaMin[kPicoPIDMax];
    float mTOFDeltaOneOverBetaMax[kPicoPIDMax];
    float mTOFDeltaOneOverBetaMin[kPicoPIDMax];

    float mPtotRangeTOF[kPicoPIDMax][2];         // momentum range [min,max], where TOF PID is applied
    float mPtotRangeHybridTOF[kPicoPIDMax][2];   // momentum range [min,max], where Hybrid TOF PID is applied

    ClassDef(StPicoCutsBase,1)
};

inline void StPicoCutsBase::setBadRunListFileName(const char* fileName) { mBadRunListFileName = fileName; }
inline void StPicoCutsBase::addTriggerId(unsigned int triggerId) {mVecTriggerIdList.push_back(triggerId);}

inline void StPicoCutsBase::setnMatchedFast(int i)              { mnMatchedFast    = i; }
inline void StPicoCutsBase::setCutVzMax(float f)              { mVzMax            = f; }
inline void StPicoCutsBase::setCutVzVpdVzMax(float f)         { mVzVpdVzMax       = f; }

inline float StPicoCutsBase::getCutVzMax()              { return mVzMax; }
inline float StPicoCutsBase::getCutVzVpdVzMax()         { return mVzVpdVzMax; }

inline void StPicoCutsBase::setCheckHotSpot(bool b)         { mOnlyHotSpot       = b; }

inline void StPicoCutsBase::setCutNHitsFitMin(int i)          { mNHitsFitMin      = i; }
inline void StPicoCutsBase::setCutRequireHFT(bool b)          { mRequireHFT       = b; }
inline void StPicoCutsBase::setCutNHitsFitnHitsMax(float f)   { mNHitsFitnHitsMax = f; }

inline void StPicoCutsBase::setCutPrimaryDCAtoVtxMax(float f) { mPrimaryDCAtoVtxMax = f; }

inline void StPicoCutsBase::setCutPtRange(float min, float max, int pidFlag)            { mPtRange[pidFlag][0] = min;
    mPtRange[pidFlag][1] = max; }

inline void StPicoCutsBase::setCutPtMin(float min)            { mPtMin = min;}
inline void StPicoCutsBase::setCutEtaMax(float max)            { mEtaMax = max;}
inline void StPicoCutsBase::setCutEtaMin(float min)            { mEtaMin = min;}

inline void StPicoCutsBase::setCutDcaMin(float min, int pidFlag)                        { mDcaMin[pidFlag] = min; }
inline void StPicoCutsBase::setCutDcaMinTertiary(float min, int pidFlag)                { mDcaMinTertiary[pidFlag] = min; }

inline void StPicoCutsBase::setCutTPCNSigmaMax(float f, int pidFlag)                       { mTPCNSigmaMax[pidFlag] = f; }
inline void StPicoCutsBase::setCutTPCNSigmaMin(float f, int pidFlag)                       { mTPCNSigmaMin[pidFlag] = f; }
inline void StPicoCutsBase::setCutTOFNSigmaMax(float f, int pidFlag)     { 
    if (pidFlag == kPion) mTOFResolution = 0.011;
    if (pidFlag == kKaon) mTOFResolution = 0.012;
    mTOFDeltaOneOverBetaMax[pidFlag] = f*mTOFResolution;
    }
inline void StPicoCutsBase::setCutTOFNSigmaMin(float f, int pidFlag)     { 
    if (pidFlag == kPion) mTOFResolution = 0.011;
    if (pidFlag == kKaon) mTOFResolution = 0.012;
    mTOFDeltaOneOverBetaMin[pidFlag] = f*mTOFResolution;
    }
inline void StPicoCutsBase::setCutTOFDeltaOneOverBeta(float f, int pidFlag)             { mTOFDeltaOneOverBetaMax[pidFlag] = f;}
inline void StPicoCutsBase::setCutPtotRangeTOF(float min, float max, int pidFlag)       { mPtotRangeTOF[pidFlag][0] = min;
    mPtotRangeTOF[pidFlag][1] = max; }
inline void StPicoCutsBase::setCutPtotRangeHybridTOF(float min, float max, int pidFlag) { mPtotRangeHybridTOF[pidFlag][0] = min;
    mPtotRangeHybridTOF[pidFlag][1] = max; }

inline void StPicoCutsBase::setCutPionPtRange(float min, float max)              { setCutPtRange(min, max, StPicoCutsBase::kPion); }
inline void StPicoCutsBase::setCutPionDcaMin(float min)                          { setCutDcaMin(min, StPicoCutsBase::kPion); }
inline void StPicoCutsBase::setCutPionDcaMinTertiary(float min)                  { setCutDcaMinTertiary(min, StPicoCutsBase::kPion); }
inline void StPicoCutsBase::setCutTPCNSigmaPionMax(float f)                         { setCutTPCNSigmaMax(f, StPicoCutsBase::kPion); }
inline void StPicoCutsBase::setCutTPCNSigmaPionMin(float f)                         { setCutTPCNSigmaMin(f, StPicoCutsBase::kPion); }
inline void StPicoCutsBase::setCutTOFNSigmaPionMax(float f)                         { setCutTOFNSigmaMax(f, StPicoCutsBase::kPion); }
inline void StPicoCutsBase::setCutTOFNSigmaPionMin(float f)                         { setCutTOFNSigmaMin(f, StPicoCutsBase::kPion); }
inline void StPicoCutsBase::setCutTOFDeltaOneOverBetaPion(float f)               { setCutTOFDeltaOneOverBeta(f, StPicoCutsBase::kPion); }
inline void StPicoCutsBase::setCutPionPtotRangeTOF(float min, float max)         { setCutPtotRangeTOF(min, max, StPicoCutsBase::kPion); }
inline void StPicoCutsBase::setCutPionPtotRangeHybridTOF(float min, float max)   { setCutPtotRangeHybridTOF(min, max, StPicoCutsBase::kPion); }

inline void StPicoCutsBase::setCutKaonPtRange(float min, float max)              { setCutPtRange(min, max, StPicoCutsBase::kKaon); }
inline void StPicoCutsBase::setCutKaonDcaMin(float min)                          { setCutDcaMin(min, StPicoCutsBase::kKaon); }
inline void StPicoCutsBase::setCutKaonDcaMinTertiary(float min)                  { setCutDcaMinTertiary(min, StPicoCutsBase::kKaon); }
inline void StPicoCutsBase::setCutTPCNSigmaKaonMax(float f)                         { setCutTPCNSigmaMax(f, StPicoCutsBase::kKaon); }
inline void StPicoCutsBase::setCutTPCNSigmaKaonMin(float f)                         { setCutTPCNSigmaMin(f, StPicoCutsBase::kKaon); }
inline void StPicoCutsBase::setCutTOFNSigmaKaonMax(float f)                         { setCutTOFNSigmaMax(f, StPicoCutsBase::kKaon); }
inline void StPicoCutsBase::setCutTOFNSigmaKaonMin(float f)                         { setCutTOFNSigmaMin(f, StPicoCutsBase::kKaon); }
inline void StPicoCutsBase::setCutTOFDeltaOneOverBetaKaon(float f)               { setCutTOFDeltaOneOverBeta(f, StPicoCutsBase::kKaon); }
inline void StPicoCutsBase::setCutKaonPtotRangeTOF(float min, float max)         { setCutPtotRangeTOF(min, max, StPicoCutsBase::kKaon); }
inline void StPicoCutsBase::setCutKaonPtotRangeHybridTOF(float min, float max)   { setCutPtotRangeHybridTOF(min, max, StPicoCutsBase::kKaon); }

inline void StPicoCutsBase::setCutProtonPtRange(float min, float max)            { setCutPtRange(min, max, StPicoCutsBase::kProton); }
inline void StPicoCutsBase::setCutProtonDcaMin(float min)                        { setCutDcaMin(min, StPicoCutsBase::kProton); }
inline void StPicoCutsBase::setCutProtonDcaMinTertiary(float min)                { setCutDcaMinTertiary(min, StPicoCutsBase::kProton); }
inline void StPicoCutsBase::setCutTPCNSigmaProtonMax(float f)                       { setCutTPCNSigmaMax(f, StPicoCutsBase::kProton); }
inline void StPicoCutsBase::setCutTPCNSigmaProtonMin(float f)                       { setCutTPCNSigmaMin(f, StPicoCutsBase::kProton); }
inline void StPicoCutsBase::setCutTOFNSigmaProtonMax(float f)                       { setCutTOFNSigmaMax(f, StPicoCutsBase::kProton); }
inline void StPicoCutsBase::setCutTOFNSigmaProtonMin(float f)                       { setCutTOFNSigmaMin(f, StPicoCutsBase::kProton); }
inline void StPicoCutsBase::setCutTOFDeltaOneOverBetaProton(float f)             { setCutTOFDeltaOneOverBeta(f, StPicoCutsBase::kProton); }
inline void StPicoCutsBase::setCutProtonPtotRangeTOF(float min, float max)       { setCutPtotRangeTOF(min, max, StPicoCutsBase::kProton); }
inline void StPicoCutsBase::setCutProtonPtotRangeHybridTOF(float min, float max) { setCutPtotRangeHybridTOF(min, max, StPicoCutsBase::kProton); }

inline void StPicoCutsBase::setHybridTof(bool t) {mHybridTof = t;}
inline void StPicoCutsBase::setHybridTofKaon(bool t) {mHybridTofKaon = t;}
inline void StPicoCutsBase::setHybridTofPion(bool t) {mHybridTofPion = t;}
inline void StPicoCutsBase::setHybridTofBetterBetaCuts(bool t) {mHybridTofBetterBetaCuts = t;}
inline void StPicoCutsBase::setHybridTofBetterBetaCutsKaon(bool t) {mHybridTofBetterBetaCutsKaon = t;}
inline void StPicoCutsBase::setHybridTofBetterBetaCutsPion(bool t) {mHybridTofBetterBetaCutsPion = t;}

inline void StPicoCutsBase::setHybridTofWithBEMC(bool t) {mHybridTofWithBEMC = t;}



inline const float&    StPicoCutsBase::getHypotheticalMass(int pidFlag)        const { return mHypotheticalMass[pidFlag]; }

// -- check for good hadrons in TPC - in ptRange
inline bool StPicoCutsBase::isTPCPion(StPicoTrack const * const trk)   const {return isTPCHadron(trk, StPicoCutsBase::kPion); }
inline bool StPicoCutsBase::isTPCKaon(StPicoTrack const * const trk)   const {return isTPCHadron(trk, StPicoCutsBase::kKaon); }
inline bool StPicoCutsBase::isTPCProton(StPicoTrack const * const trk) const {return isTPCHadron(trk, StPicoCutsBase::kProton); }

inline bool StPicoCutsBase::isTOFPion(StPicoTrack const *trk)   const { float tofBeta = getTofBeta(trk);
    return isTOFHadron(trk, tofBeta, StPicoCutsBase::kPion); }
inline bool StPicoCutsBase::isTOFKaon(StPicoTrack const *trk)   const { float tofBeta = getTofBeta(trk);
    return isTOFHadron(trk, tofBeta, StPicoCutsBase::kKaon); }
inline bool StPicoCutsBase::isTOFProton(StPicoTrack const *trk) const { float tofBeta = getTofBeta(trk);
    return isTOFHadron(trk, tofBeta, StPicoCutsBase::kProton); }

inline bool StPicoCutsBase::isTOFPion(StPicoTrack const *trk,   float const & tofBeta) const { return isTOFHadron(trk, tofBeta, StPicoCutsBase::kPion); }
inline bool StPicoCutsBase::isTOFKaon(StPicoTrack const *trk,   float const & tofBeta) const { return isTOFHadron(trk, tofBeta, StPicoCutsBase::kKaon); }
inline bool StPicoCutsBase::isTOFProton(StPicoTrack const *trk, float const & tofBeta) const { return isTOFHadron(trk, tofBeta, StPicoCutsBase::kProton); }

inline bool StPicoCutsBase::isHybridTOFPion(StPicoTrack const *trk)   const { float tofBeta = getTofBeta(trk);
    return isHybridTOFHadron(trk, tofBeta, StPicoCutsBase::kPion); }
inline bool StPicoCutsBase::isHybridTOFKaon(StPicoTrack const *trk)   const { float tofBeta = getTofBeta(trk);
    return isHybridTOFHadron(trk, tofBeta, StPicoCutsBase::kKaon); }
inline bool StPicoCutsBase::isHybridTOFProton(StPicoTrack const *trk) const { float tofBeta = getTofBeta(trk);
    return isHybridTOFHadron(trk, tofBeta, StPicoCutsBase::kProton); }

inline bool StPicoCutsBase::isHybridTOFPion(StPicoTrack const *trk,   float const & tofBeta) const { return isHybridTOFHadron(trk, tofBeta, StPicoCutsBase::kPion); }
inline bool StPicoCutsBase::isHybridTOFKaon(StPicoTrack const *trk,   float const & tofBeta) const { return isHybridTOFHadron(trk, tofBeta, StPicoCutsBase::kKaon); }
inline bool StPicoCutsBase::isHybridTOFProton(StPicoTrack const *trk, float const & tofBeta) const { return isHybridTOFHadron(trk, tofBeta, StPicoCutsBase::kProton); }

inline bool StPicoCutsBase::isTOFKaonBetterCuts(StPicoTrack const *trk)   const { float tofBeta = getTofBeta(trk);
    return isTOFBetterKaon(trk, tofBeta, StPicoCutsBase::kKaon); }
inline bool StPicoCutsBase::isTOFKaonBetterCuts(StPicoTrack const *trk,   float const & tofBeta) const { return isTOFBetterKaon(trk, tofBeta, StPicoCutsBase::kKaon); }
inline bool StPicoCutsBase::isHybridTOFKaonBetterCuts(StPicoTrack const *trk)   const { float tofBeta = getTofBeta(trk);
    return isHybridTOFBetterKaon(trk, tofBeta, StPicoCutsBase::kKaon); }
inline bool StPicoCutsBase::isHybridTOFKaonBetterCuts(StPicoTrack const *trk,   float const & tofBeta) const { return isHybridTOFBetterKaon(trk, tofBeta, StPicoCutsBase::kKaon); }

inline bool StPicoCutsBase::isTOFPionBetterCuts(StPicoTrack const *trk)   const { float tofBeta = getTofBeta(trk);
    return isTOFBetterPion(trk, tofBeta, StPicoCutsBase::kPion); }
inline bool StPicoCutsBase::isTOFPionBetterCuts(StPicoTrack const *trk,   float const & tofBeta) const { return isTOFBetterPion(trk, tofBeta, StPicoCutsBase::kPion); }
inline bool StPicoCutsBase::isHybridTOFPionBetterCuts(StPicoTrack const *trk)   const { float tofBeta = getTofBeta(trk);
    return isHybridTOFBetterPion(trk, tofBeta, StPicoCutsBase::kPion); }
inline bool StPicoCutsBase::isHybridTOFPionBetterCuts(StPicoTrack const *trk,   float const & tofBeta) const { return isHybridTOFBetterPion(trk, tofBeta, StPicoCutsBase::kPion); }

/*inline void StPicoCutsBase::addTriggerByName(std::string name) {
    if (name == "y15ppall") {
        mVecTriggerIdList.insert(mVecTriggerIdList.end(), {470202, 480202, 490202, 470404, 480404, 490404, 470401, 480401, 490401, 470405, 480405, 490405, 470402, 480402, 490402});
    } else if (name == "y15ppht1") {
        mVecTriggerIdList.insert(mVecTriggerIdList.end(), {470202, 480202, 490202, 500202, 510202}); //BHT1*VPDMB-30
    } else if (name == "y15ppjp1") {
        mVecTriggerIdList.insert(mVecTriggerIdList.end(), {470404, 480404, 480414, 490404}); //JP1
    } else if (name == "y15ppjp2") {
        mVecTriggerIdList.insert(mVecTriggerIdList.end(), {470401, 480401, 480411, 490401, 500401, 500411}); //JP2
    } else if (name == "y15ppjp2l2") {
        mVecTriggerIdList.insert(mVecTriggerIdList.end(), {470405, 480405, 480415, 490405, 500405, 500415}); //JP2*L2JetHigh
    } else if (name == "y15ppjp2bsmd") {
        mVecTriggerIdList.insert(mVecTriggerIdList.end(), {470402, 480402, 490402, 500402, 500412}); //JP2-bsmd
    } else if (name == "y15ppmb") {
        mVecTriggerIdList.insert(mVecTriggerIdList.end(), {470011, 470021, 480001, 490001, 500001}); //VPDMB-5-ssd; st_ssdmb
    }
}*/

inline void StPicoCutsBase::showTriggers() {
    std::cout << "Selected trigger ids: " << std::endl;
    for(std::vector<unsigned int>::const_iterator iter = mVecTriggerIdList.begin(); iter != mVecTriggerIdList.end(); ++iter) {
        std::cout << *iter << " ";
    }
    std::cout << std::endl;
}



#endif
