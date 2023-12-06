#ifndef StPicoD0AnaMaker_h
#define StPicoD0AnaMaker_h

#include "StPicoHFMaker/StPicoHFMaker.h"
#include "TNtuple.h"
#include "TVector3.h"

//#include "StRefMultCorr/StRefMultCorr.h"
#include "TH2F.h"
//#include "StPicoD0AnaHists.h"
#include <vector>
#include "TClonesArray.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoHFMaker/StPicoHFEvent.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoHFMaker/StHFPair.h"
#include "StPicoHFMaker/StHFTriplet.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "phys_constants.h"

#include "TH1F.h"
#include "TH3F.h"
#include "TProfile.h"
#include <ctime>

class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoHFEvent;

class StHFPair;
class StHFTriplet;
class StHFCuts;

class StPicoD0AnaMaker : public StPicoHFMaker
{
public:
    StPicoD0AnaMaker(char const*, StPicoDstMaker*, char const*);
    virtual ~StPicoD0AnaMaker();

    virtual Int_t InitHF();
    virtual Int_t MakeHF();
    virtual void  ClearHF(Option_t *opt);
    virtual Int_t FinishHF();
    void workWithRefit(bool);

protected:
    std::vector<unsigned short> mIdxPicoPions;
    std::vector<unsigned short> mIdxPicoKaons;
    std::vector<int> tracksToRemove;
    std::vector<int> primaryTracks;
    std::vector<std::vector<float>> vect_daughter_id;
    std::vector<float> daughter_id;

private:
    int nGoodTracks;

    int createCandidates();
    int analyzeCandidates();
    TVector3 refitVertex(bool);
    bool isVectorInVectorOfVectors(const std::vector<std::vector<float>> vectorOfVectors, const std::vector<float> targetVector);
    int ConvertRunID(int runID);


    TNtuple *ntp_DMeson_UnlikeSign;
    TNtuple *ntp_DMeson_Rotated;
    TNtuple *ntp_DMeson_LikeSign;

    TNtuple *ntp_DstarMeson_RightSign;
    TNtuple *ntp_DstarMeson_WrongSign;
    TNtuple *ntp_DstarMeson_SideBand;

    int RunId;
    std::vector<int> RunNumberVector;


    int mRunNumber;
    Int_t nPrimary;
    TString mOutFileBaseName;
    bool mSwitchRefit;

    TFile* mOutFile;

    TMVA::Reader *reader[3];
    Float_t k_dca[3], pi1_dca[3], dcaDaughters[3], cosTheta[3], D_decayL[3], dcaD0ToPv[3], thetaStar[3];

    ClassDef(StPicoD0AnaMaker, 1) //set to 1
};

inline void StPicoD0AnaMaker::workWithRefit(bool b) { mSwitchRefit = b; }

#endif
