#include <limits>

#include "TNtuple.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TList.h"

#include "StPicoMixedEventMaker/StPicoEventMixer.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoDst.h"

#include "StPicoMixedEventMaker/StPicoMixedEventMaker.h"
#include "StPicoMixedEventMaker/StMixerEvent.h"
#include "StPicoMixedEventMaker/StMixerHists.h"
#include "StPicoHFMaker/StHFCuts.h"
#include "StPicoCutsBase/StPicoCutsBase.h"
#include "StPicoMixedEventMaker/StMixerPair.h"


using namespace std;

ClassImp(StPicoEventMixer)


//-----------------------------------------------------------
StPicoEventMixer::StPicoEventMixer(char* category):
//        mEvents(),
//        mHists(NULL),
//        mHFCuts(NULL),
//        mEventsBuffer(5),
//        filledBuffer(0),
        mSETupleSig(NULL),
        mSETupleBack(NULL),
        mMETupleSig(NULL),
        mMETupleBack(NULL)
{
//    mHists = new StMixerHists(category);
}
//-----------------------------------------------------------
StPicoEventMixer::~StPicoEventMixer()
{
    for(unsigned int i =0 ; i<mEvents.size() ; i++){
        delete mEvents.at(i);
    }
}
//-----------------------------------------------------------
void StPicoEventMixer::finish() {
//    mHists->closeFile();
}

//StPicoCutsBase* picoCutsBase = new StPicoCutsBase("picoCutsBase");


//-----------------------------------------------------------
bool StPicoEventMixer::addPicoEvent(StPicoDst const* const picoDst, float weight) {

    if (!mHFCuts){
        cout<< "Nejsou mHFCuts (addPico)"<<endl;
        }
    unsigned int nTracks = picoDst->numberOfTracks();
    TVector3 pVertex = picoDst->event()->primaryVertex();
    StMixerEvent* event = new StMixerEvent(pVertex, picoDst->event()->bField());
    event->addPicoEvent(*(picoDst->event()));


    cout << "Testík 5" <<endl;

    for(unsigned int iTrk = 0; iTrk < nTracks; ++iTrk) {
        StPicoTrack const* trk = picoDst->track(iTrk);
        bool saveTrack = false;

        mEventsBuffer=10;
        filledBuffer=12;
        cout << "Tady jsem"<<endl;
        cout << "Event Buffer je " << mEventsBuffer <<endl;
        cout << "Filled buffer je " << filledBuffer <<endl;

        cout << "Testík 6" <<endl;

        if(mHFCuts->isGoodPion(trk)) {
            cout << "Testík 6.1" <<endl;

            saveTrack = true;
            event->addPion(event->getNoTracks());

        }
        cout << "Testík 7" <<endl;

        if(mHFCuts->isGoodKaon(trk)) {
            saveTrack = true;
            event->addKaon(event->getNoTracks());
//            cout<<mHFCuts->getOneOverBeta(trk, mHFCuts->getTofBetaBase(trk), StPicoCutsBase::kKaon)<<endl;
        }
        cout << "Testík 8" <<endl;

        if(saveTrack){
            cout << "Test 9"<<endl;
            event->addTrack(*trk);
            cout << "Test v posledním ifu"<<endl;
        } cout << "Test za posledním ifem" << endl;
    }
    cout <<"Poslední cout"<<endl;

    cout <<"Pocet pionu: "<<event->getNoPions()<<" Pocet kaonu: "<< event->getNoKaons()<<endl;

    if ( event->getNoPions() > 0 ||  event->getNoKaons() > 0) {
        mEvents.push_back(event);
        filledBuffer+=1;
        cout << "Předposlední if: Filled buffer je " << filledBuffer << " a Event Buffer" << mEventsBuffer <<endl;

    }




    //Returns true if need to do mixing, false if buffer has space still
    if ( filledBuffer == mEventsBuffer) {
        cout << "Poslední if: Filled buffer je " << filledBuffer << " a Event Buffer" << mEventsBuffer <<endl;
        cout << "True" <<endl;
        return true;
    }
    return false;
}
//-----------------------------------------------------------
void StPicoEventMixer::mixEvents() {
    cout << "Test mixEvents" <<endl;


    size_t const nEvent = mEvents.size();
    cout << "Pocet tracku" <<nEvent<<endl;

    int const nTracksEvt1 = mEvents.at(0)->getNoPions();
    cout << "Pocet pionu" <<nTracksEvt1<<endl;


    // Go through the event buffer
    for( size_t iEvt2 = 0; iEvt2 < nEvent; ++iEvt2) {
        int const nTracksEvt2 = mEvents.at(iEvt2)->getNoKaons();
        cout << "Pocet kaonu" <<nTracksEvt2<<endl;

        // evts trk loops
        for(int iTrk1 = 0; iTrk1 < nTracksEvt1; ++iTrk1) { //pions
            for( int iTrk2 = 0; iTrk2 < nTracksEvt2; ++iTrk2) { //kaons
                // check if the tracks are the same
                if(iEvt2 == 0) {
                    if(mEvents.at(0)->pionId(iTrk1) == mEvents.at(iEvt2)->kaonId(iTrk2)){
//                        cout<<"TRACKS ARE THE SAME"<<endl;
                        continue;
                    }
                }

                StPicoTrack const pion = mEvents.at(0)->pionAt(iTrk1);
                StPicoTrack const kaon = mEvents.at(iEvt2)->kaonAt(iTrk2);

                StMixerPair *pair = new StMixerPair(pion, kaon,
                                 mHFCuts->getHypotheticalMass(StHFCuts::kPion),
                                 mHFCuts->getHypotheticalMass(StHFCuts::kKaon),
                                 mEvents.at(0)->vertex(), mEvents.at(iEvt2)->vertex(),
                                 mEvents.at(0)->field() );

                if (!isCloseMixerPair(pair)) continue;

                int ii=0;
                float ntVar[21];
                ntVar[ii++] = pion.gPt();
                ntVar[ii++] = pair->particle1Dca();
                ntVar[ii++] = pion.nSigmaPion();
                ntVar[ii++] = pion.nHitsFit();
                ntVar[ii++] = mHFCuts->getOneOverBeta(&pion, mHFCuts->getTofBetaBase(&pion), StPicoCutsBase::kPion); //probably not working now, you are in weong event, so calc. for TOF traits is not right
                ntVar[ii++] = kaon.gPt();
                ntVar[ii++] = pair->particle2Dca();
                ntVar[ii++] = kaon.nSigmaKaon();
                ntVar[ii++] = kaon.nHitsFit();
                ntVar[ii++] = mHFCuts->getOneOverBeta(&kaon, mHFCuts->getTofBetaBase(&kaon), StPicoCutsBase::kKaon);
                ntVar[ii++] = pair->dcaDaughters();
                ntVar[ii++] = pair->rapidity();
                ntVar[ii++] = pair->pointingAngle();
                ntVar[ii++] = cos(pair->pointingAngle());
                ntVar[ii++] = pair->decayLength();
                ntVar[ii++] = pair->DcaToPrimaryVertex(); //(pair->decayLength())*sin(pair->pointingAngle());
                ntVar[ii++] = pair->phi();
                ntVar[ii++] = pair->eta();
                ntVar[ii++] = pair->cosThetaStar();
                ntVar[ii++] = pair->pt();
                ntVar[ii++] = pair->m();



                int charge = mEvents.at(0)->pionAt(iTrk1).charge() +  mEvents.at(iEvt2)->kaonAt(iTrk2).charge(); // 0 = signal

                if(iEvt2 == 0)
                    fillNtpSameEvtPair(ntVar, charge );
                else
                    fillNtpMixedEvtPair(ntVar, charge);
            }
        } //second track track loop
    } // the first track track loop

--filledBuffer;
delete mEvents.at(0);
mEvents.erase(mEvents.begin());
}

void StPicoEventMixer::fillNtpSameEvtPair(float ntVar[], int charge) {
    if(charge == 0 )
        mSETupleSig -> Fill(ntVar);
    else
        mSETupleBack-> Fill(ntVar);
    return;
}


void StPicoEventMixer::fillNtpMixedEvtPair(float ntVar[], int charge) {
    if(charge == 0 )
        mMETupleSig -> Fill(ntVar);
    else
        mMETupleBack-> Fill(ntVar);
    return;
}

//// _________________________________________________________
//void StPicoEventMixer::fillTracks(StMixerEvent* evt, bool isSameEvt, int pidFlag)
//{
////    if (!fillSinglePartHists)
////        return;
//
//    // get the corresponting histograms and track vectors
//    const std::string evtName = isSameEvt ? "SE" : "ME";
//    std::string particleName;
//    int nTracks;
//    switch (pidFlag)
//    {
//        case StHFCuts::kProton:
//            particleName = "p";
//            nTracks = evt->getNoProtons();
//            break;
//        case StHFCuts::kKaon:
//            particleName = "K";
//            nTracks = evt->getNoKaons();
//            break;
//        case StHFCuts::kPion:
//            particleName = "pi";
//            nTracks = evt->getNoPions();
//            break;
//        default:
//            cerr << "StPicoEventMixer::fillTracks: unknown pidFlag ... exiting" << endl;
//            throw;
//    }
////    TH2D *etaPhiHist = static_cast<TH2D*>(mSingleParticleList->FindObject(Form("%sEtaPhi%s",particleName.data(), evtName.data())));
////    TH2D *phiPtHist  = static_cast<TH2D*>(mSingleParticleList->FindObject(Form("%sPhiPt%s",particleName.data(), evtName.data())));
////    TH1D *dcaHist = static_cast<TH1D*>(mSingleParticleList->FindObject(Form("%sDCA%s",particleName.data(), evtName.data())));
////    TH1D *nTracksHist = static_cast<TH1D*>(mSingleParticleList->FindObject(Form("%stracks%s",particleName.data(), evtName.data())));
////    nTracksHist->Fill(nTracks);
//
//    // particle loop
////    const float weight = evt->weight();
//    for (int i = 0; i < nTracks; ++i)
//    {
//        StPicoTrack trk;
////        StMixerTrack trk;
//        switch (pidFlag)
//        {
//            case StHFCuts::kProton: trk = evt->protonAt(i);
//                break;
//            case StHFCuts::kKaon:   trk = evt->kaonAt(i);
//                break;
//            case StHFCuts::kPion:   trk = evt->pionAt(i);
//                break;
//        }
//        const float eta = trk.gMom().pseudoRapidity();
//        const float phi = trk.gMom().phi();
//        const float pt = trk.gMom().perp();
//        const float dca  = (trk.origin() - evt->vertex()).Mag();
//
////        etaPhiHist->Fill(phi,eta,weight);
////        phiPtHist->Fill(pt,phi, weight);
////        dcaHist->Fill(dca, weight);
//    }
//
//}
// _________________________________________________________
bool StPicoEventMixer::isMixerPion(StMixerTrack const& track) {
    short info = track.getTrackInfo();
    //TPC pion
    if( (info & (1 << kPionTPCbit)) >> kPionTPCbit != 1) return false;
    //TOF pion
    if( (info & (1 << kPionTOFbit)) >> kPionTOFbit != 1) return false;
    return true;
}
// _________________________________________________________
bool StPicoEventMixer::isMixerKaon(StMixerTrack const& track) {
    short info = track.getTrackInfo();
    //TPC Kaon
    if( (info & (1 << kKaonTPCbit)) >> kKaonTPCbit != 1) return false;
    //TOF Kaon
    if( (info & (1 << kKaonTOFbit)) >> kKaonTOFbit != 1) return false;
    return true;
}
//-----------------------------------------------------------
bool StPicoEventMixer::isMixerProton(StMixerTrack const& track) {
    short info = track.getTrackInfo();
    //TPC Proton
    if( (info & (1 << kProtonTPCbit)) >> kProtonTPCbit != 1) return false;
    //TOF Proton
    if( (info & (1 << kProtonTOFbit)) >> kProtonTOFbit != 1) return false;
    return true;
}
//-----------------------------------------------------------
/*bool StPicoEventMixer::isGoodEvent(StPicoDst const * const picoDst)
{
    return (mHFCuts->isGoodEvent(picoDst));
}
//-----------------------------------------------------------
bool StPicoEventMixer::isTpcPion(StPicoTrack const * const trk)
{
    return( isTPCHadron(trk, StPicoCutsBase::kPion));
}
//-----------------------------------------------------------
bool StPicoEventMixer::isTpcKaon(StPicoTrack const * const trk)
{
    return( isTPCHadron(trk, StPicoCutsBase::kKaon));
}
//-----------------------------------------------------------
bool StPicoEventMixer::isTpcProton(StPicoTrack const * const trk)
{
    return( isTPCHadron(trk, StPicoCutsBase::kProton));
}
//-----------------------------------------------------------
bool StPicoEventMixer::isTPCHadron(StPicoTrack const * const trk, int pidFlag)
{
    return( mHFCuts->isTPCHadron(trk, pidFlag));
}
//-----------------------------------------------------------
bool StPicoEventMixer::isGoodTrack(StPicoTrack const * const trk)
{
    return (mHFCuts->isGoodTrack(trk));
}
//-----------------------------------------------------------
bool StPicoEventMixer::isGoodTrigger(StPicoEvent const * const mPicoEvent) const
{
    return mHFCuts->isGoodTrigger(mPicoEvent);
}
 */
//-----------------------------------------------------------
bool StPicoEventMixer::isCloseMixerPair(StMixerPair const * pair) const {
    return ( std::cos(pair->pointingAngle()) > mHFCuts->cutSecondaryPairCosThetaMin() &&
             pair->decayLength() > mHFCuts->cutSecondaryPairDecayLengthMin() && pair->decayLength() < mHFCuts->cutSecondaryPairDecayLengthMax() &&
             pair->dcaDaughters() < mHFCuts->cutSecondaryPairDcaDaughtersMax() &&
             pair->m() > mHFCuts->cutSecondaryPairMassMin() && pair->m() < mHFCuts->cutSecondaryPairMassMax());
}


