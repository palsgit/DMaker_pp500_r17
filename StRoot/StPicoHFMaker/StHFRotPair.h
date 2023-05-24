#ifndef StHFRotPair_hh
#define StHFRotPair_hh

/* **************************************************
 *  Generic class calculating and storing pairs in HF analysis
 *  Allows to combine:
 *  - two particles, using
 *      StHFRotPair(StPicoTrack const * particle1, StPicoTrack const * particle2, ...
 *  - a particle and another pair, using
 *      StHFRotPair(StPicoTrack const * particle1, StHFRotPair * particle2, ...
 *    - in the current implementation the incoming pair is seen as having charge = 0
 *    - after determining the vertex of particle and incoming pair, the 
 *      decay vertex (tertiary vertex) of incoming particle can be updated
 *    - straight line approximation is the default, but full helix can be used
 *
 * **************************************************
 *
 *  Initial Authors: 
 *            Xin Dong        (xdong@lbl.gov)
 *            Mustafa Mustafa (mmustafa@lbl.gov)
 *            Jochen Thaeder  (jmthader@lbl.gov)
 *          **Michael Lomnitz (mrlomnitz@lbl.gov)
 *
 *  ** Code Maintainer 
 *
 * **************************************************
 */

#include "TObject.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "StPicoEvent/StPicoEvent.h"

class StPicoTrack;

class StHFRotPair : public TObject
{
public:
	StHFRotPair();
	StHFRotPair(StHFRotPair const *);

	StHFRotPair(StPicoTrack const * particle1, StPicoTrack const * particle2,
			 float p1MassHypo, float p2MassHypo,
			 unsigned short p1Idx, unsigned short p2Idx,
			 TVector3 const & vtx, float bField, bool useStraightLine = true);

	~StHFRotPair();

	TLorentzVector const & lorentzVector() const;
	TVector3 const & decayVertex() const;
	float rapidity()    const;
	float m()    const;
	float pt()   const;
	float eta()  const;
	float phi()  const;
	float pointingAngle() const;
	float pointingAngle(TVector3 const & vtx2) const;
	float decayLength() const;
	float decayLength(TVector3 const & vtx2) const;
	float particle1Dca() const;
	float particle1Dca(StPicoTrack const * p1track, TVector3 const & vtx2, float bField) const;
	float particle2Dca() const;
	float particle2Dca(StPicoTrack const * p1track, TVector3 const & vtx2, float bField) const;
	unsigned short particle1Idx() const;
	unsigned short particle2Idx() const;
	float dcaDaughters() const;
	float cosThetaStar() const;
	float v0x() const;
	float v0y() const;
	float v0z() const;
	float px() const;
	float py() const;
	float pz() const;
	float DcaToPrimaryVertex() const;

private:
	StHFRotPair(StHFRotPair const &);
	StHFRotPair& operator=(StHFRotPair const &);
	TLorentzVector mLorentzVector;
	TVector3   mDecayVertex;

	float mPointingAngle;
	float mDecayLength;
	float mParticle1Dca;
	float mParticle2Dca;

	unsigned short  mParticle1Idx; // index of track in StPicoDstEvent
	unsigned short  mParticle2Idx; // index of track in StPicoDstEvent for particle, idx in tertiary vertex array for pair

	float mDcaDaughters;
	float mCosThetaStar;

	ClassDef(StHFRotPair,1)
};
inline TLorentzVector const & StHFRotPair::lorentzVector() const { return mLorentzVector;}
inline float StHFRotPair::rapidity()    const { return mLorentzVector.Rapidity();}
inline float StHFRotPair::m()    const { return mLorentzVector.M();}
inline float StHFRotPair::pt()   const { return mLorentzVector.Pt();}
inline float StHFRotPair::eta()  const { return mLorentzVector.Eta();}
inline float StHFRotPair::phi()  const { return mLorentzVector.Phi();}
inline float StHFRotPair::px()   const { return mLorentzVector.Px();}
inline float StHFRotPair::py()   const { return mLorentzVector.Py();}
inline float StHFRotPair::pz()   const { return mLorentzVector.Pz();}
inline float StHFRotPair::pointingAngle() const { return mPointingAngle;}
inline float StHFRotPair::decayLength()   const { return mDecayLength;}
inline float StHFRotPair::particle1Dca()  const { return mParticle1Dca;}
inline float StHFRotPair::particle2Dca()  const { return mParticle2Dca;}
inline unsigned short StHFRotPair::particle1Idx() const { return mParticle1Idx;}
inline unsigned short StHFRotPair::particle2Idx() const { return mParticle2Idx;}
inline float StHFRotPair::dcaDaughters() const { return mDcaDaughters;}
inline float StHFRotPair::cosThetaStar() const { return mCosThetaStar;}
inline TVector3 const & StHFRotPair::decayVertex() const { return mDecayVertex;}
inline float StHFRotPair::v0x() const { return mDecayVertex.x();}
inline float StHFRotPair::v0y() const { return mDecayVertex.y();}
inline float StHFRotPair::v0z() const { return mDecayVertex.z();}
inline float StHFRotPair::DcaToPrimaryVertex() const { return mDecayLength*std::sin(mPointingAngle); }
#endif

