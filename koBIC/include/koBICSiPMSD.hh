#ifndef koBICSiPMSD_h
#define koBICSiPMSD_h 1

#include "koBICSiPMHit.hh"
#include "koBICInterface.h"

#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"

class koBICSiPMSD : public G4VSensitiveDetector {
public:
  koBICSiPMSD(const G4String& name, const G4String& hitsCollectionName, const G4int& isLeft, koBICInterface::koBICModuleProperty ModuleProp, G4int modNum);
  virtual ~koBICSiPMSD();

  virtual void Initialize(G4HCofThisEvent* HCE);
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent* HCE);

private:
  koBICSiPMHitsCollection* fHitCollection;
  G4int fHCID;
  G4int fWavBin;
  G4int fTimeBin;
  G4float fWavlenStart;
  G4float fWavlenEnd;
  G4float fTimeStart;
  G4float fTimeEnd;
  G4float fWavlenStep;
  G4float fTimeStep;

  G4int fModuleNum;
  G4int fisLeft;
  koBICInterface::hitXY fTowerXY;

  G4double wavToE(G4double wav) { return h_Planck*c_light/wav; }

  koBICInterface::hitRange findWavRange(G4double en);
  koBICInterface::hitRange findTimeRange(G4double stepTime);
  koBICInterface::hitXY findSiPMXY(G4int SiPMnum, koBICInterface::hitXY towerXY);
};

#endif
