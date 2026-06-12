#ifndef koBICEventAction_h
#define koBICEventAction_h 1

#include "koBICInterface.h"
#include "koBICSiPMHit.hh"

#include "G4UserEventAction.hh"
#include "G4HCofThisEvent.hh"
#include "G4Event.hh"

class koBICEventAction : public G4UserEventAction {
public:

  koBICEventAction();
  virtual ~koBICEventAction();

  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  void fillEdeps(koBICInterface::koBICEdepData edepData);
  void fillLeaks(koBICInterface::koBICLeakageData leakData);

private:
  void clear();
  void fillHits(koBICSiPMHit* hit);
  void fillPtcs(G4PrimaryVertex* vtx, G4PrimaryParticle* ptc);
  void queue();

  koBICInterface::koBICEventData* fEventData;
  std::map<int, koBICInterface::koBICTowerData> fTowerMap;
  std::map<int, koBICInterface::koBICEdepData> fEdepMap;

  std::vector<G4int> fSiPMCollID;
};

#endif
