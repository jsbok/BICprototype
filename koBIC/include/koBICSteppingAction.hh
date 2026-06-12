#ifndef koBICSteppingAction_h
#define koBICSteppingAction_h 1

#include "koBICInterface.h"
#include "koBICEventAction.hh"

#include "G4UserSteppingAction.hh"
#include "G4LogicalVolume.hh"
#include "G4Step.hh"

using namespace std;

class koBICSteppingAction : public G4UserSteppingAction {
public:
  koBICSteppingAction(koBICEventAction* eventAction);
  virtual ~koBICSteppingAction();
  virtual void UserSteppingAction(const G4Step*);

private:
  koBICEventAction* fEventAction;
  koBICInterface::koBICEdepData fEdep;
  koBICInterface::koBICLeakageData fLeak;
  
  G4VPhysicalVolume* GetMotherTower(G4TouchableHandle touchable) { return touchable->GetVolume(touchable->GetHistoryDepth()-1); }

  G4int GetModuleNum(G4String towerName) {
    return 0;
  }
};

#endif
