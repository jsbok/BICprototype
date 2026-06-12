#ifndef koBICRunAction_h
#define koBICRunAction_h 1

#include "RootInterface.h"
#include "koBICInterface.h"

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class koBICRunAction : public G4UserRunAction {
public:
  koBICRunAction(G4int seed, G4String filename, G4bool useHepMC);
  virtual ~koBICRunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

  static RootInterface<koBICInterface::koBICEventData>* sRootIO;
  static int sNumEvt;

private:
  G4int fSeed;
  G4String fFilename;
  G4bool fUseHepMC;
};

#endif
