#ifndef koBICActionInitialization_h
#define koBICActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "globals.hh"

class G4GenericMessenger;

class koBICActionInitialization : public G4VUserActionInitialization {
public:
  koBICActionInitialization(G4int seed, G4String filename);
  virtual ~koBICActionInitialization();

  virtual void BuildForMaster() const;
  virtual void Build() const;

private:
  void DefineCommands();

  G4GenericMessenger* fMessenger;
  G4int fSeed;
  G4String fFilename;
  G4bool fUseHepMC;
  G4bool fUseCalib;
  G4bool fUseGPS;
};

#endif
