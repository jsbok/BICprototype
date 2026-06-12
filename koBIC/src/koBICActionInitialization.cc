#include "koBICActionInitialization.hh"
#include "koBICPrimaryGeneratorAction.hh"
#include "koBICRunAction.hh"
#include "koBICEventAction.hh"
#include "koBICSteppingAction.hh"

#include "G4GenericMessenger.hh"

using namespace std;
koBICActionInitialization::koBICActionInitialization(G4int seed, G4String filename)
: G4VUserActionInitialization()
{
  fSeed = seed;
  fFilename = filename;

  DefineCommands();
}

koBICActionInitialization::~koBICActionInitialization() {
  if (fMessenger) delete fMessenger;
}

void koBICActionInitialization::BuildForMaster() const {
  SetUserAction(new koBICRunAction(fSeed,fFilename,fUseHepMC));
}

void koBICActionInitialization::Build() const {
  SetUserAction(new koBICPrimaryGeneratorAction(fSeed,fUseHepMC,fUseCalib,fUseGPS));
  SetUserAction(new koBICRunAction(fSeed,fFilename,fUseHepMC));

  koBICEventAction* eventAction = new koBICEventAction();
  SetUserAction(eventAction);

  SetUserAction(new koBICSteppingAction(eventAction));
}

void koBICActionInitialization::DefineCommands() {
  fMessenger = new G4GenericMessenger(this, "/koBIC/action/", "action initialization control");
  G4GenericMessenger::Command& ioCmd = fMessenger->DeclareProperty("useHepMC",fUseHepMC,"use HepMC");
  ioCmd.SetParameterName("useHepMC",true);
  ioCmd.SetDefaultValue("False");

  G4GenericMessenger::Command& calibCmd = fMessenger->DeclareProperty("useCalib",fUseCalib,"use Calib");
  calibCmd.SetParameterName("useCalib",true);
  calibCmd.SetDefaultValue("False");

  G4GenericMessenger::Command& gpsCmd = fMessenger->DeclareProperty("useGPS",fUseGPS,"use GPS");
  gpsCmd.SetParameterName("useGPS",true);
  gpsCmd.SetDefaultValue("False");
}
