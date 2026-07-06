#include "koBICRunAction.hh"
#include "koBICEventAction.hh"
#include "G4AutoLock.hh"
#include "G4Threading.hh"
#include "G4Run.hh"
#include <vector>

using namespace std;

namespace { G4Mutex koBICRunActionMutex = G4MUTEX_INITIALIZER; }
RootInterface<koBICInterface::koBICEventData>* koBICRunAction::sRootIO = 0;
int koBICRunAction::sNumEvt = 0;

koBICRunAction::koBICRunAction(G4int seed, G4String filename, G4bool useHepMC)
: G4UserRunAction()
{
  fSeed = seed;
  fFilename = filename;
  fUseHepMC = useHepMC;

  G4AutoLock lock(&koBICRunActionMutex);

  if (!sRootIO) {
    sRootIO = new RootInterface<koBICInterface::koBICEventData>(filename+"_"+std::to_string(fSeed)+".root", true);
    sRootIO->create("koBIC","koBICEventData");
  }
}

koBICRunAction::~koBICRunAction() {
  G4AutoLock lock(&koBICRunActionMutex);
  if (IsMaster() && sRootIO) {
    sRootIO->close();
    delete sRootIO;
    sRootIO = nullptr;
  }
}

void koBICRunAction::BeginOfRunAction(const G4Run* aRun) {
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

void koBICRunAction::EndOfRunAction(const G4Run* aRun) {
  if (IsMaster()) {
    G4AutoLock lock(&koBICRunActionMutex);
    if (sRootIO) {
      G4cout << "--- Run " << aRun->GetRunID() << " ended. Writing data to ROOT file..." << G4endl;
      
      sRootIO->write(); 
      sRootIO->close(); 
      
      delete sRootIO;
      sRootIO = nullptr;
      G4cout << "--- ROOT file closed and deleted safely before geometry destruction." << G4endl;
    }
  }
}


