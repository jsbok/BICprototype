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

  // 파일 생성: 파일명 뒤에 시드값이 붙어 생성됩니다 (예: filename_0.root)
  if (!sRootIO) {
    sRootIO = new RootInterface<koBICInterface::koBICEventData>(filename+"_"+std::to_string(fSeed)+".root", true);
    sRootIO->create("koBIC","koBICEventData");
  }
}

koBICRunAction::~koBICRunAction() {
  G4AutoLock lock(&koBICRunActionMutex);
  if (IsMaster() && sRootIO) {
    // EndOfRunAction에서 이미 write()를 했으므로 
    // 여기서는 close()만 확실히 하고 메모리만 비웁니다.
    sRootIO->close();
    // delete sRootIO; // 때에 따라 여기서 delete가 충돌을 일으키면 주석처리해도 무방합니다.
    sRootIO = 0;
  }
}
void koBICRunAction::BeginOfRunAction(const G4Run* aRun) {
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
}

void koBICRunAction::EndOfRunAction(const G4Run* aRun) {
  // Run이 끝날 때마다(beamOn이 끝날 때마다) 데이터를 파일에 씁니다.
  // 이렇게 해야 로컬에서 중도 종료해도 데이터가 보존됩니다.
  if (IsMaster()) {
    G4AutoLock lock(&koBICRunActionMutex);
    if (sRootIO) {
      G4cout << "--- Run " << aRun->GetRunID() << " ended. Writing data to ROOT file..." << G4endl;
      sRootIO->write(); 
      // 로컬 테스트 시 매번 닫고 싶다면 close()를 호출할 수 있으나, 
      // 통상적으로 write()만 해도 TTree 버퍼가 비워지며 파일에 기록됩니다.
    }
  }
}
