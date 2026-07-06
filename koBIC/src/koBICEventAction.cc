#include "koBICEventAction.hh"
#include "koBICRunAction.hh"
#include "koBICPrimaryGeneratorAction.hh"
#include "koBICDetectorConstruction.hh"

#include "G4PrimaryVertex.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4AutoLock.hh"
#include "G4Threading.hh"

namespace {
  G4Mutex koBICEventActionMutex = G4MUTEX_INITIALIZER;
  G4Condition koBICEventActionCV = G4CONDITION_INITIALIZER;
}

koBICEventAction::koBICEventAction()
: G4UserEventAction()
{
  // set printing per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);
}

koBICEventAction::~koBICEventAction() {}

void koBICEventAction::BeginOfEventAction(const G4Event*) {
	clear();

  G4SDManager* sdManager = G4SDManager::GetSDMpointer();
  for (int i = 0; i < 28; i++) {
    fSiPMCollID.push_back(sdManager->GetCollectionID("ModuleC"+std::to_string(2*i)));
    fSiPMCollID.push_back(sdManager->GetCollectionID("ModuleC"+std::to_string(2*i+1)));
  }

  fEventData = new koBICInterface::koBICEventData();
}

void koBICEventAction::clear() {
  fSiPMCollID.clear();
  fTowerMap.clear();
  fEdepMap.clear();
}

void koBICEventAction::EndOfEventAction(const G4Event* event) {
  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  if (!hce) {
    G4ExceptionDescription msg;
    msg << "No hits collection of this event found." << G4endl;
    G4Exception("koBICEventAction::EndOfEventAction()",
    "koBICCode001", JustWarning, msg);
    return;
  }

  G4int totSDNum = hce->GetNumberOfCollections();

for (int iSD = 0; iSD < (int)fSiPMCollID.size(); iSD++) {
    G4int hcID = fSiPMCollID[iSD];
    if (hcID < 0) continue; 

    koBICSiPMHitsCollection* sipmHC = (koBICSiPMHitsCollection*)(hce->GetHC(hcID));

    if (sipmHC) {
      G4int SiPMs = sipmHC->entries();
      for (G4int iHC = 0; iHC < SiPMs; iHC++) {
        fillHits((*sipmHC)[iHC]);
      }
    }
  }

  for (const auto& towerMap : fTowerMap) {
    fEventData->towers.push_back(towerMap.second);
  }

  for (const auto& edepMap : fEdepMap) {
    fEventData->Edeps.push_back(edepMap.second);
  }

  for (int iVtx = 0; iVtx < event->GetNumberOfPrimaryVertex(); iVtx++) {
    G4PrimaryVertex* vtx = event->GetPrimaryVertex(iVtx);

    for (int iPtc = 0; iPtc < vtx->GetNumberOfParticle(); iPtc++) {
      G4PrimaryParticle* ptc = vtx->GetPrimary(iPtc);
      fillPtcs(vtx,ptc);
    }
  }

  fEventData->event_number = koBICPrimaryGeneratorAction::sIdxEvt;

  queue();

  delete fEventData;
}

void koBICEventAction::fillHits(koBICSiPMHit* hit) {
  koBICInterface::koBICSiPMData sipmData;
  sipmData.count = hit->GetPhotonCount();
  sipmData.SiPMnum = hit->GetSiPMnum();
  sipmData.x = hit->GetSiPMXY().first;
  sipmData.y = hit->GetSiPMXY().second;
  sipmData.isleft = hit->GetisLeft();
  sipmData.pos = std::make_tuple(hit->GetSiPMpos().x(),hit->GetSiPMpos().y(),hit->GetSiPMpos().z());
  sipmData.timeStruct = hit->GetTimeStruct();
  sipmData.wavlenSpectrum = hit->GetWavlenSpectrum();

  auto towerIter = fTowerMap.find(hit->GetModuleNum());

  if ( towerIter==fTowerMap.end() ) {
    koBICInterface::koBICTowerData towerData;
    towerData.ModuleNum = hit->GetModuleNum();
    towerData.numx = hit->GetTowerXY().first;
    towerData.numy = hit->GetTowerXY().second;
    // towerData.innerR = hit->GetTowerInnerR();
    // towerData.towerH = hit->GetTowerH();
    towerData.SiPMs.push_back(sipmData);

    fTowerMap.insert(std::make_pair(hit->GetModuleNum(),towerData));
  } else {
    towerIter->second.SiPMs.push_back(sipmData);
  }
}

void koBICEventAction::fillPtcs(G4PrimaryVertex* vtx, G4PrimaryParticle* ptc) {
  koBICInterface::koBICGenData GenData;
  GenData.E = ptc->GetTotalEnergy();
  GenData.px = ptc->GetPx();
  GenData.py = ptc->GetPy();
  GenData.pz = ptc->GetPz();
  GenData.pdgId = ptc->GetPDGcode();
  GenData.vx = vtx->GetX0();
  GenData.vy = vtx->GetY0();
  GenData.vz = vtx->GetZ0();
  GenData.vt = vtx->GetT0();

  fEventData->GenPtcs.push_back(GenData);
}

void koBICEventAction::fillEdeps(koBICInterface::koBICEdepData edepData) {
  auto towerIter = fEdepMap.find(edepData.ModuleNum);

  if ( towerIter==fEdepMap.end() ) {
    fEdepMap.insert(std::make_pair(edepData.ModuleNum,edepData));
  } else {
    towerIter->second.Edep += edepData.Edep;
    towerIter->second.EdepEle += edepData.EdepEle;
    towerIter->second.EdepGamma += edepData.EdepGamma;
    towerIter->second.EdepCharged += edepData.EdepCharged;
  }
}

void koBICEventAction::fillLeaks(koBICInterface::koBICLeakageData leakData) {
  fEventData->leaks.push_back(leakData);
}

void koBICEventAction::queue() {
  while ( koBICRunAction::sNumEvt != koBICPrimaryGeneratorAction::sIdxEvt ) {
    G4AutoLock lock(&koBICEventActionMutex);
    if ( koBICRunAction::sNumEvt == koBICPrimaryGeneratorAction::sIdxEvt ) break;
    G4CONDITIONWAIT(&koBICEventActionCV, &lock);
  }
  G4AutoLock lock(&koBICEventActionMutex);
  koBICRunAction::sRootIO->fill(fEventData);
  koBICRunAction::sNumEvt++;
  G4CONDITIONBROADCAST(&koBICEventActionCV);
}
