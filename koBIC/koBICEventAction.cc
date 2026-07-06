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
  G4RunManager::GetRunManager()->SetPrintProgress(1);
}

koBICEventAction::~koBICEventAction() {}

void koBICEventAction::BeginOfEventAction(const G4Event*) {
  // 1. 매 이벤트 시작 시 전역/멤버 맵 및 벡터 클리어 보장
  clear();

  // 2. [원래 구조 복원] 원래 작동하던 방식대로 0번부터 27번까지 0점짜리 기본 방(Data)을 미리 생성
  // 이렇게 해두면 중성미자를 쏴서 물리 반응이 전혀 없어도 항상 28개의 빈 데이터가 보장됩니다.
  for (int i = 0; i < 28; i++) {
    // 에너지(Edep) 맵 초기화
    koBICInterface::koBICEdepData defaultEdep;
    defaultEdep.ModuleNum = i;
    defaultEdep.EdepCore = 0.0;
    defaultEdep.Edep = 0.0;
    defaultEdep.EdepEle = 0.0;
    defaultEdep.EdepGamma = 0.0;
    defaultEdep.EdepCharged = 0.0;
    fEdepMap[i] = defaultEdep;

    // 타워(SiPM Hit) 맵 초기화
    koBICInterface::koBICTowerData defaultTower;
    defaultTower.ModuleNum = i;
    defaultTower.numx = 0; // 기본값 0, fillHits가 호출되면 업데이트됨
    defaultTower.numy = 0;
    // defaultTower.SiPMs 벡터는 비어있는 상태(0개)로 둠
    fTowerMap[i] = defaultTower;
  }

  // 3. Collection ID는 런타임 중 최초 1회만 안전하게 로드하도록 변경 (무한 push_back 차단)
  if (fSiPMCollID.empty()) {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    for (int i = 0; i < 28; i++) {
      G4int id1 = sdManager->GetCollectionID("ModuleC" + std::to_string(2 * i));
      G4int id2 = sdManager->GetCollectionID("ModuleC" + std::to_string(2 * i + 1));
      
      if (id1 >= 0) fSiPMCollID.push_back(id1);
      if (id2 >= 0) fSiPMCollID.push_back(id2);
    }
  }

  fEventData = new koBICInterface::koBICEventData();
}

void koBICEventAction::clear() {
  // fSiPMCollID는 static 성격으로 유지하므로 여기서는 clear하지 않습니다.
  fTowerMap.clear();
  fEdepMap.clear();
}

void koBICEventAction::EndOfEventAction(const G4Event* event) {
  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  if (!hce) {
    G4ExceptionDescription msg;
    msg << "No hits collection of this event found." << G4endl;
    msg << "koBICEventAction::EndOfEventAction()", "koBICCode001", JustWarning, msg;
    return;
  }

  // 데이터 추출 및 fillHits 처리 (반응이 있는 모듈만 맵의 데이터를 업데이트함)
  for (int iSD = 0; iSD < (int)fSiPMCollID.size(); iSD++) {
    G4int hcID = fSiPMCollID[iSD];
    if (hcID < 0) continue; 

    koBICSiPMHitsCollection* sipmHC = (koBICSiPMHitsCollection*)(hce->GetHC(hcID));

    if (sipmHC) {
      G4int SiPMs = sipmHC->entries();
      for (G4int iHC = 0; iHC < SiPMs; iHC++) {
        // 유효성 체크 방어 코드 추가
        if ((*sipmHC)[iHC]) {
          fillHits((*sipmHC)[iHC]);
        }
      }
    }
  }

  // 🔥 [핵심 수정] 0번부터 27번까지 빠짐없이 "순서대로" EventData 구조체 벡터에 꽉 채워 담기
  // 반응이 없던 모듈(중성미자 포함)은 Begin에서 채워둔 0점짜리 기본 데이터가 그대로 들어갑니다.
  for (int i = 0; i < 28; i++) {
    fEventData->towers.push_back(fTowerMap[i]);
    fEventData->Edeps.push_back(fEdepMap[i]);
  }

  // Primary Particle 정보 저장
  for (int iVtx = 0; iVtx < event->GetNumberOfPrimaryVertex(); iVtx++) {
    G4PrimaryVertex* vtx = event->GetPrimaryVertex(iVtx);
    if (!vtx) continue;

    for (int iPtc = 0; iPtc < vtx->GetNumberOfParticle(); iPtc++) {
      G4PrimaryParticle* ptc = vtx->GetPrimary(iPtc);
      if (ptc) {
        fillPtcs(vtx, ptc);
      }
    }
  }

  fEventData->event_number = koBICPrimaryGeneratorAction::sIdxEvt;

  // 스레드 대기 큐 진입 및 ROOT 파일 전송
  queue();

  // 사용이 끝난 메모리 동적 해제
  delete fEventData;
  fEventData = nullptr;
}

void koBICEventAction::fillHits(koBICSiPMHit* hit) {
  if (!hit) return;

  // ⚠️ 모듈 번호 범위 예외 처리 (0~27 범위를 벗어나면 에러 유발하므로 즉시 버림)
  G4int modNum = hit->GetModuleNum();
  if (modNum < 0 || modNum >= 28) return;

  koBICInterface::koBICSiPMData sipmData;
  sipmData.count = hit->GetPhotonCount();
  sipmData.SiPMnum = hit->GetSiPMnum();
  sipmData.x = hit->GetSiPMXY().first;
  sipmData.y = hit->GetSiPMXY().second;
  sipmData.isleft = hit->GetisLeft();
  sipmData.pos = std::make_tuple(hit->GetSiPMpos().x(), hit->GetSiPMpos().y(), hit->GetSiPMpos().z());
  sipmData.timeStruct = hit->GetTimeStruct();
  sipmData.wavlenSpectrum = hit->GetWavlenSpectrum();

  // BeginOfEventAction에서 이미 0~27번 방이 확실히 존재하므로 find 없이 즉시 접근하여 수정합니다.
  fTowerMap[modNum].numx = hit->GetTowerXY().first;
  fTowerMap[modNum].numy = hit->GetTowerXY().second;
  fTowerMap[modNum].SiPMs.push_back(sipmData);
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
  // ⚠️ 모듈 번호 범위 최종 확인 방어코드
  if (edepData.ModuleNum < 0 || edepData.ModuleNum >= 28) return;

  // BeginOfEventAction에서 이미 방이 만들어져 있으므로 대입(=) 대신 기존 값에 누적(+=)해 줍니다.
  fEdepMap[edepData.ModuleNum].EdepCore += edepData.EdepCore;
  fEdepMap[edepData.ModuleNum].Edep     += edepData.Edep;
  fEdepMap[edepData.ModuleNum].EdepEle  += edepData.EdepEle;
  fEdepMap[edepData.ModuleNum].EdepGamma += edepData.EdepGamma;
  fEdepMap[edepData.ModuleNum].EdepCharged += edepData.EdepCharged;
}

void koBICEventAction::fillLeaks(koBICInterface::koBICLeakageData leakData) {
  if (fEventData) {
    fEventData->leaks.push_back(leakData);
  }
}

void koBICEventAction::queue() {
  while (koBICRunAction::sNumEvt != koBICPrimaryGeneratorAction::sIdxEvt) {
    G4AutoLock lock(&koBICEventActionMutex);
    if (koBICRunAction::sNumEvt == koBICPrimaryGeneratorAction::sIdxEvt) break;
    G4CONDITIONWAIT(&koBICEventActionCV, &lock);
  }
  G4AutoLock lock(&koBICEventActionMutex);
  
  if (fEventData && koBICRunAction::sRootIO) {
    koBICRunAction::sRootIO->fill(fEventData);
  }
  
  koBICRunAction::sNumEvt++;
  G4CONDITIONBROADCAST(&koBICEventActionCV);
}
