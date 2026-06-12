#include "koBICSteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VProcess.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpticalPhoton.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VSensitiveDetector.hh"

koBICSteppingAction::koBICSteppingAction(koBICEventAction* eventAction)
: G4UserSteppingAction(), fEventAction(eventAction)
{}

koBICSteppingAction::~koBICSteppingAction() {}

void koBICSteppingAction::UserSteppingAction(const G4Step* step) {
  G4Track* track = step->GetTrack();

static long long nCreated = 0;
static long long lastPrint = 0;

    const auto* secondaries =
        step->GetSecondaryInCurrentStep();

static double totEdep = 0;
static long long totPhoton = 0;

totEdep += step->GetTotalEnergyDeposit();

for(const auto* sec : *secondaries)
{
    if(sec->GetDefinition()==G4OpticalPhoton::Definition())
        totPhoton++;
        /*
                G4cout
        << sec->GetCreatorProcess()->GetProcessName()
        << G4endl;
*/
}
/*
if(totPhoton-lastPrint>=1000)
{
    G4cout
    << "Yield = "
    << totPhoton/(totEdep/MeV)
    << " photons/MeV"
    << G4endl;

    lastPrint = totPhoton;
}*/
/*
for(const auto* sec : *secondaries)
{
    if(sec->GetDefinition()
       == G4OpticalPhoton::Definition())
    {
        nCreated++;
    }
}

if(nCreated - lastPrint >= 100)
{
    G4cout
      << "Generated optical photons = "
      << nCreated
      << G4endl;

    lastPrint = nCreated;
}
*/
if(track->GetDefinition()==G4OpticalPhoton::Definition())
{
    if(track->GetCurrentStepNumber()==1)
    {
        static long long nPhoton=0;
        nPhoton++;
/*
        if(nPhoton%1000==0)
        {
            G4cout<<"Photon track = "<<nPhoton<<G4endl;
        }*/
    }
}

// if (track->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
 //   if (track->GetGlobalTime() > 20*ns) {
  //    track->SetTrackStatus(fStopAndKill);
//    } return;
//}

auto prePV = step->GetPreStepPoint()->GetPhysicalVolume();
auto postPV = step->GetPostStepPoint()->GetPhysicalVolume();
/*
if(prePV && prePV->GetName()=="fiberCore_Phys")
{
    coreEdep += step->GetTotalEnergyDeposit();
}
*/
if(step->GetTrack()->GetDefinition()
   == G4OpticalPhoton::Definition())
{
static int cathEnter = 0;

if(prePV &&
   postPV &&
   prePV->GetName()=="fiberCore_Phys" &&
   postPV->GetName().find("ModuleC_Cath")
   != std::string::npos)
{
    cathEnter++;

    if(cathEnter%100==0)
    {
        G4cout
        << "Cath enter = "
        << cathEnter
        << G4endl;
    }
}

    auto postPV =
      step->GetPostStepPoint()
      ->GetPhysicalVolume();

    if(!postPV) return;

    G4OpBoundaryProcess* boundary = nullptr;

    auto pm =
      step->GetTrack()
      ->GetDefinition()
      ->GetProcessManager();

    auto postStep =
      pm->GetPostStepProcessVector(typeDoIt);

    for(size_t i=0; i<postStep->size(); i++)
    {
        auto proc =
          (*postStep)[i];

        if(proc->GetProcessName()
           == "OpBoundary")
        {
            boundary =
              dynamic_cast<G4OpBoundaryProcess*>(proc);
            break;
        }
    }

static G4int detCount = 0;

if(boundary)
{
    auto status =
      boundary->GetStatus();

    if(status == Detection)
    {
        detCount++;

        auto postLV =
          postPV->GetLogicalVolume();

        auto sd =
          postLV->GetSensitiveDetector();
/*
        G4cout
        << "DETECTION #"
        << detCount
        << " at "
        << postPV->GetName()
        << " SD=" << sd
        << G4endl;
*/
        if(sd)
        {
            sd->Hit(
              const_cast<G4Step*>(step)
            );
        }
    }
}
}
  if (step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) return;

//  G4Track* track = step->GetTrack();
  G4ParticleDefinition* particle = track->GetDefinition();
  G4int pdgID = particle->GetPDGEncoding();

  G4StepPoint* presteppoint = step->GetPreStepPoint();
  G4StepPoint* poststeppoint = step->GetPostStepPoint();
  G4LogicalVolume* preVol = presteppoint->GetPhysicalVolume()->GetLogicalVolume();
  G4TouchableHandle theTouchable = presteppoint->GetTouchableHandle();
  G4String matName = preVol->GetMaterial()->GetName();
  G4VPhysicalVolume* motherTower = GetMotherTower(theTouchable);


//  if (poststeppoint->GetStepStatus() == fWorldBoundary) {
    if ( matName=="Polystyrene" ) {
    fLeak.E = track->GetTotalEnergy();
    fLeak.px = track->GetMomentum().x();
    fLeak.py = track->GetMomentum().y();
    fLeak.pz = track->GetMomentum().z();
    fLeak.vx = presteppoint->GetPosition().x();
    fLeak.vy = presteppoint->GetPosition().y();
    fLeak.vz = presteppoint->GetPosition().z();
    fLeak.vt = presteppoint->GetGlobalTime();
    fLeak.pdgId = track->GetDefinition()->GetPDGEncoding();
    fLeak.EdepCore = step->GetTotalEnergyDeposit();
    fLeak.ModuleNum = GetModuleNum(motherTower->GetName());

    fEventAction->fillLeaks(fLeak);
  }

if (step->GetTrack()->GetDefinition() ==
    G4OpticalPhoton::Definition())

//  G4String matName = preVol->GetMaterial()->GetName();

  if ( matName=="G4_Galactic" || matName=="Air" ) return;

if (step->GetTrack()->GetDefinition() ==
    G4OpticalPhoton::Definition())
{
    auto proc =
      poststeppoint->GetProcessDefinedStep();

    if(proc &&
       proc->GetProcessName()=="OpBoundary")
    {
        auto boundary =
        (G4OpBoundaryProcess*)proc;

        G4cout
        << boundary->GetStatus()
        << G4endl;
    }
}

  G4ThreeVector globalPos = presteppoint->GetPosition();
  
  // 상황에 따라 Depth(보통 1 또는 2)를 조절하여 최상위 ModuleLogical 기준의 로컬 좌표를 얻습니다.
  G4ThreeVector localPos = theTouchable->GetHistory()->GetTransform(1).TransformPoint(globalPos); 
  G4double x = localPos.x();
  G4double y = localPos.y();
  G4double z = localPos.z(); // -totalLength/2 ~ +totalLength/2 범위

  // 2. 보내주신 코드의 지오메트리 상수 계산 정의
  G4double rmin = 904.485 * mm; 
  G4double rmax = 1037.675 * mm + (20.889 * 7) * mm; 
  G4double totalLength = rmax - rmin; // 전체 Z축 길이
  G4double pDz = totalLength / 2.0;   // Z Half-length

  G4double dPhi = (360. / 48.) * deg;
  G4double pDx1 = rmin * std::tan(dPhi/2.); // 로컬 z = -pDz 에서의 x half-width
  G4double pDx2 = rmax * std::tan(dPhi/2.); // 로컬 z = +pDz 에서의 x half-width

  // 3. Z축 로컬 좌표 기준 가상 레이어 인덱스(0 ~ 9) 판별
  // 로컬 z는 -pDz ~ +pDz 이므로, 계산 편의를 위해 0 ~ totalLength 범위로 시프트
  G4double z_shifted = z + pDz; 
  
  G4int layerIdx = -1;
  G4double zStart = 0, zEnd = 0, thick = 0;
  G4double curZ = 0;

  // [A] SFIL 구역: 3개 가상 레이어 반복 스캔
  for(int r = 0; r < 3; r++) {
      curZ += 17.0*mm; // 앞선 Vac 레이어 두께 점프
      thick = 21.73*mm; // 가상 Pb 레이어 두께
      if (z_shifted >= curZ && z_shifted <= curZ + thick) {
          layerIdx = r;
          zStart = curZ;
          zEnd = zStart + thick;
          break;
      }
      curZ += thick;
  }

  // [B] BULK 구역: 7개 가상 레이어 반복 스캔 (SFIL 끝난 지점 이후)
  if (layerIdx == -1) {
      curZ += 17.0*mm; // SFIL 영역의 마지막 Vac 레이어 점프
      for(int r = 0; r < 7; r++) {
          thick = 20.889*mm; // 가상 Bulk 레이어 두께
          if (z_shifted >= curZ && z_shifted <= curZ + thick) {
              layerIdx = 3 + r;
              zStart = curZ;
              zEnd = zStart + thick;
              break;
          }
          curZ += thick;
      }
  }

  // 만약 입자가 가상 레이어 경계 밖(진공 레이어 구간 등)에 찍혔다면 리턴 처리
  if (layerIdx == -1) return;

  // 4. X축 좌표를 기준으로 지그재그 5개 사다리꼴 섹션(0 ~ 4) 판별
  // 해당 가상 레이어의 시작점(zStart)과 끝점(zEnd)에서의 전체 엄마 사다리꼴 너비 계산
  G4double fullDx1 = pDx1 + (pDx2 - pDx1) * (zStart / totalLength);
  G4double fullDx2 = pDx1 + (pDx2 - pDx1) * (zEnd / totalLength);
  
  G4double delta = fullDx2 - fullDx1;
  G4double w_avg = (fullDx1 + fullDx2) / 10.0; // 5분할 조각의 평균 Half-width

  // 현재 가상 레이어 내부에서의 상대적인 Z축 위치 비율 (0.0 ~ 1.0)
  G4double z_rel = (z_shifted - zStart) / thick;
  
  // 180도 지그재그 회전 배치를 모사하기 위한 빗각 경계면 오프셋
  G4double offset = delta * (z_rel - 0.5);

  // DetectorConstruction의 5분할 간격식(posX = (2*s - 4)*w_avg)과 
  // 지그재그 교차 치수 대입법을 역산한 4개의 X 경계선
  G4double bnd1 = -3.0 * w_avg + offset;
  G4double bnd2 = -1.0 * w_avg - offset;
  G4double bnd3 =  1.0 * w_avg + offset;
  G4double bnd4 =  3.0 * w_avg - offset;

  G4int sectionIdx = 0;
  if      (x < bnd1) sectionIdx = 0;
  else if (x < bnd2) sectionIdx = 1;
  else if (x < bnd3) sectionIdx = 2;
  else if (x < bnd4) sectionIdx = 3;
  else               sectionIdx = 4;

  // 5. 최종 가상 모듈 번호 (0 ~ 49) 저장
  fEdep.ModuleNum = (layerIdx * 5) + sectionIdx;

  G4double pdgCharge = particle->GetPDGCharge();

  fEdep.Edep = step->GetTotalEnergyDeposit();
  fEdep.EdepEle = (std::abs(pdgID)==11) ? fEdep.Edep : 0.;
  fEdep.EdepGamma = (std::abs(pdgID)==22) ? fEdep.Edep : 0.;
  fEdep.EdepCharged = ( std::round(std::abs(pdgCharge)) != 0. ) ? fEdep.Edep : 0.;

  if ( fEdep.Edep > 0. ) fEventAction->fillEdeps(fEdep);
//  }
  return;
}
