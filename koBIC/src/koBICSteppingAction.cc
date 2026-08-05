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

//  if ( matName=="G4_Galactic" || matName=="Air" ) return;

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

G4String volName = preVol->GetName();
  G4ThreeVector globalPos = presteppoint->GetPosition();
  
  G4ThreeVector localPos = theTouchable->GetHistory()->GetTransform(1).TransformPoint(globalPos); 
  G4double x = localPos.x();
  G4double y = localPos.y();
  G4double z = localPos.z(); // -totalLength/2 ~ +totalLength/2 Range
  
G4int finalModuleNum = -1;


G4int layers = 3;
/*
G4double rmin = 825.805 * mm;
G4double rmax = 1036.455 * mm + 180.0 * mm + (38.73 * layers) * mm; 
G4double totalLength = rmax - rmin;
G4double pDz = totalLength / 2.0;
G4double z_shifted = z + pDz; // 
G4double boxStartZ = (totalLength - 180.0) * mm;// 

if (volName.contains("Module") || volName.contains("Box") || volName.contains("Pb")) {
  */  
G4double rmin = 825.805 * mm;
G4double rmax = 1036.455 * mm + 180.0 * mm + (38.73 * layers) * mm; 
G4double totalLength = rmax - rmin;
G4double pDz = totalLength / 2.0;
G4double z_shifted = z + pDz; // 

G4double boxStartZ = (totalLength - 180.0) * mm;// 

G4double zStart = pDz - 180.0 * mm;
G4double zEnd   = pDz;

if (z >= zStart && z <= zEnd) {
    G4int layerIdx = std::floor((z - zStart) / (30.0 * mm));
    if (layerIdx < 0) layerIdx = 0;
    if (layerIdx > 5) layerIdx = 5;


    G4int colIdx = 1; // (0mm)
    if (x < -15.0 * mm) {
        colIdx = 0; // -30mm 
    } else if (x > 15.0 * mm) {
	    colIdx = 2; // +30mm 
    }
    finalModuleNum = 10 + 5*layers + (layerIdx * 3) + colIdx;
}
//}

// =========================================================================
// SFIL (1~5)
// =========================================================================
// 💡 1. "Pb_SFIL"만 엄격하게 검사 (Vac_SFIL 진공 볼륨 제외)
if (volName.contains("Pb_SFIL")) {
    if (z_shifted < boxStartZ) {
        
        G4double refZ = totalLength - (30.0 * 6.0 * mm) - (38.73 * 3.0 * mm); 
        G4double thick = 21.24 * mm; // SFIL 두께

        G4double zStart_SFIL3 = refZ;
        G4double zStart_SFIL2 = zStart_SFIL3 - 38.73 * mm;
        G4double zStart_SFIL1 = zStart_SFIL2 - 38.24 * mm; // 17mm 간격 반영

        G4int layerIdx = -1;
        G4double zStart = 0.0;

        // 💡 2. 각 레이어의 실제 두께 범위 [zStart, zStart + thick] 내에 들어왔는지 엄격히 검증
        if (z_shifted >= zStart_SFIL1 && z_shifted <= zStart_SFIL1 + thick) {
            layerIdx = 0; zStart = zStart_SFIL1;
        } else if (z_shifted >= zStart_SFIL2 && z_shifted <= zStart_SFIL2 + thick) {
            layerIdx = 1; zStart = zStart_SFIL2;
        } else if (z_shifted >= zStart_SFIL3 && z_shifted <= zStart_SFIL3 + thick) {
            layerIdx = 2; zStart = zStart_SFIL3;
        } else if (z_shifted >= (zStart_SFIL3 + 38.73*mm) && z_shifted <= (zStart_SFIL3 + 38.73*mm) + thick) {
            layerIdx = 3; zStart = zStart_SFIL3 + 38.73*mm;
        } else if (z_shifted >= (zStart_SFIL3 + 77.46*mm) && z_shifted <= (zStart_SFIL3 + 77.46*mm) + thick) {
            layerIdx = 4; zStart = zStart_SFIL3 + 77.46*mm;
        }

        // 💡 레이어 범위 밖(레이어 사이의 진공 Gap 등)이라면 즉시 탈락
        if (layerIdx == -1) return;

        // --- 이하 X축 5등분 사다리꼴 매핑 수식 (이전과 동일) ---
        G4double dPhi = (360. / 48.) * deg;
        G4double pDx1 = rmin * std::tan(dPhi/2.);
        G4double pDx2 = rmax * std::tan(dPhi/2.);
        
        G4double zEnd = zStart + thick;

        G4double fullDx1 = pDx1 + (pDx2 - pDx1) * (zStart / totalLength);
        G4double fullDx2 = pDx1 + (pDx2 - pDx1) * (zEnd / totalLength);
        
        G4double delta = fullDx2 - fullDx1;
        G4double w_avg = (fullDx1 + fullDx2) / 10.0;
        
        G4double z_rel_sfil = (z_shifted - zStart) / thick;
        G4double offset = delta * (z_rel_sfil - 0.5);

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

        finalModuleNum = (layerIdx * 5) + sectionIdx;
    }
}

if (finalModuleNum < 0 || finalModuleNum > 42) return; 

// =========================================================================
// [진공/여백 오매핑 검증 및 시뮬레이션 중단 로직]
// =========================================================================
/*
// 1. 현재 입자가 위치한 영역의 실제 물리적 물질(Material) 
G4Material* currentMat = presteppoint->GetMaterial();
G4String matName2 = currentMat->GetName();

// 2. 진공(Vacuum/Air) 물질인지 판별 
bool isVacuum = (matName2.contains("AIR"));

// 3. 모듈 번호(0~42)는 정상 할당되었으나, 실제 물질이 진공인 경우 -> 시뮬레이션 즉시 폭파
if (finalModuleNum >= 0 && finalModuleNum <= 42 && isVacuum) {
    
    // 터미널에 출력할 에러 메시지 
    std::string errorMsg = "\n=========================================================\n";
    errorMsg += " [GEOMETRY MISMATCH DETECTED - SIMULATION ABORTED]\n";
    errorMsg += "  A Vacuum region was wrongly mapped as a valid Module ID!\n";
    errorMsg += "---------------------------------------------------------\n";
    errorMsg += "  * Mapped Module Num : " + std::to_string(finalModuleNum) + "\n";
    errorMsg += "  * Actual Volume Name: " + std::string(volName) + "\n";
    errorMsg += "  * Actual Material   : " + std::string(matName2) + "\n";
    errorMsg += "  * Step Position (X) : " + std::to_string(x/mm) + " mm\n";
    errorMsg += "  * Step Position (Y) : " + std::to_string(y/mm) + " mm\n";
    errorMsg += "  * Position (z_shift): " + std::to_string(z_shifted/mm) + " mm\n";
    errorMsg += "=========================================================\n";

    // G4Exception을 호출하여 터미널을 에러 로그와 함께 즉시 강제 종료.
    G4Exception("SteppingAction::UserSteppingAction()", 
                "ERR_VACUUM_MAPPED_AS_MODULE", 
                FatalException, 
                errorMsg.c_str());
}
*/
fEdep.ModuleNum = finalModuleNum;

G4double pdgCharge = particle->GetPDGCharge();
fEdep.Edep = step->GetTotalEnergyDeposit();
fEdep.EdepEle = (std::abs(pdgID) == 11) ? fEdep.Edep : 0.;
fEdep.EdepGamma = (std::abs(pdgID) == 22) ? fEdep.Edep : 0.;
fEdep.EdepCharged = (std::round(std::abs(pdgCharge)) != 0.) ? fEdep.Edep : 0.;

if (fEdep.Edep > 0.) {
    fEventAction->fillEdeps(fEdep);
}
return;
}
