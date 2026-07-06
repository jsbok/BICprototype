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

G4String volName = preVol->GetName();
  G4ThreeVector globalPos = presteppoint->GetPosition();
  
  G4ThreeVector localPos = theTouchable->GetHistory()->GetTransform(1).TransformPoint(globalPos); 
  G4double x = localPos.x();
  G4double y = localPos.y();
  G4double z = localPos.z(); // -totalLength/2 ~ +totalLength/2 범위
  
G4int finalModuleNum = -1;


//if (volName.contains("Module") || volName.contains("Box") || volName.contains("Pb")) {
    
G4double rmin = 825.805 * mm;
G4double rmax = 1036.455 * mm + 180.0 * mm; 
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
    finalModuleNum = 10 + (layerIdx * 3) + colIdx;
}
//}

// =========================================================================
// SFIL (0~9)
// =========================================================================
if (volName.contains("Pb_SFIL") || volName.contains("SFIL")) {
    if (z_shifted < boxStartZ) {
        
        G4int layerIdx = (z < -30.0 * mm) ? 0 : 1; 

        G4double dPhi = (360. / 48.) * deg;
        G4double pDx1 = rmin * std::tan(dPhi/2.);
        G4double pDx2 = rmax * std::tan(dPhi/2.);
        G4double thick = 21.73 * mm;
        
        G4double zStart = (layerIdx == 1) ? (totalLength - 180.0 - 17.0 - 21.73) * mm 
                                          : (totalLength - 180.0 - 17.0 - 21.73 - 17.0 - 21.73) * mm; 
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

if (finalModuleNum < 0 || finalModuleNum > 27) return; 

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
