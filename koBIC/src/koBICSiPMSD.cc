#include "koBICSiPMSD.hh"
#include "koBICSiPMHit.hh"
#include "koBICDetectorConstruction.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

using namespace std;

koBICSiPMSD::koBICSiPMSD(const G4String& name, const G4String& hitsCollectionName, const G4int& isLeft, koBICInterface::koBICModuleProperty ModuleProp, G4int modNum)
: G4VSensitiveDetector(name), fHitCollection(0), fHCID(-1), fWavBin(60), fTimeBin(700),
fModuleNum(-1), fWavlenStart(900.), fWavlenEnd(300.), fTimeStart(0.), fTimeEnd(70.)
{
  collectionName.insert(hitsCollectionName);
  fWavlenStep = (fWavlenStart-fWavlenEnd)/(float)fWavBin;
  fTimeStep = (fTimeEnd-fTimeStart)/(float)fTimeBin;

  fModuleNum = ModuleProp.ModuleNum;
//  fTowerXY = ModuleProp.towerXY;
  fisLeft = isLeft;
}

koBICSiPMSD::~koBICSiPMSD() {}

void koBICSiPMSD::Initialize(G4HCofThisEvent* hce) {
  fHitCollection = new koBICSiPMHitsCollection(SensitiveDetectorName,collectionName[0]);
  if (fHCID<0) { fHCID = GetCollectionID(0); }
  hce->AddHitsCollection(fHCID,fHitCollection);
}

G4bool koBICSiPMSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
 /*   G4cout << "SiPM HIT : "
           << step->GetTrack()->GetDefinition()->GetParticleName()
           << " in "
           << step->GetPreStepPoint()
                   ->GetTouchableHandle()
                   ->GetVolume()
                   ->GetName()
           << G4endl;
*/
  if(step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) return false;
  G4VPhysicalVolume* currentVol = step->GetPostStepPoint()->GetTouchable()->GetVolume(0);
    G4String volName = currentVol->GetName();
  G4StepPoint* postStep = step->GetPostStepPoint();
G4int copy0 = postStep->GetTouchable()->GetCopyNumber(0);
G4int copy1 = postStep->GetTouchable()->GetCopyNumber(1);

/*
if (volName.find("Cell") == std::string::npos) {
        return false; // Cell이 아니면(Envelope이면) 기록하지 않고 종료
    }
*/
/*
G4cout << "Hit! Volume: " << volName 
       << " | Copy0: " << copy0 
       << " | Copy1: " << copy1 
       << " | ModuleID: " << fModuleNum << G4endl;
  */
  auto touchable = step->GetPostStepPoint()->GetTouchable();
  G4int SiPMnum = step->GetPostStepPoint()->GetTouchable()->GetVolume(0)->GetCopyNo();

  G4int sectionIdx = touchable->GetReplicaNumber(1); 
  G4int layerIdx = touchable->GetReplicaNumber(2); 

  G4int realModuleNum = (layerIdx * 5) + sectionIdx;

  G4int nofHits = fHitCollection->entries();
  G4double hitTime = step->GetPostStepPoint()->GetGlobalTime();
  G4double energy = step->GetTrack()->GetTotalEnergy();

  koBICSiPMHit* hit = NULL;

  for (G4int i = 0; i < nofHits; i++) {

    // G4cout << " Hit iteration : " 
    //        << i << " " 
    //        << SiPMnum << " " 
    //        << (*fHitCollection)[i]->GetSiPMnum() << " " 
    //        << fModuleNum << " " 
    //        << (*fHitCollection)[i]->GetModuleNum() << " "
    //        << G4endl;
    
if ( (*fHitCollection)[i]->GetSiPMnum() == SiPMnum && (*fHitCollection)[i]->GetModuleNum() == fModuleNum && (*fHitCollection)[i]->GetisLeft() == fisLeft) {
      hit = (*fHitCollection)[i];
      break;
  }
  }

  if (hit==NULL) {
    hit = new koBICSiPMHit(fWavBin,fTimeBin);
    hit->SetSiPMnum(SiPMnum);
    hit->SetModuleNum(fModuleNum);
    hit->SetisLeft(fisLeft);
 //   hit->SetTowerXY(fTowerXY);
 //   hit->SetSiPMXY(findSiPMXY(SiPMnum,fTowerXY));
    hit->SetSiPMpos(step->GetPostStepPoint()->GetTouchableHandle()->GetHistory()->GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0.,0.,0.)));

    fHitCollection->insert(hit);
  }

  hit->photonCount();
  ///G4cout<<"photoncount" <<hit->GetPhotonCount()<< G4endl;

  koBICInterface::hitRange wavRange = findWavRange(energy);
  hit->CountWavlenSpectrum(wavRange);

  koBICInterface::hitRange timeRange = findTimeRange(hitTime);
  hit->CountTimeStruct(timeRange);

  return true;
}

void koBICSiPMSD::EndOfEvent(G4HCofThisEvent*) {
  if ( verboseLevel>1 ) {
    G4int nofHits = fHitCollection->entries();
    G4cout
    << G4endl
    << "-------->Hits Collection: in this event they are " << nofHits
    << " hits in the tracker chambers: " << G4endl;
    for ( G4int i=0; i<nofHits; i++ ) (*fHitCollection)[i]->Print();
  }
}

koBICInterface::hitRange koBICSiPMSD::findWavRange(G4double en) {
  int i = 0;
  for ( ; i < fWavBin+1; i++) {
    if ( en < wavToE( (fWavlenStart - (float)i*fWavlenStep)*nm ) ) break;
  }

  if (i==0) return std::make_pair(fWavlenStart,99999.);
  else if (i==fWavBin+1) return std::make_pair(0.,fWavlenEnd);

  return std::make_pair( fWavlenStart-(float)i*fWavlenStep, fWavlenStart-(float)(i-1)*fWavlenStep );
}

koBICInterface::hitRange koBICSiPMSD::findTimeRange(G4double stepTime) {
  int i = 0;
  for ( ; i < fTimeBin+1; i++) {
    if ( stepTime < ( (fTimeStart + (float)i*fTimeStep)*ns ) ) break;
  }

  if (i==0) return std::make_pair(0.,fTimeStart);
  else if (i==fTimeBin+1) return std::make_pair(fTimeEnd,99999.);

  return std::make_pair( fTimeStart+(float)(i-1)*fTimeStep, fTimeStart+(float)i*fTimeStep );
}

koBICInterface::hitXY koBICSiPMSD::findSiPMXY(G4int SiPMnum, koBICInterface::hitXY towerXY) {
  int x = 0; //SiPMnum/(towerXY.second*2-1)*2 + (SiPMnum%(towerXY.second*2-1))/towerXY.second;
  int y = 0; //(SiPMnum%(towerXY.second*2-1))%towerXY.second;

  return std::make_pair(x,y);
}
