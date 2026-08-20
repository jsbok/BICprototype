#include "koBICPrimaryGeneratorAction.hh"
#include "koBICRunAction.hh"

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"
#include "Randomize.hh"
// #include "PhysicalConstants.h"
#include <cmath>

namespace { G4Mutex koBICPrimaryGeneratorMutex = G4MUTEX_INITIALIZER; }
int koBICPrimaryGeneratorAction::sNumEvt = 0;
G4ThreadLocal int koBICPrimaryGeneratorAction::sIdxEvt = 0;

using namespace std;
koBICPrimaryGeneratorAction::koBICPrimaryGeneratorAction(G4int seed, G4bool useHepMC, G4bool useCalib, G4bool useGPS)
: G4VUserPrimaryGeneratorAction()
{
  fSeed = seed;
  fUseHepMC = useHepMC;
  fUseCalib = useCalib;
  fUseGPS = useGPS;

  if (!fUseHepMC) {
    if (fUseGPS) initGPS();
    else initPtcGun();
  }
}

void koBICPrimaryGeneratorAction::initPtcGun() {
  fTheta = -0.01111;
  fPhi = 0.;
  fRandX = 0.*mm;
  fRandZ = 10.*mm;
  fX_0 = 0.;
  fY_0 = 0.;
  fZ_0 = -1000.;
  fParticleGun = new G4ParticleGun(1);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  fElectron = particleTable->FindParticle(particleName="e-");
  fPositron = particleTable->FindParticle(particleName="e+");
  fMuon = particleTable->FindParticle(particleName="mu+");
  fPion = particleTable->FindParticle(particleName="pi+");
  fKaon0L = particleTable->FindParticle(particleName="kaon0L");
  fProton = particleTable->FindParticle(particleName="proton");
  fOptGamma = particleTable->FindParticle(particleName="opticalphoton");

  // define commands for this class
  DefineCommands();
}

void koBICPrimaryGeneratorAction::initGPS() {
  fGPS = new G4GeneralParticleSource();
}

koBICPrimaryGeneratorAction::~koBICPrimaryGeneratorAction() {
  if (!fUseHepMC) {
    if (fUseGPS) delete fGPS;
    else {
      if (fParticleGun) delete fParticleGun;
      if (fMessenger) delete fMessenger;
    }
  }
}

void koBICPrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {

  if (fUseGPS) {
    G4AutoLock lock(&koBICPrimaryGeneratorMutex);
    fGPS->GeneratePrimaryVertex(event);
    sIdxEvt = sNumEvt;
    sNumEvt++;

    return;
  }

  G4double x = (G4UniformRand()-0.5)*fRandX + fX_0;
  G4double y = (G4UniformRand()-0.5)*fRandZ + fY_0;
  G4double z = 0 + fZ_0;
  fOrg.set(x,y,z);

  fParticleGun->SetParticlePosition(fOrg); // http://www.apc.univ-paris7.fr/~franco/g4doxy/html/classG4VPrimaryGenerator.html

  fDirection.setREtaPhi(1.,0.,0.);
  fDirection.rotateY( -M_PI * ((90 - fTheta)/180) );
  fDirection.rotateX( M_PI * (fPhi/180.) );

  fParticleGun->SetParticleMomentumDirection(fDirection);

  G4AutoLock lock(&koBICPrimaryGeneratorMutex);
  fParticleGun->GeneratePrimaryVertex(event);
  sIdxEvt = sNumEvt;
  sNumEvt++;
}

void koBICPrimaryGeneratorAction::DefineCommands() {
  // Define /koBIC/generator command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/koBIC/generator/", "Primary generator control");

  G4GenericMessenger::Command& etaCmd = fMessenger->DeclareMethodWithUnit("theta","rad",&koBICPrimaryGeneratorAction::SetTheta,"theta of beam");
  etaCmd.SetParameterName("theta",true);
  etaCmd.SetDefaultValue("0.");

  G4GenericMessenger::Command& phiCmd = fMessenger->DeclareMethodWithUnit("phi","rad",&koBICPrimaryGeneratorAction::SetPhi,"phi of beam");
  phiCmd.SetParameterName("phi",true);
  phiCmd.SetDefaultValue("0.");

  G4GenericMessenger::Command& x0Cmd = fMessenger->DeclareMethodWithUnit("x0","cm",&koBICPrimaryGeneratorAction::SetX0,"x_0 of beam");
  x0Cmd.SetParameterName("x0",true);
  x0Cmd.SetDefaultValue("0.");

  G4GenericMessenger::Command& y0Cmd = fMessenger->DeclareMethodWithUnit("y0","cm",&koBICPrimaryGeneratorAction::SetY0,"y_0 of beam");
  y0Cmd.SetParameterName("y0",true);
  y0Cmd.SetDefaultValue("0.");

  G4GenericMessenger::Command& z0Cmd = fMessenger->DeclareMethodWithUnit("z0","cm",&koBICPrimaryGeneratorAction::SetZ0,"z_0 of beam");
  z0Cmd.SetParameterName("z0",true);
  z0Cmd.SetDefaultValue("0.");

  G4GenericMessenger::Command& randxCmd = fMessenger->DeclareMethodWithUnit("randx","mm",&koBICPrimaryGeneratorAction::SetRandX,"x width of beam");
  randxCmd.SetParameterName("randx",true);
  randxCmd.SetDefaultValue("0.");

  G4GenericMessenger::Command& randzCmd = fMessenger->DeclareMethodWithUnit("randz","mm",&koBICPrimaryGeneratorAction::SetRandZ,"z width of beam");
  randzCmd.SetParameterName("randz",true);
  randzCmd.SetDefaultValue("0.");
}
