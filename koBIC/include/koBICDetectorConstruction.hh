#ifndef koBICDetectorConstruction_h
#define koBICDetectorConstruction_h 1

#include "koBICMagneticField.hh"
#include "koBICMaterials.hh"
#include "koBICSiPMHit.hh"

#include "G4Box.hh"
#include "G4FieldManager.hh"
#include "G4GenericMessenger.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VSolid.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4Region.hh"

#include "dimensionCalc.hh"

using namespace std;

class koBICMagneticField;

class koBICDetectorConstruction : public G4VUserDetectorConstruction {
public:
  koBICDetectorConstruction();
  virtual ~koBICDetectorConstruction();

  virtual G4VPhysicalVolume *Construct();
  virtual void ConstructSDandField();

  static int fNofModules;
  static int fNofRow;
  static int fNofCol;

private:
  void DefineCommands();
  void DefineMaterials();
  G4Material *FindMaterial(G4String matName) {
    return fMaterials->GetMaterial(matName);
  }
  G4OpticalSurface *FindSurface(G4String surfName) {
    return fMaterials->GetOpticalSurface(surfName);
  }

  void ModuleBuild(
      G4LogicalVolume *ModuleLogical_[], G4LogicalVolume *PMTGLogical_[],
      G4LogicalVolume *PMTfilterLogical_[], G4LogicalVolume *PMTcellLogical_[],
      G4LogicalVolume *PMTcathLogical_[],
      std::vector<koBICInterface::koBICModuleProperty> &towerProps_);

  void FiberImplement(G4int i, 
  		      std::vector<G4LogicalVolume*> logicPbSFIL, 
                      std::vector<G4LogicalVolume*> logicPbBulk
                      );

  G4bool checkOverlaps;
  G4GenericMessenger *fMessenger;
  koBICMaterials *fMaterials;

  static G4ThreadLocal koBICMagneticField *fMagneticField;
  static G4ThreadLocal G4FieldManager *fFieldMgr;

  G4VisAttributes *fVisAttrOrange;
  G4VisAttributes *fVisAttrBlue;
  G4VisAttributes *fVisAttrGray;
  G4VisAttributes *fVisAttrGreen;
  G4VisAttributes *fVisAttrSkyBlue;

  G4Region* fScintRegion;
  G4Region* fCerenRegion;

  G4double fFrontL;
  G4double fTowerDepth;
  G4double fModuleH;
  G4double fModuleW;
  G4double fFiberUnitH;
  G4int fRandomSeed;
  G4double fFiber_vert_dis;
  G4double fFiber_hori_dis;

  G4double PMTT;
  G4double filterT;

  G4bool doFiber;
  G4bool doPMT;
  G4bool doGlue;

  dimensionCalc *dimCalc;

  char name[20];
  G4String moduleName;
  G4Box *module;
  G4Box *pmtg;
  G4Box *pmtg1;
  G4Box *pmtcath;

  G4Box *fiberUnit;
  G4Tubs *fiberClad_SFIL;
  G4Tubs *fiberCore_SFIL;
  G4Tubs *fiberGlue_SFIL;
  G4Tubs *fiberClad_Bulk;
  G4Tubs *fiberCore_Bulk;
  G4Tubs *fiberGlue_Bulk;
  G4Box *gluebox;

  G4VSolid *tfiberUnitIntersection;
  G4VSolid *tfiberCladIntersection;
  G4VSolid *tfiberCoreIntersection;
  G4VSolid *tGlueSubtraction;
  G4VSolid *tGlueIntersection;

  G4LogicalVolume *ModuleLogical[100];

  G4LogicalVolume *PMTGLogical[100];
  G4LogicalVolume *PMTcathLogical[100];
  G4LogicalVolume *PMTcellLogical[100];
  G4LogicalVolume *PMTfilterLogical[100];

  vector<G4LogicalVolume *> fiberUnitIntersection[100];
  vector<G4LogicalVolume *> fiberCladIntersection[100];
  vector<G4LogicalVolume *> fiberCoreIntersection[100];
  vector<G4LogicalVolume *> glueIntersection__[100];

  koBICInterface::hitXY fTowerXY;
  std::vector<koBICInterface::koBICModuleProperty> fModuleProp;

  // G4double clad_S_rMin;
  G4double clad_S_rMax;
  // G4double clad_S_Dz  ;
  // G4double clad_S_Sphi;
  // G4double clad_S_Dphi;

  // G4double core_S_rMin;
  G4double core_S_rMax;
  // G4double core_S_Dz  ;
  // G4double core_S_Sphi;
  // G4double core_S_Dphi;
  G4double glue_S_rMax;
  
  std::vector<G4float> fFiberX;
  std::vector<G4float> fFiberY;
  std::vector<G4float> fFiberZ;
  std::vector<G4bool> fFiberWhich;

  G4LogicalVolume *worldLogical;

  G4String setModuleName(int i) { return std::to_string(i); }
};

#endif
