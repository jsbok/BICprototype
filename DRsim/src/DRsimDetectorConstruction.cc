#include "DRsimDetectorConstruction.hh"
#include "DRsimCellParameterisation.hh"
#include "DRsimFilterParameterisation.hh"
#include "DRsimMirrorParameterisation.hh"
#include "DRsimSiPMSD.hh"

#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"

#include "G4GeometryManager.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"

#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <string>

using namespace std;

G4ThreadLocal DRsimMagneticField *DRsimDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager *DRsimDetectorConstruction::fFieldMgr = 0;

int DRsimDetectorConstruction::fNofRow = 4;
int DRsimDetectorConstruction::fNofCol = 4;
int DRsimDetectorConstruction::fNofModules = fNofRow * fNofCol;

DRsimDetectorConstruction::DRsimDetectorConstruction()
    : G4VUserDetectorConstruction(), fMessenger(0), fMaterials(NULL) {
  DefineCommands();
  DefineMaterials();

  //   clad_S_rMin = 0.485 * mm;
  clad_S_rMax = 0.50 * mm;
  // clad_S_Dz   = 2.5*m;
  // clad_S_Sphi = 0.;
  // clad_S_Dphi = 2.*M_PI;

  //   core_S_rMin = 0. * mm;
  core_S_rMax = 0.485 * mm;
  // core_S_Dz   = 2.5*m;
  // core_S_Sphi = 0.;
  // core_S_Dphi = 2.*M_PI;

  PMTT = 0.3 * mm;
  // filterT = 0.01*mm;

  fVisAttrOrange = new G4VisAttributes(G4Colour(1.0, 0.5, 0., 0.7));
  fVisAttrOrange->SetForceSolid(true);
  fVisAttrOrange->SetVisibility(true);
  fVisAttrBlue = new G4VisAttributes(G4Colour(0., 0., 1.0, 0.7));
  fVisAttrBlue->SetForceSolid(true);
  fVisAttrBlue->SetVisibility(true);
  fVisAttrGray = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3, 0.7));
  fVisAttrGray->SetVisibility(true);
  fVisAttrGreen = new G4VisAttributes(G4Colour(0.3, 0.7, 0.3, 0.7));
  fVisAttrGreen->SetVisibility(true);
  // fVisAttrSkyBlue = new G4VisAttributes(G4Colour(0.5, 0.8, 0.9, 0.7));
  // fVisAttrSkyBlue->SetVisibility(true);
}

DRsimDetectorConstruction::~DRsimDetectorConstruction() {
  delete fMessenger;
  delete fMaterials;

  delete fVisAttrOrange;
  delete fVisAttrBlue;
  delete fVisAttrGray;
  delete fVisAttrGreen;
  // delete fVisAttrSkyBlue;
}

void DRsimDetectorConstruction::DefineMaterials() {
  fMaterials = DRsimMaterials::GetInstance();
}

G4VPhysicalVolume *DRsimDetectorConstruction::Construct() {
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  checkOverlaps = false;

  G4VSolid *worldSolid = new G4Box("worldBox", 10. * m, 10. * m, 10. * m);
  worldLogical = new G4LogicalVolume(worldSolid, FindMaterial("G4_Galactic"),
                                     "worldLogical");
  G4VPhysicalVolume *worldPhysical =
      new G4PVPlacement(0, G4ThreeVector(), worldLogical, "worldPhysical", 0,
                        false, 0, checkOverlaps);

  fFrontL =
      1000.; // NOTE :: Length from the center of world box to center of module
  fTowerDepth = 320.;
  fModuleH = 30;
  fModuleW = 30;
  fFiberUnitH = 1.;
  fFiber_vert_dis = 1.22;
  fFiber_hori_dis = 1.35;

  G4double fGlue_thickness = 0.233; // Expected value as in Glue-X. TBD more
                                    // realistic (0.15 ~ 0.30, late June 2024).

  doFiber = true;
  doPMT = true;
  doGlue = true;

  fiberUnit = new G4Box("fiber_SQ", (fFiberUnitH / 2) * mm, (1. / 2) * mm,
                        (fTowerDepth / 2) * mm);
  fiberClad = new G4Tubs("fiber", 0, clad_S_rMax, (fTowerDepth + 6) / 2.,
                         0 * deg, 360. * deg);
  fiberCoreS = new G4Tubs("fiberS", 0, core_S_rMax, (fTowerDepth + 6) / 2.,
                          0 * deg, 360. * deg);
  gluebox = new G4Box("gluebox", (fGlue_thickness / 2) * mm,
                      (fFiber_hori_dis / 2) * mm, fTowerDepth / 2.);
  tGlueSubtraction = new G4SubtractionSolid("glueCladSubt", gluebox, fiberClad, 0, G4ThreeVector(.0, .0, .0));

  dimCalc = new dimensionCalc();
  dimCalc->SetFrontL(fFrontL);
  dimCalc->SetTower_height(fTowerDepth);
  dimCalc->SetPMTT(PMTT);
  dimCalc->SetNofModules(fNofModules);
  dimCalc->SetNofRow(fNofRow);
  dimCalc->SetNofCol(fNofCol);
  dimCalc->SetModule_height(fModuleH);
  dimCalc->SetModule_width(fModuleW);

  ModuleBuild(ModuleLogical, PMTGLogical, PMTfilterLogical, PMTcellLogical,
              PMTcathLogical, fiberUnitIntersection, fiberCladIntersection,
              fiberCoreIntersection, fModuleProp);

  delete dimCalc;
  return worldPhysical;
}

void DRsimDetectorConstruction::ConstructSDandField() {
  G4SDManager *SDman = G4SDManager::GetSDMpointer();
  G4String SiPMName = "SiPMSD";

  // ! Not a memory leak - SDs are deleted by G4SDManager. Deleting them
  // manually will cause double delete!
  if (doPMT) {
    for (int i = 0; i < fNofModules; i++) {
      DRsimSiPMSD *SiPMSDmodule_left = new DRsimSiPMSD(
          "Module" + std::to_string(2 * i), "ModuleC" + std::to_string(2 * i),
          0, fModuleProp.at(i));
      DRsimSiPMSD *SiPMSDmodule_right = new DRsimSiPMSD(
          "Module" + std::to_string(2 * i + 1),
          "ModuleC" + std::to_string(2 * i + 1), 1, fModuleProp.at(i));
      SDman->AddNewDetector(SiPMSDmodule_left);
      SDman->AddNewDetector(SiPMSDmodule_right);
      PMTcathLogical[2 * i]->SetSensitiveDetector(SiPMSDmodule_left);
      PMTcathLogical[2 * i + 1]->SetSensitiveDetector(SiPMSDmodule_right);
    }
  }
}

void DRsimDetectorConstruction::ModuleBuild(
    G4LogicalVolume *ModuleLogical_[], G4LogicalVolume *PMTGLogical_[],
    G4LogicalVolume *PMTfilterLogical_[], G4LogicalVolume *PMTcellLogical_[],
    G4LogicalVolume *PMTcathLogical_[],
    std::vector<G4LogicalVolume *> fiberUnitIntersection_[],
    std::vector<G4LogicalVolume *> fiberCladIntersection_[],
    std::vector<G4LogicalVolume *> fiberCoreIntersection_[],
    std::vector<DRsimInterface::DRsimModuleProperty> &ModuleProp_) {

  for (int i = 0; i < fNofModules; i++) {
    moduleName = setModuleName(i);

    dimCalc->SetisModule(true);
    module = new G4Box("Module", (fModuleH / 2.) * mm, (fModuleW / 2.) * mm,
                       (fTowerDepth / 2.) * mm);
    ModuleLogical_[i] =
        new G4LogicalVolume(module, FindMaterial("Lead"), moduleName);
    new G4PVPlacement(0, dimCalc->GetOrigin(i), ModuleLogical_[i], moduleName,
                      worldLogical, false, 0, checkOverlaps);

    if (doPMT) {
      dimCalc->SetisModule(false);
      pmtg = new G4Box("PMTG", (fModuleH / 2.) * mm, (fModuleW / 2.) * mm,
                       PMTT / 2. * mm);
      PMTGLogical_[2 * i] =
          new G4LogicalVolume(pmtg, FindMaterial("G4_AIR"), moduleName);
      PMTGLogical_[2 * i + 1] =
          new G4LogicalVolume(pmtg, FindMaterial("G4_AIR"), moduleName);
      new G4PVPlacement(0, dimCalc->GetOrigin_PMTG(2 * i), PMTGLogical_[2 * i],
                        moduleName, worldLogical, false, 0, checkOverlaps);
      new G4PVPlacement(0, dimCalc->GetOrigin_PMTG(2 * i + 1),
                        PMTGLogical_[2 * i + 1], moduleName, worldLogical,
                        false, 0, checkOverlaps);
    }

    FiberImplement(i, ModuleLogical_, fiberUnitIntersection_,
                   fiberCladIntersection_, fiberCoreIntersection_);
    // GlueImplement(i,ModuleLogical_,fiberUnitIntersection_,fiberCladIntersection_,fiberCoreIntersection_);

    DRsimInterface::DRsimModuleProperty ModulePropSingle;
    ModulePropSingle.towerXY = fTowerXY;
    ModulePropSingle.ModuleNum = i;
    ModuleProp_.push_back(ModulePropSingle);

    if (doPMT) {
      G4VSolid *SiPMlayerSolid =
          new G4Box("SiPMlayerSolid", (fModuleH / 2.) * mm,
                    (fModuleW / 2.) * mm, (PMTT / 2.) * mm);
      G4LogicalVolume *SiPMlayerLogical_left = new G4LogicalVolume(
          SiPMlayerSolid, FindMaterial("G4_AIR"), "SiPMlayerLogical");
      G4LogicalVolume *SiPMlayerLogical_right = new G4LogicalVolume(
          SiPMlayerSolid, FindMaterial("G4_AIR"), "SiPMlayerLogical");
      new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), SiPMlayerLogical_left,
                        "SiPMlayerPhysical", PMTGLogical_[2 * i], false, 0,
                        checkOverlaps);
      new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), SiPMlayerLogical_right,
                        "SiPMlayerPhysical", PMTGLogical_[2 * i + 1], false, 0,
                        checkOverlaps);

      G4VSolid *PMTcellSolid = new G4Box("PMTcellSolid", 1.1 / 2. * mm,
                                         1.1 / 2. * mm, PMTT / 2. * mm);
      PMTcellLogical_[2 * i] = new G4LogicalVolume(
          PMTcellSolid, FindMaterial("Glass"), "PMTcellLogical_");
      PMTcellLogical_[2 * i + 1] = new G4LogicalVolume(
          PMTcellSolid, FindMaterial("Glass"), "PMTcellLogical_");

      DRsimCellParameterisation *PMTcellParam = new DRsimCellParameterisation(
          fTowerXY.first, fTowerXY.second, fModuleH, fModuleW);
      G4PVParameterised *PMTcellPhysical_left = new G4PVParameterised(
          "PMTcellPhysical_left", PMTcellLogical_[2 * i], SiPMlayerLogical_left,
          kXAxis, fTowerXY.first * fTowerXY.second - fTowerXY.first / 2,
          PMTcellParam);
      G4PVParameterised *PMTcellPhysical_right = new G4PVParameterised(
          "PMTcellPhysical_right", PMTcellLogical_[2 * i + 1],
          SiPMlayerLogical_right, kXAxis,
          fTowerXY.first * fTowerXY.second - fTowerXY.first / 2, PMTcellParam);

      G4VSolid *PMTcathSolid = new G4Box("PMTcathSolid", 1.1 / 2. * mm,
                                         1.1 / 2. * mm, PMTT / 2. * mm);
      PMTcathLogical_[2 * i] = new G4LogicalVolume(
          PMTcathSolid, FindMaterial("Silicon"), "PMTcathLogical_");
      PMTcathLogical_[2 * i + 1] = new G4LogicalVolume(
          PMTcathSolid, FindMaterial("Silicon"), "PMTcathLogical_");
      new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), PMTcathLogical_[2 * i],
                        "PMTcathPhysical", PMTcellLogical_[2 * i], false, 0,
                        checkOverlaps);
      new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
                        PMTcathLogical_[2 * i + 1], "PMTcathPhysical",
                        PMTcellLogical_[2 * i + 1], false, 0, checkOverlaps);
      new G4LogicalSkinSurface("Photocath_surf_left", PMTcathLogical_[2 * i],
                               FindSurface("SiPMSurf"));
      new G4LogicalSkinSurface("Photocath_surf_right",
                               PMTcathLogical_[2 * i + 1],
                               FindSurface("SiPMSurf"));

      PMTcathLogical_[2 * i]->SetVisAttributes(fVisAttrGreen);
      PMTcathLogical_[2 * i + 1]->SetVisAttributes(fVisAttrGreen);
    }
  }
}

void DRsimDetectorConstruction::DefineCommands() {}

void DRsimDetectorConstruction::FiberImplement(
    G4int i, G4LogicalVolume *ModuleLogical__[],
    std::vector<G4LogicalVolume *> fiberUnitIntersection__[],
    std::vector<G4LogicalVolume *> fiberCladIntersection__[],
    std::vector<G4LogicalVolume *> fiberCoreIntersection__[]) {

  fFiberX.clear();
  fFiberY.clear();
  fFiberWhich.clear();

  int NofPlate = fModuleH / (fFiber_vert_dis);
  int NofFiber = fModuleW / (fFiber_hori_dis);
  fTowerXY = std::make_pair(NofPlate, NofFiber);

  G4bool fWhich = false;
  for (int k = 0; k < NofPlate; k++) {
    for (int j = 0; j < NofFiber; j++) {
      /*
        ? fX : # of plate , fY : # of fiber in the plate
      */
      if (fWhich && j == NofFiber - 1)
        break;
      G4float fX = -fModuleH * mm / 2 + k * fFiber_vert_dis * mm +
                   fFiber_vert_dis / 2 * mm;
      G4float fY = -fModuleW * mm / 2 + j * fFiber_hori_dis * mm +
                   fFiber_hori_dis / 2 * mm;
      if (fWhich)
        fY += fFiber_hori_dis / 2 * mm;
      fFiberX.push_back(fX);
      fFiberY.push_back(fY);
      fFiberWhich.push_back(fWhich);
    }
    fWhich = !fWhich;
  }

  if (doFiber) {
    for (unsigned int fiberId = 0; fiberId < fFiberX.size(); fiberId++) {

      tfiberCladIntersection = new G4IntersectionSolid(
          "fiberClad", fiberClad, module, 0,
          G4ThreeVector(-fFiberX.at(fiberId), -fFiberY.at(fiberId), 0.));
      fiberCladIntersection__[i].push_back(new G4LogicalVolume(
          tfiberCladIntersection, FindMaterial("PMMA"), name));
      new G4PVPlacement(
          0, G4ThreeVector(fFiberX.at(fiberId), fFiberY.at(fiberId), 0),
          fiberCladIntersection__[i].at(fiberId), name, ModuleLogical__[i],
          false, fiberId, checkOverlaps);

      tfiberCoreIntersection = new G4IntersectionSolid(
          "fiberCore", fiberCoreS, module, 0,
          G4ThreeVector(-fFiberX.at(fiberId), -fFiberY.at(fiberId), 0.));
      fiberCoreIntersection__[i].push_back(new G4LogicalVolume(
          tfiberCoreIntersection, FindMaterial("Polystyrene"), name));
      new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
                        fiberCoreIntersection__[i].at(fiberId), name,
                        fiberCladIntersection__[i].at(fiberId), false, fiberId,
                        checkOverlaps);

      fiberCladIntersection__[i].at(fiberId)->SetVisAttributes(fVisAttrGray);
      fiberCoreIntersection__[i].at(fiberId)->SetVisAttributes(fVisAttrOrange);

      if (doGlue) {
        // Create Glue Intersection
        tGlueIntersection = new G4IntersectionSolid(
            "glue", tGlueSubtraction, module, 0,
            G4ThreeVector(-fFiberX.at(fiberId), -fFiberY.at(fiberId), 0.));
        glueIntersection__[i].push_back(new G4LogicalVolume(
            tGlueIntersection, FindMaterial("G4_Galactic"),
            std::string(name) + "_glue_" + std::to_string(fiberId)));
        new G4PVPlacement(
            0, G4ThreeVector(fFiberX.at(fiberId), fFiberY.at(fiberId), 0),
            glueIntersection__[i].at(fiberId),
            std::string(name) + "_glue_" + std::to_string(fiberId),
            ModuleLogical__[i], false, fiberId, checkOverlaps);

        // Set Glue Visualization Attributes
        glueIntersection__[i].at(fiberId)->SetVisAttributes(fVisAttrBlue);
      }
    }
  }
}
