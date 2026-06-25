#include "koBICDetectorConstruction.hh"
#include "koBICCellParameterisation.hh"
#include "koBICFilterParameterisation.hh"
#include "koBICMirrorParameterisation.hh"
#include "koBICSiPMSD.hh"
#include "G4Trd.hh"
#include "G4VisExtent.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "FastOpTransportModel.hh"
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

G4ThreadLocal koBICMagneticField *koBICDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager *koBICDetectorConstruction::fFieldMgr = 0;

koBICDetectorConstruction::koBICDetectorConstruction()
    : G4VUserDetectorConstruction(), fMessenger(0), fMaterials(NULL) {
  DefineCommands();
  DefineMaterials();

  PMTT = 0.3 * mm;

  fVisAttrOrange = new G4VisAttributes(G4Colour(1.0, 0.5, 0., 0.7));
  fVisAttrOrange->SetForceSolid(true);
  fVisAttrOrange->SetVisibility(true);
  fVisAttrBlue = new G4VisAttributes(G4Colour(0., 0., 1.0, 0.7));
  fVisAttrBlue->SetForceSolid(true);
  fVisAttrBlue->SetVisibility(true);
  fVisAttrGray = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3, 0.2));
  fVisAttrGray->SetVisibility(true);
  fVisAttrGreen = new G4VisAttributes(G4Colour(0.3, 0.7, 0.3, 0.7));
  fVisAttrGreen->SetVisibility(true);
  fVisAttrSkyBlue = new G4VisAttributes(G4Colour(0.5, 0.8, 0.9, 0.7));
  fVisAttrSkyBlue->SetVisibility(true);
}

koBICDetectorConstruction::~koBICDetectorConstruction() {
  delete fMessenger;
  delete fMaterials;

  delete fVisAttrOrange;
  delete fVisAttrBlue;
  delete fVisAttrGray;
  delete fVisAttrGreen;
  delete fVisAttrSkyBlue;
}

void koBICDetectorConstruction::DefineMaterials() {
  fMaterials = koBICMaterials::GetInstance();
}

G4VPhysicalVolume *koBICDetectorConstruction::Construct() {
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  checkOverlaps = false; 

  G4VSolid *worldSolid = new G4Box("worldBox", 10. * m, 10. * m, 10. * m);
  worldLogical = new G4LogicalVolume(worldSolid, FindMaterial("G4_Galactic"), "worldLogical");
  G4VPhysicalVolume *worldPhysical =
      new G4PVPlacement(0, G4ThreeVector(), worldLogical, "worldPhysical", 0, false, 0, checkOverlaps);

  fFrontL = 1000.; 
  fTowerDepth = 700.; 
  
  fFiber_vert_dis = 1.22; // spacing_z
  fFiber_hori_dis = 1.35; // spacing_x
  
  // 5 * (17.0 + 21.730) + 1 * 17.0 + 6 * 20.889 = 335.984 mm
  fModuleH = 133.19 + 146.223; // = 279.413 mm
  fModuleW = 710.0; // 700mm + spare space
  fFiberUnitH = 1.;

  G4double rmin = 904.485 * mm; 
  G4double rmax = 1037.675 * mm + (20.889 * 7) * mm; 
  G4double totalLength = rmax - rmin * mm;

  G4double dPhi = (360. / 48.) * deg;
  G4double dx1 = rmin * std::tan(dPhi/2.); // Inner x
  G4double dx2 = rmax * std::tan(dPhi/2.); // Outer x
  G4double dy  = fModuleW / 2.0;           // Width

  // G4Trd (name, x_half_small, x_half_large, y_half_small, y_half_large, z_half)
  G4VSolid* moduleSolid = new G4Trd("moduleSolid", dx1, dx2, dy, dy, totalLength/2.);

  // Vacant Space : AIR or Galactic  
  ModuleLogical[0] = new G4LogicalVolume(moduleSolid, FindMaterial("G4_Galactic"), "ModuleLogical");
  new G4PVPlacement(0, G4ThreeVector(0,0,0), ModuleLogical[0], "ModulePhys", worldLogical, false, 0, checkOverlaps);

  doFiber = true;
  doPMT = true;
  doGlue = true;

  // Fiber Dimension
  clad_S_rMax = 0.50 * mm; // EcalBarrel_FiberRadius
  core_S_rMax = 0.485 * mm; // 0.50 - 0.04 (CladdingThickness)
  glue_S_rMax = 0.553 * mm;
  fiberUnit = new G4Box("fiber_SQ", (fFiberUnitH / 2) * mm, (1. / 2) * mm, (fTowerDepth / 2) * mm);
  fiberClad_SFIL = new G4Tubs("fiberClad_SFIL", 0, clad_S_rMax, 700./ 2., 0 * deg, 360. * deg);
  fiberCore_SFIL = new G4Tubs("fiberCore_SFIL", 0, core_S_rMax, 700. / 2., 0 * deg, 360. * deg);
  fiberGlue_SFIL = new G4Tubs("fiberGlue_SFIL", 0, glue_S_rMax, 700./ 2., 0*deg, 360*deg);

  fiberClad_Bulk = new G4Tubs("fiberClad_Bulk", 0, clad_S_rMax, 300. / 2., 0 * deg, 360. * deg);
  fiberCore_Bulk = new G4Tubs("fiberCore_Bulk", 0, core_S_rMax, 300. / 2., 0 * deg, 360. * deg);
  fiberGlue_Bulk = new G4Tubs("fiberGlue_Bulk", 0, glue_S_rMax, 300./ 2., 0*deg, 360*deg);

  fCerenRegion = new G4Region("cerenRegion");
  fScintRegion = new G4Region("scintRegion");

  dimCalc = new dimensionCalc();
  dimCalc->SetFrontL(fFrontL);
  dimCalc->SetTower_height(fTowerDepth);
  dimCalc->SetPMTT(PMTT);
  dimCalc->SetNofModules(1);
  dimCalc->SetNofRow(1);
  dimCalc->SetNofCol(1);
  dimCalc->SetModule_height(fModuleH);
  dimCalc->SetModule_width(fModuleW);

  ModuleBuild(ModuleLogical, PMTGLogical, PMTfilterLogical, PMTcellLogical, PMTcathLogical, fModuleProp);

  delete dimCalc;
  return worldPhysical;

}

void koBICDetectorConstruction::ConstructSDandField() {
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();

    if (doPMT) {
        for (int i = 0; i < 50; i++) {
            G4String nameL = "ModuleSD_L_" + std::to_string(i);
            G4String nameR = "ModuleSD_R_" + std::to_string(i);
            G4String collL = "ModuleC" + std::to_string(2 * i);
            G4String collR = "ModuleC" + std::to_string(2 * i + 1);

            int idL = 2 * i;
            int idR = 2 * i + 1;

            koBICSiPMSD* SiPMSD_L = new koBICSiPMSD(nameL, collL, 0, fModuleProp.at(i), idL);
            koBICSiPMSD* SiPMSD_R = new koBICSiPMSD(nameR, collR, 1, fModuleProp.at(i), idR);

            SDman->AddNewDetector(SiPMSD_L);
            SDman->AddNewDetector(SiPMSD_R);

            // Cathode(Silicon) 
            G4String lvNameL = "ModuleC_Cath" + std::to_string(2 * i);
            G4String lvNameR = "ModuleC_Cath" + std::to_string(2 * i + 1);
            G4LogicalVolume* targetLV_L = lvStore->GetVolume(lvNameL);
            G4LogicalVolume* targetLV_R = lvStore->GetVolume(lvNameR);

            // Cell(Glass) 
            G4String lvNameL_Cell = "ModuleC_Cell" + std::to_string(2 * i);
            G4String lvNameR_Cell = "ModuleC_Cell" + std::to_string(2 * i + 1);
            G4LogicalVolume* targetCell_L = lvStore->GetVolume(lvNameL_Cell);
            G4LogicalVolume* targetCell_R = lvStore->GetVolume(lvNameR_Cell);

            // Left Assignmnet
            if (targetLV_L) {
                targetLV_L->SetSensitiveDetector(SiPMSD_L);
                if (targetCell_L) targetCell_L->SetSensitiveDetector(SiPMSD_L);
            }
            // Right Assignmnet
            if (targetLV_R) {
                targetLV_R->SetSensitiveDetector(SiPMSD_R);
                if (targetCell_R) targetCell_R->SetSensitiveDetector(SiPMSD_R);
            }
        }
    }
  FastOpTransportModel* cerenModel = new FastOpTransportModel("fastOpTransportCeren",fCerenRegion);
  FastOpTransportModel* scintModel = new FastOpTransportModel("fastOpTransportScint",fScintRegion);
  cerenModel->SetFiberLength(fTowerDepth);
  cerenModel->SetCoreMaterial(FindMaterial("PMMA"));
  scintModel->SetFiberLength(fTowerDepth);
  scintModel->SetCoreMaterial(FindMaterial("Polystyrene"));
}

void koBICDetectorConstruction::ModuleBuild(
    G4LogicalVolume *ModuleLogical_[], 
    G4LogicalVolume *PMTGLogical_[],
    G4LogicalVolume *PMTfilterLogical_[], 
    G4LogicalVolume *PMTcellLogical_[],
    G4LogicalVolume *PMTcathLogical_[],
    std::vector<koBICInterface::koBICModuleProperty> &ModuleProp_) {
    G4Material* pbMat = FindMaterial("Lead");
    G4Material* vacMat = FindMaterial("G4_Galactic");

    // Trapezoid dimensions (dx1, dx2, dy, dz)
    G4Trd* parentTrd = (G4Trd*)ModuleLogical_[0]->GetSolid();
    G4double pDx1 = parentTrd->GetXHalfLength1();
    G4double pDx2 = parentTrd->GetXHalfLength2();
    G4double pDy  = parentTrd->GetYHalfLength1() - 5;
    G4double pDz  = parentTrd->GetZHalfLength(); // 500mm

    G4double currentZ = -pDz; // -500mm
    int copyNo = 0;

    // Layer production
    auto PlaceTrdLayer = [&](G4String name, G4double thick, G4Material* mat, G4VisAttributes* vis, G4double customDy = -1.0) -> G4LogicalVolume* {
    G4double zStart = currentZ + pDz; 
    G4double zEnd   = zStart + thick;

    G4double layerDx1 = pDx1 + (pDx2 - pDx1) * (zStart / (2 * pDz));
    G4double layerDx2 = pDx1 + (pDx2 - pDx1) * (zEnd / (2 * pDz));

    G4double currentDy = (customDy > 0) ? customDy : pDy;

    G4VSolid* layerSolid = new G4Trd(name + "_sol", layerDx1, layerDx2, currentDy, currentDy, thick / 2.0);
    G4LogicalVolume* logVol = new G4LogicalVolume(layerSolid, mat, name + "_log");
    
    if(vis) logVol->SetVisAttributes(vis);

    new G4PVPlacement(0, G4ThreeVector(0, 0, currentZ + thick / 2.0), 
                      logVol, name + "_phys", ModuleLogical_[0], false, copyNo++, checkOverlaps);
    
    currentZ += thick;
    return logVol;
};


std::vector<G4LogicalVolume*> logicPbSFIL;
std::vector<G4LogicalVolume*> logicPbBulk;

for(int r = 0; r < 3; r++) {
    fVisAttrBlue->SetVisibility(false);
    PlaceTrdLayer("Vac_SFIL", 17.0*mm, vacMat, fVisAttrBlue);
    
    fVisAttrGray->SetColour(G4Colour(0.5, 0.5, 0.5, 0.5)); 
    fVisAttrGray->SetForceSolid(true);
    
    G4LogicalVolume* tmpPbSFIL = PlaceTrdLayer("Pb_SFIL", 21.73*mm, pbMat, fVisAttrGray);
    logicPbSFIL.push_back(tmpPbSFIL); 
}

fVisAttrBlue->SetVisibility(false);
PlaceTrdLayer("Vac_SFIL_End", 17.0*mm, vacMat, fVisAttrBlue);

G4double bulkHalfHeight = (300.0 * mm) / 2.0; 

// --- Bulk  ---
for(int r = 0; r < 7; r++) {
    fVisAttrGray->SetVisibility(true);
    fVisAttrGray->SetForceSolid(true);
    
    G4LogicalVolume* tmpPbBulk = PlaceTrdLayer("Pb_Bulk", 20.889*mm, pbMat, fVisAttrGray, bulkHalfHeight);
    logicPbBulk.push_back(tmpPbBulk);
}

    // Fiber arrangement
    FiberImplement(0, logicPbSFIL, logicPbBulk);

    // =========================================================
    // 3. PMT
    // =========================================================
if (doPMT) {
    G4double sipmSize = 16.0 * mm;
    G4double sipmThick = PMTT;      // 0.3mm
    G4double cathThick = PMTT;      

    std::vector<G4double> layerThicks = {
        17.0*mm, 21.73*mm, 17.0*mm, 21.73*mm, 17.0*mm, 21.73*mm, 17.0*mm,
        20.889*mm, 20.889*mm, 20.889*mm, 20.889*mm, 20.889*mm, 20.889*mm, 20.889*mm
    };

    G4double runningZ = -pDz;
    int pbLayerIdx = 0;


  fVisAttrSkyBlue = new G4VisAttributes(G4Colour(0.5, 0.8, 0.9, 0.5));
  fVisAttrSkyBlue->SetVisibility(true);
  fVisAttrSkyBlue->SetForceSolid(true);

    G4VSolid* sipmCellSolid = new G4Box("SiPMCellSolid", sipmSize/2., sipmThick/2., sipmSize/2.);
    G4VSolid* sipmCathSolid = new G4Box("SiPMCathSolid", sipmSize/2., cathThick/2., sipmSize/2.);

    for (int i = 0; i < 14; i++) {

        G4double thick = layerThicks[i];
        G4double zCenter = runningZ + thick/2.0;
        bool isPbLayer = (i == 1 || i == 3 || i == 5 || i >= 7);

if (isPbLayer) {
        G4double currentDy = (i >= 7) ? (150.0 * mm) : pDy;
        G4double currentDx = pDx1 + (pDx2 - pDx1) * ((zCenter + pDz) / (2.0 * pDz));
        G4double stepX = (currentDx * 2.0) / 5.0;

        for (int xIdx = 0; xIdx < 5; xIdx++) {
            int s = pbLayerIdx * 5 + xIdx;

            // SiPM name
            G4String nameCellL = "ModuleC_Cell" + std::to_string(2 * s);
            G4String nameCellR = "ModuleC_Cell" + std::to_string(2 * s + 1);
            G4String nameCathL = "ModuleC_Cath" + std::to_string(2 * s);
            G4String nameCathR = "ModuleC_Cath" + std::to_string(2 * s + 1);

            // 1. Glass Cell 
            PMTcellLogical_[2 * s] = new G4LogicalVolume(sipmCellSolid, FindMaterial("Glass"), nameCellL);
            PMTcellLogical_[2 * s + 1] = new G4LogicalVolume(sipmCellSolid, FindMaterial("Glass"), nameCellR);
            
            // 2. Silicon Cathode
            PMTcathLogical_[2 * s] = new G4LogicalVolume(sipmCathSolid, FindMaterial("Silicon"), nameCathL);
            PMTcathLogical_[2 * s + 1] = new G4LogicalVolume(sipmCathSolid, FindMaterial("Silicon"), nameCathR);

            new G4PVPlacement(0, G4ThreeVector(), PMTcathLogical_[2 * s], nameCathL + "_phys", PMTcellLogical_[2 * s], false, 2 * s, checkOverlaps);
            new G4PVPlacement(0, G4ThreeVector(), PMTcathLogical_[2 * s + 1], nameCathR + "_phys", PMTcellLogical_[2 * s + 1], false, 2 * s + 1, checkOverlaps);

            // 3. Glass Cell Arrangement
            G4double xPos = -currentDx + (xIdx + 0.5) * stepX;
            G4double yPos_L = -(currentDy + sipmThick / 2.);
            G4double yPos_R =  (currentDy + sipmThick / 2.);

            new G4PVPlacement(0, G4ThreeVector(xPos, yPos_L, zCenter), PMTcellLogical_[2 * s], nameCellL + "_phys", ModuleLogical_[0], false, 2 * s, checkOverlaps);
            new G4PVPlacement(0, G4ThreeVector(xPos, yPos_R, zCenter), PMTcellLogical_[2 * s + 1], nameCellR + "_phys", ModuleLogical_[0], false, 2 * s + 1, checkOverlaps);

            // 4. Optical Surface
            new G4LogicalSkinSurface(nameCathL + "_surf", PMTcathLogical_[2 * s], FindSurface("SiPMSurf"));
            new G4LogicalSkinSurface(nameCathR + "_surf", PMTcathLogical_[2 * s + 1], FindSurface("SiPMSurf"));

            PMTcathLogical_[2 * s]->SetVisAttributes(fVisAttrGreen);
            PMTcathLogical_[2 * s + 1]->SetVisAttributes(fVisAttrGreen);
            
        }
        pbLayerIdx++;
    }
        runningZ += thick;
    }
}

for (int s = 0; s < 50; s++) {
    koBICInterface::koBICModuleProperty ModulePropSingle;
    ModulePropSingle.towerXY = fTowerXY;
    ModulePropSingle.ModuleNum = s; // 0~49
    ModuleProp_.push_back(ModulePropSingle);
    }
}

void koBICDetectorConstruction::DefineCommands() {}
void koBICDetectorConstruction::FiberImplement(
    G4int i, 
    std::vector<G4LogicalVolume*> logicPbSFIL, 
    std::vector<G4LogicalVolume*> logicPbBulk) {

    if (!doFiber) return;

    fVisAttrSkyBlue = new G4VisAttributes(G4Colour(0.5, 0.8, 0.9, 0.4)); 
    fVisAttrSkyBlue->SetVisibility(true);
    fVisAttrSkyBlue->SetForceSolid(true);

    G4RotationMatrix* fiberRot = new G4RotationMatrix();
    fiberRot->rotateX(90 * deg); 

    G4double pitch_horizontal = 1.35 * mm; 
    G4double pitch_vertical = 1.22 * mm;   
    int fiberId = 0;

    auto FillFibersInMother = [&](G4LogicalVolume* motherLog, G4double thickness, G4LogicalVolume* targetGlueLog, G4double zLength, G4double glue_thickness) {

        G4Trd* motherSolid = (G4Trd*)motherLog->GetSolid();
        G4double h1 = motherSolid->GetXHalfLength1();
        G4double h2 = motherSolid->GetXHalfLength2();
        G4double zHalf = motherSolid->GetZHalfLength(); 

        int numRows = std::floor(thickness / pitch_vertical);
        G4double startV = -((numRows - 1) * pitch_vertical) / 2.0;

        G4double maxHalfWidth = std::max(h1, h2);
        int fixedNumCols = std::floor((2.0 * maxHalfWidth - 1.0 * mm) / pitch_horizontal);
        G4double startH = -((fixedNumCols - 1) * pitch_horizontal) / 2.0;

        // Fiber Row Loop
        for (int k = 0; k < numRows; k++) {

            G4double localV = startV + k * pitch_vertical;
            G4double currentHalfWidth = h1 + (h2 - h1) * (localV + zHalf) / (2.0 * zHalf);
            bool isShifted = (k % 2 != 0);

            // -----------------------------------------------------------------
            G4VSolid* rowGlueTrdSolid = new G4Trd("RowGlueTrdBase", zLength / 2.0, zLength / 2.0, glue_thickness / 2.0, glue_thickness / 2.0, zHalf);

            G4IntersectionSolid* trimmedRowGlueSolid = new G4IntersectionSolid(
                "TrimmedRowGlue", rowGlueTrdSolid, motherSolid, nullptr, G4ThreeVector(0., 0., -localV)
            );

	    if (doGlue) {
                // - (X): Row Width (currentHalfWidth)
                // - (Y): Glue Thickness (glue_thickness / 2.0)
                // - (Z): Fiber Length (zLength / 2.0) 
                G4VSolid* rowGlueSolid = new G4Box("RowGlueBox", currentHalfWidth, glue_thickness / 2.0, zLength / 2.0);

              /*
                G4double sizeX = currentHalfWidth * 2.0;   //  width
                G4double sizeY = glue_thickness;          // 0.233 mm
                G4double sizeZ = zLength;                 // 300 mm or 700 mm (Fiber Length)

                G4cout << "[Layer " << i << " - Row " << k << "] Size X: " << sizeX/mm << " mm | Size Y: " << sizeY/mm << " mm | Size Z: " << sizeZ/mm << " mm" << G4endl;
*/
                G4LogicalVolume* rowGlueLog = new G4LogicalVolume(rowGlueSolid, FindMaterial("Glue"), "rowGlue_Log");
                rowGlueLog->SetVisAttributes(fVisAttrSkyBlue);
                glueIntersection__[i].push_back(rowGlueLog);

                new G4PVPlacement(fiberRot, G4ThreeVector(0, 0, localV), rowGlueLog, "rowGlue_Phys", motherLog, false, k, false);
            }
            // -----------------------------------------------------------------

            // Fiber Arrangement Loop
            for (int j = 0; j < fixedNumCols; j++) {
                G4double localH = startH + j * pitch_horizontal;
                if (isShifted) localH += pitch_horizontal / 2.0;

                if (std::abs(localH) + 0.553 * mm > currentHalfWidth - 0.1 * mm) continue;

                new G4PVPlacement(fiberRot, G4ThreeVector(localH, 0, localV), targetGlueLog, "fiberGlue_Phys", motherLog, false, fiberId, false);
                fiberId++;
            }
        }
    };

    // [SFIL Layers Loop]
    for (auto pbLog : logicPbSFIL) {
        G4LogicalVolume* glueLog_SFIL = new G4LogicalVolume(fiberGlue_SFIL, FindMaterial("Glue"), "fGlue_SFIL_L");
        G4LogicalVolume* cladLog_SFIL = new G4LogicalVolume(fiberClad_SFIL, FindMaterial("PMMA"), "fClad_SFIL_L");
        G4LogicalVolume* coreLog_SFIL = new G4LogicalVolume(fiberCore_SFIL, FindMaterial("Polystyrene"), "fCore_SFIL_L");

        new G4PVPlacement(0, G4ThreeVector(0,0,0), coreLog_SFIL, "fCore_P", cladLog_SFIL, false, 0, false);
        new G4PVPlacement(0, G4ThreeVector(0,0,0), cladLog_SFIL, "fClad_P", glueLog_SFIL, false, 0, false);

        glueLog_SFIL->SetVisAttributes(fVisAttrSkyBlue);
        cladLog_SFIL->SetVisAttributes(fVisAttrGray);
        coreLog_SFIL->SetVisAttributes(fVisAttrOrange);

        FillFibersInMother(pbLog, 21.73 * mm, glueLog_SFIL, 700.0 * mm, 0.233 * mm);
    }

    // [Bulk Layers Loop]
    for (auto pbLog : logicPbBulk) {
        G4LogicalVolume* glueLog_Bulk = new G4LogicalVolume(fiberGlue_Bulk, FindMaterial("Glue"), "fGlue_Bulk_L");
        G4LogicalVolume* cladLog_Bulk = new G4LogicalVolume(fiberClad_Bulk, FindMaterial("PMMA"), "fClad_Bulk_L");
        G4LogicalVolume* coreLog_Bulk = new G4LogicalVolume(fiberCore_Bulk, FindMaterial("Polystyrene"), "fCore_Bulk_L");

        new G4PVPlacement(0, G4ThreeVector(0,0,0), coreLog_Bulk, "fCore_P", cladLog_Bulk, false, 0, false);
        new G4PVPlacement(0, G4ThreeVector(0,0,0), cladLog_Bulk, "fClad_P", glueLog_Bulk, false, 0, false);

        glueLog_Bulk->SetVisAttributes(fVisAttrSkyBlue);
        cladLog_Bulk->SetVisAttributes(fVisAttrGray);
        coreLog_Bulk->SetVisAttributes(fVisAttrOrange);

        FillFibersInMother(pbLog, 20.889 * mm, glueLog_Bulk, 300.0 * mm, 0.233 * mm); 
    }
}
