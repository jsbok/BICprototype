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
#include "G4MultiUnion.hh"
#include "G4Transform3D.hh"
#include "G4Colour.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <string>

using namespace std;

G4ThreadLocal koBICMagneticField *koBICDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager *koBICDetectorConstruction::fFieldMgr = 0;

int koBICDetectorConstruction::fNofRow = 1;
int koBICDetectorConstruction::fNofCol = 1;
int koBICDetectorConstruction::fNofModules = fNofRow * fNofCol;

koBICDetectorConstruction::koBICDetectorConstruction()
    : G4VUserDetectorConstruction(), fMessenger(0), fMaterials(NULL) {
  DefineCommands();
  DefineMaterials();

  SiPMT = 0.3 * mm;

  fVisAttrOrange = new G4VisAttributes(G4Colour(1.0, 0.5, 0., 0.7));
  fVisAttrOrange->SetForceSolid(true);
  fVisAttrOrange->SetVisibility(true);
  fVisAttrBlue = new G4VisAttributes(G4Colour(0., 0., 1.0, 0.7));
  fVisAttrBlue->SetForceSolid(false);
  fVisAttrBlue->SetVisibility(false);
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
  fTowerDepth = 320.; 
  G4double fDepth = 710.;
  
  fFiber_vert_dis = 1.2826; // spacing_z
  fFiber_hori_dis = 1.35; // spacing_x
  
  G4double fGlue_thickness = 0.380;
  
  // 5 * (17.0 + 21.730) + 1 * 17.0 + 6 * 20.889 = 335.984 mm
  fModuleH = 30; // = 279.413 mm
  fModuleW = 30;
  G4double fWidth = 900.0; // 700mm + spare space
  fFiberUnitH = 1.;

  G4double rmin = 825.805 * mm;
  G4double rmax = 1036.455 * mm + (30.0 * 6) * mm + (38.73 * 3) * mm ; 
  G4double totalLength = rmax - rmin;

  G4double dPhi = (360. / 48.) * deg;
  G4double dx1 = rmin * std::tan(dPhi/2.); // Inner x
  G4double dx2 = rmax * std::tan(dPhi/2.); // Outer x
  G4double dy  = fWidth / 2.0;           // Width

  // G4Trd (name, x_half_small, x_half_large, y_half_small, y_half_large, z_half)
  G4VSolid* moduleSolid = new G4Trd("moduleSolid", dx1, dx2, dy, dy, totalLength/2.);

  // Vacant Space : AIR or Galactic  
  Envelope = new G4LogicalVolume(moduleSolid, FindMaterial("G4_Galactic"), "ModuleLogical");
  new G4PVPlacement(0, G4ThreeVector(0,0,0), Envelope, "ModulePhys", worldLogical, false, 0, checkOverlaps);

  doFiber = true;
  doSiPM = true;
  doGlue = true;

  // Fiber Dimension
  clad_S_rMax = 0.50 * mm; // EcalBarrel_FiberRadius
  core_S_rMax = 0.48 * mm; // 0.50 - 0.02 (CladdingThickness)
  clad_S_rMax2 = 0.553 * mm;
  glue_S_rMax = 0.539 * mm;
  fiberUnit = new G4Box("fiber_SQ", (fFiberUnitH / 2) * mm, (1. / 2) * mm, (fDepth / 2) * mm);
  fiberClad_SFIL = new G4Tubs("fiberClad_SFIL", 0, clad_S_rMax, 700./ 2., 0 * deg, 360. * deg);
  fiberCore_SFIL = new G4Tubs("fiberCore_SFIL", 0, core_S_rMax, 700. / 2., 0 * deg, 360. * deg);
  fiberGlue_SFIL = new G4Tubs("fiberGlue_SFIL", 0, glue_S_rMax, 700./ 2., 0*deg, 360*deg);

  fiberUnit = new G4Box("fiber_SQ", (fFiberUnitH / 2) * mm, (1. / 2) * mm, (fTowerDepth / 2) * mm);
  fiberClad = new G4Tubs("fiber", 0, clad_S_rMax, fTowerDepth  / 2., 0 * deg, 360. * deg);
  fiberClad2 = new G4Tubs("fiber2", 0, clad_S_rMax2, fTowerDepth  / 2., 0 * deg, 360. * deg);
  fiberCoreS = new G4Tubs("fiberS", 0, core_S_rMax, fTowerDepth  / 2., 0 * deg, 360. * deg);
  gluebox = new G4Box("gluebox", (fFiber_hori_dis / 2) -0.0001 * mm, (fGlue_thickness / 2) - 0.0001 * mm, (fTowerDepth / 2.) - 0.01 * mm );
  tGlueSubtraction = new G4SubtractionSolid("glueCladSubt", gluebox, fiberClad2, 0, G4ThreeVector(.0, .0, .0));

  fCerenRegion = new G4Region("cerenRegion");
  fScintRegion = new G4Region("scintRegion");
  fCerenRegion_SFIL = new G4Region("cerenRegion_SFIL");
  fScintRegion_SFIL = new G4Region("scintRegion_SFIL");

  dimCalc = new dimensionCalc();
  dimCalc->SetFrontL(fFrontL);
  dimCalc->SetTower_height(fTowerDepth);
  dimCalc->SetSiPMT(SiPMT);
  dimCalc->SetNofModules(1);
  dimCalc->SetNofRow(1);
  dimCalc->SetNofCol(1);
  dimCalc->SetModule_height(fModuleH);
  dimCalc->SetModule_width(fModuleW);

  ModuleBuild(ModuleLogical, SiPMGLogical, SiPMfilterLogical, SiPMcellLogical, SiPMcathLogical, fiberUnitIntersection, fiberCladIntersection, fiberClad2Intersection, fiberCoreIntersection, fModuleProp);

  delete dimCalc;
  return worldPhysical;

}

void koBICDetectorConstruction::ConstructSDandField() {
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();

if (doSiPM) {
        // SiPM Definition [0 - 43]
        for (int i = 0; i < 43; i++) {
            G4String nameL = "ModuleSD_L_" + std::to_string(i);
            G4String nameR = "ModuleSD_R_" + std::to_string(i);
            G4String collL = "ModuleC" + std::to_string(2 * i);
            G4String collR = "ModuleC" + std::to_string(2 * i + 1);

            int idL = 2 * i;
            int idR = 2 * i + 1;

            if (i >= (int)fModuleProp.size()) {
                G4cout << "[Warning] fModuleProp size is smaller than 43! Index: " << i << G4endl;
                break; 
            }

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

            // Left Assignment
            if (targetLV_L) {
                targetLV_L->SetSensitiveDetector(SiPMSD_L);
                if (targetCell_L) targetCell_L->SetSensitiveDetector(SiPMSD_L);
            }
            // Right Assignment
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
 
  FastOpTransportModel* cerenModel_SFIL = new FastOpTransportModel("fastOpTransportCeren_SFIL", fCerenRegion_SFIL); 
  FastOpTransportModel* scintModel_SFIL = new FastOpTransportModel("fastOpTransportScint_SFIL", fScintRegion_SFIL);
  cerenModel_SFIL->SetFiberLength(700.0 * mm); 
  cerenModel_SFIL->SetCoreMaterial(FindMaterial("PMMA"));
  scintModel_SFIL->SetFiberLength(700.0 * mm); 
  scintModel_SFIL->SetCoreMaterial(FindMaterial("Polystyrene"));
}

void koBICDetectorConstruction::ModuleBuild(
    G4LogicalVolume *ModuleLogical_[], 
    G4LogicalVolume *SiPMGLogical_[],
    G4LogicalVolume *SiPMfilterLogical_[], 
    G4LogicalVolume *SiPMcellLogical_[],
    G4LogicalVolume *SiPMcathLogical_[],
    std::vector<G4LogicalVolume *> fiberUnitIntersection_[],
    std::vector<G4LogicalVolume *> fiberCladIntersection_[],
    std::vector<G4LogicalVolume *> fiberClad2Intersection_[],
    std::vector<G4LogicalVolume *> fiberCoreIntersection_[],
    std::vector<koBICInterface::koBICModuleProperty> &ModuleProp_) {
    
    G4Material* pbMat = FindMaterial("Lead");
    G4Material* vacMat = FindMaterial("G4_Galactic");

    // Trapezoid dimensions (dx1, dx2, dy, dz)
    G4Trd* parentTrd = (G4Trd*)Envelope->GetSolid();
    G4double pDx1 = parentTrd->GetXHalfLength1();
    G4double pDx2 = parentTrd->GetXHalfLength2();
    G4double pDy  = parentTrd->GetYHalfLength1() - 5;
    G4double pDz  = parentTrd->GetZHalfLength(); // 500mm

    G4double currentZ = -pDz; // -500mm
    int copyNo = 0;

    auto PlaceTrdLayerAt = [&](G4String name, G4double thick, G4double zCenter, G4Material* mat, G4VisAttributes* vis, G4double zForSize = 9999.0) -> G4LogicalVolume* {
        
        G4double zCalc = (zForSize == 9999.0) ? zCenter : zForSize;

        G4double zStart = (zCalc - thick / 2.0) + pDz; 
        G4double zEnd   = zStart + thick;

        G4double layerDx1 = pDx1 + (pDx2 - pDx1) * (zStart / (2 * pDz));
        G4double layerDx2 = pDx1 + (pDx2 - pDx1) * (zEnd / (2 * pDz));

        G4double currentDy = 350.0 * mm;
        
        G4VSolid* layerSolid = new G4Trd(name + "_sol", layerDx1, layerDx2, currentDy, currentDy, thick / 2.0);
        G4LogicalVolume* logVol = new G4LogicalVolume(layerSolid, mat, name + "_log");
        
        if(vis) logVol->SetVisAttributes(vis);

        new G4PVPlacement(0, G4ThreeVector(0, 0, zCenter), 
                          logVol, name + "_phys", Envelope, false, copyNo++, checkOverlaps);
        
        return logVol;
    };

    std::vector<G4LogicalVolume*> logicPbSFIL; //  SFIL(1~5)

    // ==============================================================================
    // 1. [SFIL 1, 2]  Reverse Direction Arrangement (Vacuum 17.49mm -> SFIL 21.24mm)
    // ==============================================================================
    G4double currentBoundaryZ = pDz - (30.0 * 6.0 * mm) - (38.73 * 3.0 * mm);

    for(int r = 0; r < 2; r++) {
        G4double vacThick = 17.49 * mm;
        currentBoundaryZ -= vacThick; 
        G4double zCenterVac = currentBoundaryZ + (vacThick / 2.0); 
        
        PlaceTrdLayerAt("Vac_SFIL", vacThick, zCenterVac, vacMat, fVisAttrBlue);
        
        G4double sfilThick = 21.24 * mm;
        currentBoundaryZ -= sfilThick; 
        G4double zCenterSfil = currentBoundaryZ + (sfilThick / 2.0);
        
        fVisAttrGray->SetColour(G4Colour(0.5, 0.5, 0.5, 0.5)); 
        fVisAttrGray->SetForceSolid(true);
        
        G4LogicalVolume* tmpPbSFIL = PlaceTrdLayerAt("Pb_SFIL", sfilThick, zCenterSfil, pbMat, fVisAttrGray);
        logicPbSFIL.push_back(tmpPbSFIL); 
    }

    // ==============================================================================
    // 2. [SFIL 3, 4, 5] (Bulk) Forward Direction Arrangement, 17mm Interval
    // ==============================================================================
    G4double sfilThick2 = 21.24 * mm; // thickness
    G4double gapThick = 17.49 * mm;     // interval

    // Imaginary Z coordinate for calculating X Distance
    G4double originalBoundaryZ = pDz - (30.0 * 6.0 * mm) - (38.73 * 3.0 * mm);

    // Real Z coordinate
    G4double realBoundaryZ = originalBoundaryZ; 

    for(int r = 0; r < 3; r++) {
        fVisAttrGray->SetVisibility(true);
        fVisAttrGray->SetForceSolid(true);

        G4double originalZCenter = originalBoundaryZ + (sfilThick2 / 2.0); // 1. Imaginary Position

        G4double realZCenter = realBoundaryZ + (sfilThick2 / 2.0); // 2. Real Position

        G4LogicalVolume* tmpPbSFIL2 = PlaceTrdLayerAt(
            "Pb_SFIL", sfilThick2, realZCenter, pbMat, fVisAttrGray, originalZCenter
        );
        logicPbSFIL.push_back(tmpPbSFIL2); // 3. Push SFIL 3, 4, 5 to vector

        // 4. Updates for next layer
        originalBoundaryZ += sfilThick2; 
        
        realBoundaryZ += sfilThick2;

        if (r < 2) {
            G4double zCenterGap = realBoundaryZ + (gapThick / 2.0);
            PlaceTrdLayerAt("Vac_Gap", gapThick, zCenterGap, vacMat, fVisAttrBlue);
            realBoundaryZ += gapThick; // 17 mm gap except last layer
        }
    }

    
    
    G4RotationMatrix* zRot = new G4RotationMatrix;
    zRot->rotateY(M_PI/2.*rad);
    zRot->rotateZ(M_PI/2.*rad);
    
    for (int i = 0; i < fNofModules; i++) {
    moduleName = setModuleName(i);

    dimCalc->SetisModule(true);
    module = new G4Box("Module", (fTowerDepth / 2.) * mm, (fModuleH / 2.) * mm, (fModuleW / 2.) * mm);
    ModuleLogical_[i] = new G4LogicalVolume(module, FindMaterial("Lead"), moduleName);
        
        int quotient  = i / 3; // Z-Axis(0, 0, 0, 1, 1, 1, ..., 5, 5, 5) -> 6 Coloums
        int remainder = i % 3; // X-Axis (0, 1, 2, 0, 1, 2, ..., 0, 1, 2) -> 3 Rows

        // 1. X: +30, 0, -30 mm
        G4double xPos = 30.0 * mm - (30.0 * remainder * mm);
        
        // 2. Y: 0
        G4double yPos = 0.0 * mm;
        
        // 3. Z: quotient * (0~5)
        G4double zCenter = pDz - (30.0 * 6.0 * mm) + (30.0 * quotient * mm) + 15.0 * mm;
        
        new G4PVPlacement(zRot, G4ThreeVector(xPos, yPos, zCenter), ModuleLogical_[i], moduleName,
                          Envelope, false, 100+i, checkOverlaps);

    // Fiber arrangement
    FiberImplement(i, ModuleLogical_, fiberUnitIntersection_,
                   fiberCladIntersection_,fiberClad2Intersection_, fiberCoreIntersection_, logicPbSFIL);
    }

    // =========================================================
    // 3. SiPM
    // =========================================================
if (doSiPM) {
    G4double sipmSize = 13.0 * mm;
    G4double sipmThick = SiPMT;      // 0.3mm
    G4double cathThick = SiPMT;      

    G4VSolid* sipmCellSolid = new G4Box("SiPMCellSolid", sipmSize/2., sipmThick/2., sipmSize/2.);
    G4VSolid* sipmCathSolid = new G4Box("SiPMCathSolid", sipmSize/2., cathThick/2., sipmSize/2.);

    int pbLayerIdx = 0; // SFIL 1~5 Layer Index (0 ~ 4)

    // SiPM Arrangemnet Function
    auto PlaceSiPMsOnLayer = [&](G4double realZ, G4double origZ) {
        // X Distance according to Imaginary Z (origZ)
        G4double origDx = pDx1 + (pDx2 - pDx1) * ((origZ + pDz) / (2.0 * pDz));
        G4double stepX = (origDx * 2.0) / 5.0;

        for (int xIdx = 0; xIdx < 5; xIdx++) {
            int s = pbLayerIdx * 5 + xIdx;

            G4String nameCellL = "ModuleC_Cell" + std::to_string(2 * s);
            G4String nameCellR = "ModuleC_Cell" + std::to_string(2 * s + 1);
            G4String nameCathL = "ModuleC_Cath" + std::to_string(2 * s);
            G4String nameCathR = "ModuleC_Cath" + std::to_string(2 * s + 1);

            // 1. Glass Cell
            SiPMcellLogical_[2 * s]     = new G4LogicalVolume(sipmCellSolid, FindMaterial("Glass"), nameCellL);
            SiPMcellLogical_[2 * s + 1] = new G4LogicalVolume(sipmCellSolid, FindMaterial("Glass"), nameCellR);
            
            // 2. Silicon Cathode 
            SiPMcathLogical_[2 * s]     = new G4LogicalVolume(sipmCathSolid, FindMaterial("Silicon"), nameCathL);
            SiPMcathLogical_[2 * s + 1] = new G4LogicalVolume(sipmCathSolid, FindMaterial("Silicon"), nameCathR);

            new G4PVPlacement(0, G4ThreeVector(), SiPMcathLogical_[2 * s],     nameCathL + "_phys", SiPMcellLogical_[2 * s],     false, 2 * s,     checkOverlaps);
            new G4PVPlacement(0, G4ThreeVector(), SiPMcathLogical_[2 * s + 1], nameCathR + "_phys", SiPMcellLogical_[2 * s + 1], false, 2 * s + 1, checkOverlaps);

            // 3. Glass Cell Position (X : origDx, Z : realZ)
            G4double xPos = -origDx + (xIdx + 0.5) * stepX;
            G4double yPos_L = -(350.0 * mm + sipmThick / 2.0); // 70cm
            G4double yPos_R =  (350.0 * mm + sipmThick / 2.0);

            new G4PVPlacement(0, G4ThreeVector(xPos, yPos_L, realZ), SiPMcellLogical_[2 * s],     nameCellL + "_phys", Envelope, false, 2 * s,     checkOverlaps);
            new G4PVPlacement(0, G4ThreeVector(xPos, yPos_R, realZ), SiPMcellLogical_[2 * s + 1], nameCellR + "_phys", Envelope, false, 2 * s + 1, checkOverlaps);

            // 4. Optical Surface & VisAttributes
            new G4LogicalSkinSurface(nameCathL + "_surf", SiPMcathLogical_[2 * s],     FindSurface("SiPMSurf"));
            new G4LogicalSkinSurface(nameCathR + "_surf", SiPMcathLogical_[2 * s + 1], FindSurface("SiPMSurf"));

            SiPMcathLogical_[2 * s]->SetVisAttributes(fVisAttrGreen);
            SiPMcathLogical_[2 * s + 1]->SetVisAttributes(fVisAttrGreen);
        }
        pbLayerIdx++;
    };

    // ==============================================================================
    // 1. [SFIL 1, 2] SiPM Arrangement (Reverse Direction)
    // ==============================================================================
    G4double currentBoundaryZ = pDz - (30.0 * 6.0 * mm) - (38.73 * 3.0 * mm);

    for (int r = 0; r < 2; r++) {
        currentBoundaryZ -= 17.49 * mm; // Vac
        G4double sfilThick = 21.24 * mm;
        currentBoundaryZ -= sfilThick;
        G4double zCenterSfil = currentBoundaryZ + (sfilThick / 2.0);

        // SFIL 1, 2 : (realZ == origZ)
        PlaceSiPMsOnLayer(zCenterSfil, zCenterSfil);
    }

    // ==============================================================================
    // 2. [SFIL 3, 4, 5] SiPM Arrangement (Forward Direction, 17mm Gap)
    // ==============================================================================
    G4double sfilThick2 = 21.24 * mm;
    G4double gapThick = 17.49 * mm;

    G4double originalBoundaryZ = pDz - (30.0 * 6.0 * mm) - (38.73 * 3.0 * mm);
    G4double realBoundaryZ = originalBoundaryZ;

    for (int r = 0; r < 3; r++) {
        G4double originalZCenter = originalBoundaryZ + (sfilThick2 / 2.0);
        G4double realZCenter = realBoundaryZ + (sfilThick2 / 2.0);

        // SFIL 3, 4, 5 -> Z: realZCenter, X: originalZCenter
        PlaceSiPMsOnLayer(realZCenter, originalZCenter);

        originalBoundaryZ += sfilThick2;
        realBoundaryZ += sfilThick2;

        if (r < 2) {
            realBoundaryZ += gapThick; // 17.49 mm interval
        }
    }

 G4double bulkBoundaryZ = pDz - (30.0 * 6.0 * mm) ; // Bulk 

    // X (+30, 0, -30)
    std::vector<G4double> xPatterns = {30.0 * mm, 0.0 * mm, -30.0 * mm};
    
    // Z (repeat 6)
    int zRepeat = 6;
    G4double zGap = 30.0 * mm; // Z interval

    // 3x6 Z arrangement
    G4double startZ = bulkBoundaryZ +  zGap / 2.0;

    int arrayIdx = 0;
    for (int zIdx = 0; zIdx < zRepeat; zIdx++) {          // 6 Z
        for (int xIdx = 0; xIdx < 3; xIdx++) {          // 3 X (+30, 0, -30)
            
            // ID Assignment (25 ~ 42)
            int s = 25 + arrayIdx; 

            G4String nameCellL = "ModuleC_Cell" + std::to_string(2 * s);
            G4String nameCellR = "ModuleC_Cell" + std::to_string(2 * s + 1);
            G4String nameCathL = "ModuleC_Cath" + std::to_string(2 * s);
            G4String nameCathR = "ModuleC_Cath" + std::to_string(2 * s + 1);

            SiPMcellLogical_[2 * s] = new G4LogicalVolume(sipmCellSolid, FindMaterial("Glass"), nameCellL);
            SiPMcellLogical_[2 * s + 1] = new G4LogicalVolume(sipmCellSolid, FindMaterial("Glass"), nameCellR);
            SiPMcathLogical_[2 * s] = new G4LogicalVolume(sipmCathSolid, FindMaterial("Silicon"), nameCathL);
            SiPMcathLogical_[2 * s + 1] = new G4LogicalVolume(sipmCathSolid, FindMaterial("Silicon"), nameCathR);

            new G4PVPlacement(0, G4ThreeVector(), SiPMcathLogical_[2 * s], nameCathL + "_phys", SiPMcellLogical_[2 * s], false, 2 * s, checkOverlaps);
            new G4PVPlacement(0, G4ThreeVector(), SiPMcathLogical_[2 * s + 1], nameCathR + "_phys", SiPMcellLogical_[2 * s + 1], false, 2 * s + 1, checkOverlaps);

            // 
            G4double xPos = xPatterns[xIdx];
            G4double yPos_L = -(160*mm + sipmThick / 2.); // 
            G4double yPos_R =  (160*mm + sipmThick / 2.); // 
            
            G4double zPos = startZ + (zIdx * zGap);

            new G4PVPlacement(0, G4ThreeVector(xPos, yPos_L, zPos), SiPMcellLogical_[2 * s], nameCellL + "_phys", Envelope, false, 2 * s, checkOverlaps);
            new G4PVPlacement(0, G4ThreeVector(xPos, yPos_R, zPos), SiPMcellLogical_[2 * s + 1], nameCellR + "_phys", Envelope, false, 2 * s + 1, checkOverlaps);

            new G4LogicalSkinSurface(nameCathL + "_surf", SiPMcathLogical_[2 * s], FindSurface("SiPMSurf"));
            new G4LogicalSkinSurface(nameCathR + "_surf", SiPMcathLogical_[2 * s + 1], FindSurface("SiPMSurf"));

            SiPMcathLogical_[2 * s]->SetVisAttributes(fVisAttrGreen);
            SiPMcathLogical_[2 * s + 1]->SetVisAttributes(fVisAttrGreen);

            arrayIdx++;
        }
    }
}


for (int s = 0; s < 43; s++) {
    koBICInterface::koBICModuleProperty ModulePropSingle;
    ModulePropSingle.towerXY = fTowerXY;
    ModulePropSingle.ModuleNum = s; // 0~42
    ModuleProp_.push_back(ModulePropSingle);
    }
}

void koBICDetectorConstruction::DefineCommands() {}
void koBICDetectorConstruction::FiberImplement(
    G4int i, G4LogicalVolume *ModuleLogical__[],
    std::vector<G4LogicalVolume *> fiberUnitIntersection__[],
    std::vector<G4LogicalVolume *> fiberCladIntersection__[],
    std::vector<G4LogicalVolume *> fiberClad2Intersection__[],
    std::vector<G4LogicalVolume *> fiberCoreIntersection__[],
    std::vector<G4LogicalVolume*> logicPbSFIL) {

  fFiberX.clear();
  fFiberY.clear();
  fFiberWhich.clear();

  int NofPlate = fModuleH / (fFiber_vert_dis) ;
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

    if (!doFiber) return;

    fVisAttrSkyBlue = new G4VisAttributes(G4Colour(0.5, 0.8, 0.9, 0.4)); 
    fVisAttrSkyBlue->SetVisibility(true);
    fVisAttrSkyBlue->SetForceSolid(true);

    G4RotationMatrix* fiberRot = new G4RotationMatrix();
    fiberRot->rotateY(90 * deg); 
    G4RotationMatrix* fiberRot2 = new G4RotationMatrix();
    fiberRot2->rotateX(90 * deg); 

    G4double pitch_horizontal = 1.35 * mm; 
    G4double pitch_vertical = 1.22 * mm;   
    int fiberId = 0;

  if (doFiber) {
    for (unsigned int fiberId = 0; fiberId < fFiberX.size(); fiberId++) {

      fiberClad2Intersection__[i].push_back(new G4LogicalVolume(
          fiberClad2, FindMaterial("G4_Galactic"), name));
      new G4PVPlacement(
          fiberRot, G4ThreeVector(0, fFiberX.at(fiberId), fFiberY.at(fiberId)), 
          fiberClad2Intersection__[i].at(fiberId), name, ModuleLogical__[i],
          false, fiberId, checkOverlaps);

      fiberCladIntersection__[i].push_back(new G4LogicalVolume(
          fiberClad, FindMaterial("PMMA"), name));
      new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
                        fiberCladIntersection__[i].at(fiberId), name,
                        fiberClad2Intersection__[i].at(fiberId), false, fiberId,
                        checkOverlaps);

      fiberCoreIntersection__[i].push_back(new G4LogicalVolume(
          fiberCoreS, FindMaterial("Polystyrene"), name));
      new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
                        fiberCoreIntersection__[i].at(fiberId), name,
                        fiberCladIntersection__[i].at(fiberId), false, fiberId,
                        checkOverlaps);

      fiberClad2Intersection__[i].at(fiberId)->SetVisAttributes(fVisAttrSkyBlue);
      fiberCladIntersection__[i].at(fiberId)->SetVisAttributes(fVisAttrGray);
      fiberCoreIntersection__[i].at(fiberId)->SetVisAttributes(fVisAttrOrange);

      fScintRegion->AddRootLogicalVolume(fiberCladIntersection__[i].at(fiberId));
      fScintRegion->AddRootLogicalVolume(fiberCoreIntersection__[i].at(fiberId));
      fiberCladIntersection__[i].at(fiberId)->SetRegion(fScintRegion);
      fiberCoreIntersection__[i].at(fiberId)->SetRegion(fScintRegion);

//    G4VSolid*module2 = new G4Box("Module2", (fTowerDepth / 2.) -1 * mm, (fModuleH / 2.) -1 * mm, (fModuleW / 2.) -1 * mm);

      if (doGlue) {
        // 1. Glue Intersection 
        G4RotationMatrix* glueRotB = new G4RotationMatrix(fiberRot->inverse());
        G4ThreeVector glueTransB(0, -fFiberX.at(fiberId), -fFiberY.at(fiberId));
        glueTransB.transform(*glueRotB); // 

        tGlueIntersection = new G4IntersectionSolid(
            "glue", tGlueSubtraction, module, glueRotB, glueTransB);
            
        glueIntersection__[i].push_back(new G4LogicalVolume(
            tGlueIntersection, FindMaterial("G4_Galactic"),
            std::string(name) + "_glue_" + std::to_string(fiberId)));
            
        new G4PVPlacement(
            fiberRot, G4ThreeVector(0, fFiberX.at(fiberId), fFiberY.at(fiberId)),
            glueIntersection__[i].at(fiberId),
            std::string(name) + "_glue_" + std::to_string(fiberId),
            ModuleLogical__[i], false, fiberId, checkOverlaps);
       }


   }
        G4float y1 = fModuleW * mm / 2;
        G4float y2 = - fModuleW * mm / 2 + 22 * fFiber_hori_dis * mm  ;
        G4float y3 = 14.85; G4float y4 = 14.1125; G4float y5 = -14.4125;

        for (int k2=1; k2 < 12; k2++  )  {
          G4float x1 = -fModuleH * mm / 2 + k2 * 2 * fFiber_vert_dis * mm -
                       fFiber_vert_dis / 2 * mm;

          // --- Air1 Intersection & Placement ---
          G4RotationMatrix* airRotB1 = new G4RotationMatrix(fiberRot->inverse());
          G4ThreeVector airTransB1(0, -x1, -y1); 
          airTransB1.transform(*airRotB1);

          G4IntersectionSolid* tAirIntersection1 = new G4IntersectionSolid(
              "Air1", fiberClad, module, airRotB1, airTransB1);
              
          G4LogicalVolume* logicAir1 = new G4LogicalVolume(
              tAirIntersection1, FindMaterial("G4_Galactic"), "Air"+std::to_string(i)+std::to_string(k2));
              
          new G4PVPlacement(
              fiberRot, G4ThreeVector(0, x1, -y1),
              logicAir1, "AirPhysical"+std::to_string(i)+std::to_string(k2), ModuleLogical__[i],
              false, k2, false);

          // --- Air2 Intersection & Placement ---
          G4RotationMatrix* airRotB2 = new G4RotationMatrix(fiberRot->inverse());
          G4ThreeVector airTransB2(0, -x1, y2);
          airTransB2.transform(*airRotB2);

          G4IntersectionSolid* tAirIntersection2 = new G4IntersectionSolid(
              "Air2", fiberClad, module, airRotB2, airTransB2);
              
          G4LogicalVolume* logicAir2 = new G4LogicalVolume(
              tAirIntersection2, FindMaterial("G4_Galactic"), "Air"+std::to_string(i)+std::to_string(k2));
              
          new G4PVPlacement(
              fiberRot, G4ThreeVector(0, x1, y2),
              logicAir2, "AirPhysical2"+std::to_string(i)+std::to_string(k2), ModuleLogical__[i],
              false, k2, false);

          logicAir1 -> SetVisAttributes(fVisAttrSkyBlue);
          logicAir2 -> SetVisAttributes(fVisAttrSkyBlue);     

          // --- Space3 & Space4 Placement ---
          G4VSolid *space3 = new G4Box("glue4", (0.175 / 2.) - 0.0001 * mm,
                            (0.38 / 2.) - 0.0001 * mm, (320. / 2.) - 0.0001 * mm);
                            
          G4LogicalVolume* spaceLogical3 = new G4LogicalVolume(
              space3, FindMaterial("G4_Galactic"), "spaceLogical3"+std::to_string(i));
              
          new G4PVPlacement(
              fiberRot, G4ThreeVector(0, x1, y4), spaceLogical3, 
              "spacePhysical3"+std::to_string(i), ModuleLogical__[i], false, 0, checkOverlaps);

          spaceLogical3 -> SetVisAttributes(fVisAttrSkyBlue); 
         
          G4LogicalVolume* spaceLogical4 = new G4LogicalVolume(
              space3, FindMaterial("G4_Galactic"), "spaceLogical4"+std::to_string(i));
              
          new G4PVPlacement(
              fiberRot, G4ThreeVector(0, x1, y5), spaceLogical4, 
              "spacePhysical4"+std::to_string(i), ModuleLogical__[i], false, 0, checkOverlaps);
              
          spaceLogical4 -> SetVisAttributes(fVisAttrSkyBlue);
        }

        // space 3
        for (int k2=1; k2 < 13; k2++  )  {   
          G4float x2 = -fModuleH * mm / 2 + k2 * 2 * fFiber_vert_dis * mm -
                       3 * fFiber_vert_dis / 2 * mm;

          G4VSolid *space2 = new G4Box("glue3", (0.3 / 2.) * mm,
                            (0.38 / 2.) * mm, (320. / 2.) - 0.0001 * mm );
                            
          G4LogicalVolume* spaceLogical2 = new G4LogicalVolume(
              space2, FindMaterial("G4_Galactic"), "spaceLogical2"+std::to_string(i));
              
          new G4PVPlacement(
              fiberRot, G4ThreeVector(0, x2, y3), spaceLogical2, 
              "spacePhysical2"+std::to_string(i), ModuleLogical__[i], false, 0, checkOverlaps);

          spaceLogical2 -> SetVisAttributes(fVisAttrSkyBlue);
        }

        glueIntersection__[i].at(fiberId)->SetVisAttributes(fVisAttrSkyBlue);
      
    
  }

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

            G4MultiUnion* combinedFibersSolid = nullptr;
            if (doGlue) {
                combinedFibersSolid = new G4MultiUnion("CombinedFibers");
            }

            for (int j = 0; j < fixedNumCols; j++) {
                G4double localH = startH + j * pitch_horizontal;
                if (isShifted) localH += pitch_horizontal / 2.0;

                if (std::abs(localH) + 0.553 * mm > currentHalfWidth - 0.1 * mm) continue;

                new G4PVPlacement(fiberRot2, G4ThreeVector(localH, 0, localV), targetGlueLog, "fiberGlue_Phys", motherLog, false, fiberId, false);
                fiberId++;

                if (doGlue && combinedFibersSolid != nullptr) {
                    G4RotationMatrix rot = (fiberRot2 ? *fiberRot2 : G4RotationMatrix());
                    G4Transform3D transform(rot, G4ThreeVector(localH, 0, 0));
                    
                    combinedFibersSolid->AddNode(*(targetGlueLog->GetSolid()), transform);
                }
            }

            if (doGlue && combinedFibersSolid != nullptr) {

                combinedFibersSolid->Voxelize();

                G4VSolid* rowGlueBox = new G4Box("RowGlueBox", currentHalfWidth- 0.01 * mm, (700. / 2.0) - 0.001 * mm, 0.233 / 2.0);
                
                G4VSolid* finalGlueSolid = new G4SubtractionSolid("SubtractedGlue", rowGlueBox, combinedFibersSolid, nullptr, G4ThreeVector(0, 0, 0));

                G4LogicalVolume* rowGlueLog = new G4LogicalVolume(finalGlueSolid, FindMaterial("Glue"), "rowGlue_Log");
                rowGlueLog->SetVisAttributes(fVisAttrSkyBlue);
                glueIntersection__[i].push_back(rowGlueLog);

                new G4PVPlacement(nullptr, G4ThreeVector(0, 0, localV), rowGlueLog, "rowGlue_Phys", motherLog, false, k, false);
            }
        }
        
    };

    // [SFIL Layers Loop]
    for (auto pbLog : logicPbSFIL) {
        G4LogicalVolume* glueLog_SFIL = new G4LogicalVolume(fiberGlue_SFIL, FindMaterial("Glue"), "fGlue_SFIL_L");
        G4LogicalVolume* cladLog_SFIL = new G4LogicalVolume(fiberClad_SFIL, FindMaterial("PMMA"), "fClad_SFIL_L");
        G4LogicalVolume* coreLog_SFIL = new G4LogicalVolume(fiberCore_SFIL, FindMaterial("Polystyrene"), "fCore_SFIL_L");

	fScintRegion_SFIL->AddRootLogicalVolume(cladLog_SFIL);
        fScintRegion_SFIL->AddRootLogicalVolume(coreLog_SFIL);
        cladLog_SFIL->SetRegion(fScintRegion_SFIL);
        coreLog_SFIL->SetRegion(fScintRegion_SFIL);

        new G4PVPlacement(0, G4ThreeVector(0,0,0), coreLog_SFIL, "fCore_P", cladLog_SFIL, false, 0, false);
        new G4PVPlacement(0, G4ThreeVector(0,0,0), cladLog_SFIL, "fClad_P", glueLog_SFIL, false, 0, false);

        glueLog_SFIL->SetVisAttributes(fVisAttrSkyBlue);
        cladLog_SFIL->SetVisAttributes(fVisAttrGray);
        coreLog_SFIL->SetVisAttributes(fVisAttrOrange);

        FillFibersInMother(pbLog, 21.24 * mm, glueLog_SFIL, 700.0 * mm, 0.233 * mm);
    }
}
