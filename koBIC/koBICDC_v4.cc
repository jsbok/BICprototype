#include "koBICDetectorConstruction.hh"
#include "koBICCellParameterisation.hh"
#include "koBICFilterParameterisation.hh"
#include "koBICMirrorParameterisation.hh"
#include "koBICSiPMSD.hh"
#include "G4Trd.hh"

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

G4ThreadLocal koBICMagneticField *koBICDetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager *koBICDetectorConstruction::fFieldMgr = 0;

int koBICDetectorConstruction::fNofRow = 1;
int koBICDetectorConstruction::fNofCol = 1;
int koBICDetectorConstruction::fNofModules = fNofRow * fNofCol;

koBICDetectorConstruction::koBICDetectorConstruction()
    : G4VUserDetectorConstruction(), fMessenger(0), fMaterials(NULL) {
  DefineCommands();
  DefineMaterials();

  PMTT = 0.3 * mm;
  // filterT = 0.01*mm;

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

  checkOverlaps = false; // 겹침 확인을 위해 true로 켜두는 것을 권장합니다.

  G4VSolid *worldSolid = new G4Box("worldBox", 10. * m, 10. * m, 10. * m);
  worldLogical = new G4LogicalVolume(worldSolid, FindMaterial("G4_Galactic"), "worldLogical");
  G4VPhysicalVolume *worldPhysical =
      new G4PVPlacement(0, G4ThreeVector(), worldLogical, "worldPhysical", 0, false, 0, checkOverlaps);

  fFrontL = 1000.; 
  fTowerDepth = 700.; 
  
  // XML 구조에 따른 파라미터 적용
  fFiber_vert_dis = 1.22; // spacing_z
  fFiber_hori_dis = 1.35; // spacing_x
  
  // 모듈의 총 두께(X축 방향 적층) 계산: 
  // 5 * (17.0 + 21.730) + 1 * 17.0 + 6 * 20.889 = 335.984 mm
  fModuleH = 133.19 + 146.223; // = 279.413 mm
  fModuleW = 710.0; // 가로 폭 (필요시 수정)
  fFiberUnitH = 1.;

  G4double rmin = 904.485 * mm; 
G4double rmax = 1037.675 * mm + (20.889 * 7) * mm; 
G4double totalLength = rmax - rmin * mm; // 70cm + 30cm

// 48각 중 한 조각의 너비 계산
G4double dPhi = (360. / 48.) * deg;
G4double dx1 = rmin * std::tan(dPhi/2.); // 안쪽 짧은 변 절반
G4double dx2 = rmax * std::tan(dPhi/2.); // 바깥쪽 긴 변 절반
G4double dy  = fModuleW / 2.0;           // 모듈의 두께(높이) 절반

// G4Trd 생성 (이름, x_half_small, x_half_large, y_half_small, y_half_large, z_half)
// Geant4의 Trd는 z축을 따라 단면적이 변합니다.
G4VSolid* moduleSolid = new G4Trd("moduleSolid", dx1, dx2, dy, dy, totalLength/2.);

// 재질은 내부를 채울 것이므로 AIR 또는 Galactic(진공)으로 설정
ModuleLogical[0] = new G4LogicalVolume(moduleSolid, FindMaterial("G4_AIR"), "ModuleLogical");

// World에 배치
new G4PVPlacement(0, G4ThreeVector(0,0,0), ModuleLogical[0], "ModulePhys", worldLogical, false, 0, checkOverlaps);
  G4cout << "Step 0" << G4endl;

  G4double fGlue_thickness = 0.233;

  doFiber = true;
  doPMT = true;
  doGlue = false;

  // 광섬유 반경 설정 (XML 반영)
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
  
  gluebox = new G4Box("gluebox", (fGlue_thickness / 2) * mm, (fFiber_hori_dis / 2) * mm, fTowerDepth / 2.);

  G4Tubs* fiberClad = new G4Tubs("fiberClad", 0, clad_S_rMax, 1000.*mm, 0, 360.*deg);

  tGlueSubtraction = new G4SubtractionSolid("glueSub", glueBoxSolid, fiberClad2, 0, G4ThreeVector(0., 0., 0.)
                );
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
  G4cout << "Step 0.5: Construction Done" << G4endl;
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

            // Cathode(실리콘) 찾기
            G4String lvNameL = "ModuleC_Cath" + std::to_string(2 * i);
            G4String lvNameR = "ModuleC_Cath" + std::to_string(2 * i + 1);
            G4LogicalVolume* targetLV_L = lvStore->GetVolume(lvNameL);
            G4LogicalVolume* targetLV_R = lvStore->GetVolume(lvNameR);

            // Cell(유리) 찾기
            G4String lvNameL_Cell = "ModuleC_Cell" + std::to_string(2 * i);
            G4String lvNameR_Cell = "ModuleC_Cell" + std::to_string(2 * i + 1);
            G4LogicalVolume* targetCell_L = lvStore->GetVolume(lvNameL_Cell);
            G4LogicalVolume* targetCell_R = lvStore->GetVolume(lvNameR_Cell);

            // 왼쪽 할당
            if (targetLV_L) {
                targetLV_L->SetSensitiveDetector(SiPMSD_L);
                if (targetCell_L) targetCell_L->SetSensitiveDetector(SiPMSD_L);
//                G4cout << "SD attached to Cath/Cell L: " << i << G4endl;
            }

            // 오른쪽 할당
            if (targetLV_R) {
                targetLV_R->SetSensitiveDetector(SiPMSD_R);
                if (targetCell_R) targetCell_R->SetSensitiveDetector(SiPMSD_R);
 //               G4cout << "SD attached to Cath/Cell R: " << i << G4endl;
            }
        }
    }
}




void koBICDetectorConstruction::ModuleBuild(
    G4LogicalVolume *ModuleLogical_[], 
    G4LogicalVolume *PMTGLogical_[],
    G4LogicalVolume *PMTfilterLogical_[], 
    G4LogicalVolume *PMTcellLogical_[],
    G4LogicalVolume *PMTcathLogical_[],
    std::vector<G4LogicalVolume *> fiberUnitIntersection_[],
    std::vector<G4LogicalVolume *> fiberCladIntersection_[],
    std::vector<G4LogicalVolume *> fiberCoreIntersection_[],
    std::vector<koBICInterface::koBICModuleProperty> &ModuleProp_) {
G4Material* pbMat = FindMaterial("Lead");
    G4Material* vacMat = FindMaterial("G4_Galactic");

    // 1. 엄마 사다리꼴의 치수를 그대로 가져옵니다.
    // dx1(앞면 너비), dx2(뒷면 너비), dy(높이), dz(길이 절반)
    G4Trd* parentTrd = (G4Trd*)ModuleLogical_[0]->GetSolid();
    G4double pDx1 = parentTrd->GetXHalfLength1();
    G4double pDx2 = parentTrd->GetXHalfLength2();
    G4double pDy  = parentTrd->GetYHalfLength1() - 5;
    G4double pDz  = parentTrd->GetZHalfLength(); // 보통 500mm

    G4double currentZ = -pDz; // -500mm부터 시작
    int copyNo = 0;

    // 2. 레이어 생성 헬퍼 (작은 사다리꼴 직접 생성)
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
    return logVol; // <--- 생성된 포인터를 반환!
};


std::vector<G4LogicalVolume*> logicPbSFIL;
std::vector<G4LogicalVolume*> logicPbBulk;

for(int r = 0; r < 3; r++) {
    fVisAttrBlue->SetVisibility(false);
    PlaceTrdLayer("Vac_SFIL", 17.0*mm, vacMat, fVisAttrBlue);
    
    fVisAttrGray->SetColour(G4Colour(0.5, 0.5, 0.5, 0.5)); 
    fVisAttrGray->SetForceSolid(true);
    
    // [수정] 생성된 포인터를 변수에 저장하고 벡터에 push_back 합니다.
    G4LogicalVolume* tmpPbSFIL = PlaceTrdLayer("Pb_SFIL", 21.73*mm, pbMat, fVisAttrGray);
    logicPbSFIL.push_back(tmpPbSFIL); 
}

fVisAttrBlue->SetVisibility(false);
PlaceTrdLayer("Vac_SFIL_End", 17.0*mm, vacMat, fVisAttrBlue);

G4double bulkHalfHeight = (300.0 * mm) / 2.0; 

// --- Bulk 레이어 적층 ---
for(int r = 0; r < 7; r++) {
    fVisAttrGray->SetVisibility(true);
    fVisAttrGray->SetForceSolid(true);
    
    // [수정] 마찬가지로 포인터를 저장하여 벡터에 담습니다.
    G4LogicalVolume* tmpPbBulk = PlaceTrdLayer("Pb_Bulk", 20.889*mm, pbMat, fVisAttrGray, bulkHalfHeight);
    logicPbBulk.push_back(tmpPbBulk);
}

    // 5. Fiber 배치 (Fiber도 Z축 방향으로 길게 뻗어야 함)
    FiberImplement(0,  logicPbSFIL, logicPbBulk, 
               fiberUnitIntersection_, fiberCladIntersection_, fiberCoreIntersection_);

    // =========================================================
    // 3. PMT 및 속성 설정
    // =========================================================
if (doPMT) {
    G4double sipmSize = 16.0 * mm;
    G4double sipmThick = PMTT;      // 예: 0.3mm (전체 두께)
    
    // [수정 1] 실리콘(Cathode) 두께를 Cell(전체 두께)과 완벽히 동일하게 맞춤
    // 이렇게 하면 물리적인 '유리창' 층이 사라지고, 실리콘이 100% 공간을 차지합니다.
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

            // 이름 정의 (Layer 관련 이름 제거)
            G4String nameCellL = "ModuleC_Cell" + std::to_string(2 * s);
            G4String nameCellR = "ModuleC_Cell" + std::to_string(2 * s + 1);
            G4String nameCathL = "ModuleC_Cath" + std::to_string(2 * s);
            G4String nameCathR = "ModuleC_Cath" + std::to_string(2 * s + 1);

            // 1. Glass Cell 생성
            PMTcellLogical_[2 * s] = new G4LogicalVolume(sipmCellSolid, FindMaterial("Glass"), nameCellL);
            PMTcellLogical_[2 * s + 1] = new G4LogicalVolume(sipmCellSolid, FindMaterial("Glass"), nameCellR);
            
            // 2. Silicon Cathode 생성 및 Cell 내부에 배치
            PMTcathLogical_[2 * s] = new G4LogicalVolume(sipmCathSolid, FindMaterial("Silicon"), nameCathL);
            PMTcathLogical_[2 * s + 1] = new G4LogicalVolume(sipmCathSolid, FindMaterial("Silicon"), nameCathR);

            new G4PVPlacement(0, G4ThreeVector(), PMTcathLogical_[2 * s], nameCathL + "_phys", PMTcellLogical_[2 * s], false, 2 * s, checkOverlaps);
            new G4PVPlacement(0, G4ThreeVector(), PMTcathLogical_[2 * s + 1], nameCathR + "_phys", PMTcellLogical_[2 * s + 1], false, 2 * s + 1, checkOverlaps);

            // 3. Glass Cell을 ModuleLogical_[0]에 직접 배치 (Air Layer 생략)
            G4double xPos = -currentDx + (xIdx + 0.5) * stepX;
            G4double yPos_L = -(currentDy + sipmThick / 2.);
            G4double yPos_R =  (currentDy + sipmThick / 2.);

            new G4PVPlacement(0, G4ThreeVector(xPos, yPos_L, zCenter), PMTcellLogical_[2 * s], nameCellL + "_phys", ModuleLogical_[0], false, 2 * s, checkOverlaps);
            new G4PVPlacement(0, G4ThreeVector(xPos, yPos_R, zCenter), PMTcellLogical_[2 * s + 1], nameCellR + "_phys", ModuleLogical_[0], false, 2 * s + 1, checkOverlaps);

            // 4. Optical Surface 설정
            new G4LogicalSkinSurface(nameCathL + "_surf", PMTcathLogical_[2 * s], FindSurface("SiPMSurf"));
            new G4LogicalSkinSurface(nameCathR + "_surf", PMTcathLogical_[2 * s + 1], FindSurface("SiPMSurf"));

auto surfL = G4LogicalSkinSurface::GetSurface(PMTcathLogical[2*s]);
G4cout << "L skin = " << surfL
       << " on " << PMTcathLogical[2*s]->GetName()
       << G4endl;

auto surfR = G4LogicalSkinSurface::GetSurface(PMTcathLogical[2*s+1]);
G4cout << "R skin = " << surfR
       << " on " << PMTcathLogical[2*s+1]->GetName()
       << G4endl;

            PMTcathLogical_[2 * s]->SetVisAttributes(fVisAttrGreen);
            PMTcathLogical_[2 * s + 1]->SetVisAttributes(fVisAttrGreen);
            
if(PMTcathLogical[2*s])
{
    auto skinL =
      G4LogicalSkinSurface::GetSurface(PMTcathLogical[2*s]);

    G4cout << "L skin = "
           << skinL
           << " on "
           << PMTcathLogical[2*s]->GetName()
           << G4endl;
}

if(PMTcathLogical[2*s+1])
{
    auto skinR =
      G4LogicalSkinSurface::GetSurface(PMTcathLogical[2*s+1]);

    G4cout << "R skin = "
           << skinR
           << " on "
           << PMTcathLogical[2*s+1]->GetName()
           << G4endl;
}
        }
        pbLayerIdx++;
    }
        runningZ += thick;
    }
}

for (int s = 0; s < 50; s++) {
    koBICInterface::koBICModuleProperty ModulePropSingle;
    ModulePropSingle.towerXY = fTowerXY;
    ModulePropSingle.ModuleNum = s; // 0~49번까지 각각 부여
    ModuleProp_.push_back(ModulePropSingle);
    }

G4cout << "Step 4" << G4endl;
}

void koBICDetectorConstruction::DefineCommands() {}
void koBICDetectorConstruction::FiberImplement(
    G4int i, 
    std::vector<G4LogicalVolume*> logicPbSFIL, 
    std::vector<G4LogicalVolume*> logicPbBulk,
    std::vector<G4LogicalVolume *> fiberUnitIntersection__[],
    std::vector<G4LogicalVolume *> fiberCladIntersection__[],
    std::vector<G4LogicalVolume *> fiberCoreIntersection__[]) {

    if (!doFiber) return;

  fVisAttrSkyBlue = new G4VisAttributes(G4Colour(0.5, 0.8, 0.9, 0.5));
  fVisAttrSkyBlue->SetVisibility(true);
  fVisAttrSkyBlue->SetForceSolid(true);

    G4RotationMatrix* fiberRot = new G4RotationMatrix();
    fiberRot->rotateX(90*deg); 

    G4double pitch_horizontal = 1.35 * mm; 
    G4double pitch_vertical = 1.22 * mm;   
    int fiberId = 0;

auto FillFibersInMother = [&](G4LogicalVolume* motherLog, G4double thickness, G4LogicalVolume* targetGlueLog) {

    G4Trd* motherSolid = (G4Trd*)motherLog->GetSolid();
    G4double h1 = motherSolid->GetXHalfLength1();
    G4double h2 = motherSolid->GetXHalfLength2();
    G4double zHalf = motherSolid->GetZHalfLength();

    int numRows = std::floor(thickness / pitch_vertical);
    G4double startV = -((numRows - 1) * pitch_vertical) / 2.0;

    // row마다 numCols 바뀌면 위상 깨지므로 최대폭 기준 고정
    G4double maxHalfWidth = std::max(h1, h2);
    int fixedNumCols = std::floor((2.0 * maxHalfWidth - 1.0 * mm) / pitch_horizontal);
    G4double startH = -((fixedNumCols - 1) * pitch_horizontal) / 2.0;

    for (int k = 0; k < numRows; k++) {

        G4double localV = startV + k * pitch_vertical;
        G4double currentHalfWidth = h1 + (h2 - h1) * (localV + zHalf) / (2.0 * zHalf);
        bool isShifted = (k % 2 != 0);

        for (int j = 0; j < fixedNumCols; j++) {

            G4double localH = startH + j * pitch_horizontal;
            if (isShifted) localH += pitch_horizontal / 2.0;

            // 경계 넘는 fiber 제거
            if (std::abs(localH) + 0.553 * mm > currentHalfWidth - 0.1 * mm) continue;

            new G4PVPlacement(fiberRot, G4ThreeVector(localH, 0, localV), targetGlueLog, "fiberGlue_Phys", motherLog, false, fiberId, false);

            if (doGlue) {

                tGlueIntersection = new G4IntersectionSolid("glue", tGlueSubtraction, targetGlueLog->GetSolid(), 0, G4ThreeVector(0,0,0));

		G4LogicalVolume* glueLog = new G4LogicalVolume(tGlueSubtraction, FindMaterial("Glue"), "glue_Log");
                glueIntersection__[i].push_back(glueLog);

		new G4PVPlacement(fiberRot, G4ThreeVector(localH, 0., localV), glueLog, "glue_Phys", motherLog, false, fiberId, false);

                glueLog->SetVisAttributes(fVisAttrSkyBlue);
            }

            fiberId++;
        }
    }
};
    // 레이어 순회
for (auto pbLog : logicPbSFIL) {
    // 1. 설계도 한 세트 만들기 (레이어당 1회)
    G4LogicalVolume* glueLog_SFIL = new G4LogicalVolume(fiberGlue_SFIL, FindMaterial("Glue"), "fGlue_SFIL_L");
    G4LogicalVolume* cladLog_SFIL = new G4LogicalVolume(fiberClad_SFIL, FindMaterial("PMMA"), "fClad_SFIL_L");
    G4LogicalVolume* coreLog_SFIL = new G4LogicalVolume(fiberCore_SFIL, FindMaterial("Polystyrene"), "fCore_SFIL_L");

    // 2. 계층 조립 (여기서 미리 다 박아버립니다)
    // Core -> Clad
    new G4PVPlacement(0, G4ThreeVector(0,0,0), coreLog_SFIL, "fCore_P", cladLog_SFIL, false, 0, false);
    // Clad -> Glue
    new G4PVPlacement(0, G4ThreeVector(0,0,0), cladLog_SFIL, "fClad_P", glueLog_SFIL, false, 0, false);

    // 시각화
    glueLog_SFIL->SetVisAttributes(fVisAttrSkyBlue);
    cladLog_SFIL->SetVisAttributes(fVisAttrGray);
    coreLog_SFIL->SetVisAttributes(fVisAttrOrange);

    // 3. 람다 호출 (이제 '완성된' Glue 설계도를 넘깁니다)
    FillFibersInMother(pbLog, 21.73*mm, glueLog_SFIL);
}

for (auto pbLog : logicPbBulk) {
    // 1. 벌크용 설계도 세트 만들기 (레이어당 1회 생성)
    // fiberGlue_Solid는 위에서 SFIL과 공유하거나 새로 정의한 것을 사용
    G4LogicalVolume* glueLog_Bulk = new G4LogicalVolume(fiberGlue_Bulk, FindMaterial("Glue"), "fGlue_Bulk_L");
    G4LogicalVolume* cladLog_Bulk = new G4LogicalVolume(fiberClad_Bulk, FindMaterial("PMMA"), "fClad_Bulk_L");
    G4LogicalVolume* coreLog_Bulk = new G4LogicalVolume(fiberCore_Bulk, FindMaterial("Polystyrene"), "fCore_Bulk_L");

    // 2. 벌크 계층 조립 (Core -> Clad -> Glue)
    new G4PVPlacement(0, G4ThreeVector(0,0,0), coreLog_Bulk, "fCore_P", cladLog_Bulk, false, 0, false);
    new G4PVPlacement(0, G4ThreeVector(0,0,0), cladLog_Bulk, "fClad_P", glueLog_Bulk, false, 0, false);

    // 시각화 (SFIL과 동일하게 설정하거나 구분 가능하게 변경)
    glueLog_Bulk->SetVisAttributes(fVisAttrSkyBlue);
    cladLog_Bulk->SetVisAttributes(fVisAttrGray);
    coreLog_Bulk->SetVisAttributes(fVisAttrOrange);

    // 3. 람다 호출 (벌크 두께인 20.889*mm와 조립된 glueLog_Bulk 전달)
    FillFibersInMother(pbLog, 20.889*mm, glueLog_Bulk); 
}
}
