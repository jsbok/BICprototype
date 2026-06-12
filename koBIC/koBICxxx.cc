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

koBICDetectorConstruction::~koBICDetectorConstruction() {
  delete fMessenger;
  delete fMaterials;

  delete fVisAttrOrange;
  delete fVisAttrBlue;
  delete fVisAttrGray;
  delete fVisAttrGreen;
  // delete fVisAttrSkyBlue;
}

void koBICDetectorConstruction::DefineMaterials() {
  fMaterials = koBICMaterials::GetInstance();
}

G4VPhysicalVolume *koBICDetectorConstruction::Construct() {
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  checkOverlaps = true; // 겹침 확인을 위해 true로 켜두는 것을 권장합니다.

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

  fiberUnit = new G4Box("fiber_SQ", (fFiberUnitH / 2) * mm, (1. / 2) * mm, (fTowerDepth / 2) * mm);
fiberClad_SFIL = new G4Tubs("fiberClad_SFIL", 0, clad_S_rMax, 700./ 2., 0 * deg, 360. * deg);
fiberCore_SFIL = new G4Tubs("fiberCore_SFIL", 0, core_S_rMax, 700. / 2., 0 * deg, 360. * deg);

fiberClad_Bulk = new G4Tubs("fiberClad_Bulk", 0, clad_S_rMax, 300. / 2., 0 * deg, 360. * deg);
fiberCore_Bulk = new G4Tubs("fiberCore_Bulk", 0, core_S_rMax, 300. / 2., 0 * deg, 360. * deg);
  
  gluebox = new G4Box("gluebox", (fGlue_thickness / 2) * mm, (fFiber_hori_dis / 2) * mm, fTowerDepth / 2.);

G4Tubs* fiberClad = new G4Tubs("fiberClad", 0, clad_S_rMax, 1000.*mm, 0, 360.*deg);
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
  G4cout << "Step 0.5: Construction Done" << G4endl;
  return worldPhysical;

}

void koBICDetectorConstruction::ConstructSDandField() {
    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();

    if (doPMT) {
        for (int i = 0; i < 50; i++) {
            // SD 이름과 컬렉션 이름을 고유하게 만듭니다 (Module0, Module1... / ModuleC0, ModuleC1...)
            G4String nameL = "ModuleSD_L_" + std::to_string(i);
            G4String nameR = "ModuleSD_R_" + std::to_string(i);
            G4String collL = "ModuleC" + std::to_string(2 * i);
            G4String collR = "ModuleC" + std::to_string(2 * i + 1);


	    int idL = 2 * i;      // 0, 2, 4...
            int idR = 2 * i + 1;  // 1, 3, 5...
        
            // SD 생성
            koBICSiPMSD* SiPMSD_L = new koBICSiPMSD(nameL, collL, 0, fModuleProp.at(i), idL);
            koBICSiPMSD* SiPMSD_R = new koBICSiPMSD(nameR, collR, 1, fModuleProp.at(i), idR);
            SDman->AddNewDetector(SiPMSD_L);
            SDman->AddNewDetector(SiPMSD_R);

            // 볼륨 이름은 로그에서 확인된 "ModuleC_Cath" + 숫자
            G4LogicalVolume* targetLV_L = lvStore->GetVolume("ModuleC_Cath" + std::to_string(2 * i));
            G4LogicalVolume* targetLV_R = lvStore->GetVolume("ModuleC_Cath" + std::to_string(2 * i + 1));

            if (targetLV_L) targetLV_L->SetSensitiveDetector(SiPMSD_L);
            if (targetLV_R) targetLV_R->SetSensitiveDetector(SiPMSD_R);
        }
        G4cout << ">>> SD assigned with Collection names ModuleC0~99" << G4endl;
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
    
    fVisAttrGray->SetColour(G4Colour(0.5, 0.5, 0.5, 1.0)); 
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
    G4double cathThick = PMTT-0.001*mm;      

    std::vector<G4double> layerThicks = {
        17.0*mm, 21.73*mm, 17.0*mm, 21.73*mm, 17.0*mm, 21.73*mm, 17.0*mm,
        20.889*mm, 20.889*mm, 20.889*mm, 20.889*mm, 20.889*mm, 20.889*mm, 20.889*mm
    };

    G4double runningZ = -pDz;
    int pbLayerIdx = 0;

    // 1. 솔리드 정의 (두께가 같아짐)
    G4VSolid* sipmCellSolid = new G4Box("SiPMCellSolid", sipmSize/2., sipmThick/2., sipmSize/2.);
    G4VSolid* sipmCathSolid = new G4Box("SiPMCathSolid", sipmSize/2., cathThick/2., sipmSize/2.);

    for (int i = 0; i < 14; i++) {
        G4double thick = layerThicks[i];
        G4double zCenter = runningZ + thick / 2.0;
        bool isPbLayer = (i == 1 || i == 3 || i == 5 || i >= 7);

        if (isPbLayer) {
            G4double currentDy = (i >= 7) ? (150.0 * mm) : pDy;
            G4double currentDx = pDx1 + (pDx2 - pDx1) * ((zCenter + pDz) / (2.0 * pDz));
            G4double stepX = (currentDx * 2.0) / 5.0;

            for (int xIdx = 0; xIdx < 5; xIdx++) {
                int s = pbLayerIdx * 5 + xIdx;

                G4String nameCellL = "ModuleC_Cell" + std::to_string(2 * s);
                G4String nameCellR = "ModuleC_Cell" + std::to_string(2 * s + 1);
                G4String nameCathL = "ModuleC_Cath" + std::to_string(2 * s);
                G4String nameCathR = "ModuleC_Cath" + std::to_string(2 * s + 1);

                // 2. Cell 생성 (원본의 껍데기 역할)
                G4LogicalVolume* PMTcellLogicalL = new G4LogicalVolume(sipmCellSolid, FindMaterial("Glass"), nameCellL);
                G4LogicalVolume* PMTcellLogicalR = new G4LogicalVolume(sipmCellSolid, FindMaterial("Glass"), nameCellR);

                // 3. Cathode 생성 (Silicon 덩어리)
                PMTcathLogical[2 * s]     = new G4LogicalVolume(sipmCathSolid, FindMaterial("Silicon"), nameCathL);
                PMTcathLogical[2 * s + 1] = new G4LogicalVolume(sipmCathSolid, FindMaterial("Silicon"), nameCathR);

                // 4. Cathode 배치 
                // [수정 2] yOffset을 완전히 제거하고 (0, 0, 0) 정중앙에 배치합니다.
                // 이로써 Cathode가 Cell을 빈틈없이 100% 덮어씌웁니다.
                new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), PMTcathLogical[2 * s], 
                                  nameCathL + "_phys", PMTcellLogicalL, false, 2*s, checkOverlaps);
                
                new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), PMTcathLogical[2 * s + 1], 
                                  nameCathR + "_phys", PMTcellLogicalR, false, 2*s+1, checkOverlaps);

                // 5. 좌표 계산 및 Cell을 Module에 배치 (Fiber 끝단과 딱 맞닿는 위치)
                G4double xPos = -currentDx + (xIdx + 0.5) * stepX;
                G4double yPos_L = -(currentDy + sipmThick/2.);
                G4double yPos_R =  (currentDy + sipmThick/2.);

                new G4PVPlacement(0, G4ThreeVector(xPos, yPos_L, zCenter), PMTcellLogicalL, 
                                  nameCellL + "_phys", ModuleLogical_[0], false, 2*s, checkOverlaps);
                new G4PVPlacement(0, G4ThreeVector(xPos, yPos_R, zCenter), PMTcellLogicalR, 
                                  nameCellR + "_phys", ModuleLogical_[0], false, 2*s+1, checkOverlaps);

                // 6. 광학 표면 및 시각화 설정 (Cathode에 적용)
                // [수정 3] 주석 해제. 이 코드가 켜져 있어야 metal 표면에서 광자가 멈추지 않고 즉시 소멸(Kill)하며 검출됩니다.
                new G4LogicalSkinSurface(nameCathL + "_surf", PMTcathLogical[2 * s], FindSurface("SiPMSurf"));
                new G4LogicalSkinSurface(nameCathR + "_surf", PMTcathLogical[2 * s + 1], FindSurface("SiPMSurf"));

                PMTcathLogical[2 * s]->SetVisAttributes(fVisAttrGreen);
                PMTcathLogical[2 * s + 1]->SetVisAttributes(fVisAttrGreen);
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
    G4cout << "Step 3.6" << G4endl;
/*
    if (doPMT) {
        G4VSolid *SiPMlayerSolid = new G4Box("SiPMlayerSolid", (fModuleH / 2.) * mm, (fModuleW / 2.) * mm, (PMTT / 2.) * mm);
        G4LogicalVolume *SiPMlayerLogical_left = new G4LogicalVolume(SiPMlayerSolid, FindMaterial("G4_AIR"), "SiPMlayerLogical");
        G4LogicalVolume *SiPMlayerLogical_right = new G4LogicalVolume(SiPMlayerSolid, FindMaterial("G4_AIR"), "SiPMlayerLogical");
        
        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), SiPMlayerLogical_left, "SiPMlayerPhysical", PMTGLogical_[0], false, 0, checkOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), SiPMlayerLogical_right, "SiPMlayerPhysical", PMTGLogical_[1], false, 0, checkOverlaps);

        G4VSolid *PMTcellSolid = new G4Box("PMTcellSolid", 1.1 / 2. * mm, 1.1 / 2. * mm, PMTT / 2. * mm);
        PMTcellLogical_[0] = new G4LogicalVolume(PMTcellSolid, FindMaterial("Glass"), "PMTcellLogical_");
        PMTcellLogical_[1] = new G4LogicalVolume(PMTcellSolid, FindMaterial("Glass"), "PMTcellLogical_");

        koBICCellParameterisation *PMTcellParam = new koBICCellParameterisation(fTowerXY.first, fTowerXY.second, fModuleH, fModuleW);
        new G4PVParameterised("PMTcellPhysical_left", PMTcellLogical_[0], SiPMlayerLogical_left, kXAxis, fTowerXY.first * fTowerXY.second - fTowerXY.first / 2, PMTcellParam);
        new G4PVParameterised("PMTcellPhysical_right", PMTcellLogical_[1], SiPMlayerLogical_right, kXAxis, fTowerXY.first * fTowerXY.second - fTowerXY.first / 2, PMTcellParam);

        G4VSolid *PMTcathSolid = new G4Box("PMTcathSolid", 1.1 / 2. * mm, 1.1 / 2. * mm, PMTT / 2. * mm);
        PMTcathLogical_[0] = new G4LogicalVolume(PMTcathSolid, FindMaterial("Silicon"), "PMTcathLogical_");
        PMTcathLogical_[1] = new G4LogicalVolume(PMTcathSolid, FindMaterial("Silicon"), "PMTcathLogical_");
        
        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), PMTcathLogical_[0], "PMTcathPhysical", PMTcellLogical_[0], false, 0, checkOverlaps);
        new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), PMTcathLogical_[1], "PMTcathPhysical", PMTcellLogical_[1], false, 0, checkOverlaps);
        
        new G4LogicalSkinSurface("Photocath_surf_left", PMTcathLogical_[0], FindSurface("SiPMSurf"));
        new G4LogicalSkinSurface("Photocath_surf_right", PMTcathLogical_[1], FindSurface("SiPMSurf"));

        PMTcathLogical_[0]->SetVisAttributes(fVisAttrGreen);
        PMTcathLogical_[1]->SetVisAttributes(fVisAttrGreen);
    }*/
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

    G4RotationMatrix* fiberRot = new G4RotationMatrix();
    fiberRot->rotateX(90*deg); 

    G4double pitch_horizontal = 1.35 * mm; 
    G4double pitch_vertical = 1.22 * mm;   
    int fiberId = 0;

// 1. SFIL용 Logical Volume 미리 만들기 (딱 한 번만)
/*
G4LogicalVolume* cladLog_SFIL = new G4LogicalVolume(fiberClad_SFIL, FindMaterial("PMMA"), "fiberClad_SFIL_Log");
G4LogicalVolume* coreLog_SFIL = new G4LogicalVolume(fiberCore_SFIL, FindMaterial("Polystyrene"), "fiberCore_SFIL_Log");
new G4PVPlacement(0, G4ThreeVector(0,0,0), coreLog_SFIL, "fiberCore_Phys", cladLog_SFIL, false, 0, false);
*/
// 2. Bulk용 Logical Volume 미리 만들기 (딱 한 번만)
/*
G4LogicalVolume* cladLog_Bulk = new G4LogicalVolume(fiberClad_Bulk, FindMaterial("PMMA"), "fiberClad_Bulk_Log");
G4LogicalVolume* coreLog_Bulk = new G4LogicalVolume(fiberCore_Bulk, FindMaterial("Polystyrene"), "fiberCore_Bulk_Log");
new G4PVPlacement(0, G4ThreeVector(0,0,0), coreLog_Bulk, "fiberCore_Phys", cladLog_Bulk, false, 0, false);
*/
         //   fiberCladIntersection__[i].push_back(cladLog);
         //   fiberCoreIntersection__[i].push_back(coreLog);
    /*            cladLog_SFIL->SetVisAttributes(fVisAttrGray);
                coreLog_SFIL->SetVisAttributes(fVisAttrOrange);
                cladLog_Bulk->SetVisAttributes(fVisAttrGray);
		coreLog_Bulk->SetVisAttributes(fVisAttrOrange);
	*/	
auto FillFibersInMother = [&](G4LogicalVolume* motherLog, G4double thickness, G4LogicalVolume* targetCladLog) {
    
    G4Trd* motherSolid = (G4Trd*)motherLog->GetSolid();
    
    // 사다리꼴의 치수 정보를 가져옵니다. (Geant4 Trd는 보통 X1, X2가 Half Length임)
    G4double h1 = motherSolid->GetXHalfLength1(); // 아랫면 절반 폭
    G4double h2 = motherSolid->GetXHalfLength2(); // 윗면 절반 폭
    G4double zHalf = motherSolid->GetZHalfLength(); // 두께의 절반

    // 1. 세로(두께, Z축) 방향 층수 계산
    int numRows = std::floor(thickness / pitch_vertical);
    G4double startV = -((numRows - 1) * pitch_vertical) / 2.0;

    for (int k = 0; k < numRows; k++) {
        G4double localV = startV + k * pitch_vertical; 
        
        // [중요] 해당 높이(localV)에서의 실제 사다리꼴 절반 폭 계산 (선형 보간)
        // localV가 -zHalf일 때 h1, +zHalf일 때 h2가 되는 공식입니다.
        G4double currentHalfWidth = h1 + (h2 - h1) * (localV + zHalf) / (2.0 * zHalf);
        G4double currentFullWidth = currentHalfWidth * 2.0;

        bool isShifted = (k % 2 != 0);

        // 2. 해당 층의 폭에 맞춘 가로 개수 계산
        int numCols = std::floor((currentFullWidth - 1.0*mm) / pitch_horizontal); 
        G4double startH = -((numCols - 1) * pitch_horizontal) / 2.0;

        for (int j = 0; j < numCols; j++) {
            G4double localH = startH + j * pitch_horizontal;
            if (isShifted) localH += pitch_horizontal / 2.0;

            // [안전 검사] 현재 층의 실제 폭을 넘지 않는지 확인
            if (std::abs(localH) + clad_S_rMax > currentHalfWidth - 0.5*mm) continue;
/*
            // 3. 물리 배치 (Clad 생성 후 그 안에 Core 배치)
            G4LogicalVolume* cladLog = new G4LogicalVolume(cladSolid, FindMaterial("PMMA"), "fiberClad_Log");
            G4LogicalVolume* coreLog = new G4LogicalVolume(coreSolid, FindMaterial("Polystyrene"), "fiberCore_Log");

            fiberCladIntersection__[i].push_back(cladLog);
            fiberCoreIntersection__[i].push_back(coreLog);
*/
            // X=localH(가로), Y=0(파이버 길이방향), Z=localV(두께방향)
		new G4PVPlacement(fiberRot, G4ThreeVector(localH, 0, localV), targetCladLog,
                   		   "fiberClad_Phys", motherLog, false, fiberId, false);
/*
            new G4PVPlacement(0, G4ThreeVector(0, 0, 0), coreLog,
                              "fiberCore_Phys", cladLog, false, fiberId, false);
                              
                cladLog->SetVisAttributes(fVisAttrGray);
                coreLog->SetVisAttributes(fVisAttrOrange);
*/

                // 3. Glue 로직 (기존 벡터 참조 방식 대신 현재 local 좌표 사용)
                if (doGlue) {
                    // 주의: Glue의 Mother 볼륨과 좌표계가 무엇인지에 따라 MotherLog 혹은 ModuleLogical을 선택해야 합니다.
                    // 여기서는 파이버와 동일하게 납 레이어(motherLog)에 넣는 것으로 예시를 듭니다.
                    tGlueIntersection = new G4IntersectionSolid("glue", tGlueSubtraction, targetCladLog->GetSolid(), 0, G4ThreeVector(0,0,0));
                    
                    G4LogicalVolume* glueLog = new G4LogicalVolume(tGlueIntersection, FindMaterial("G4_Galactic"), "glue_Log");
                    glueIntersection__[i].push_back(glueLog);
                    
                    new G4PVPlacement(fiberRot, G4ThreeVector(localV, 0, localH), glueLog,
                                      "glue_Phys", motherLog, false, fiberId, false);

                    glueLog->SetVisAttributes(fVisAttrBlue);
                }

                fiberId++; // 한 개의 파이버 세트 배치가 끝날 때마다 증가
            } // j loop 끝
        } // k loop 끝
    };

    // 레이어 순회
    for (auto pbLog : logicPbSFIL) {
    G4LogicalVolume* cladLog_SFIL = new G4LogicalVolume(fiberClad_SFIL, FindMaterial("PMMA"), "fiberClad_SFIL_Log");
G4LogicalVolume* coreLog_SFIL = new G4LogicalVolume(fiberCore_SFIL, FindMaterial("Polystyrene"), "fiberCore_SFIL_Log");
new G4PVPlacement(0, G4ThreeVector(0,0,0), coreLog_SFIL, "fiberCore_Phys", cladLog_SFIL, false, 0, false);
                cladLog_SFIL->SetVisAttributes(fVisAttrGray);
                coreLog_SFIL->SetVisAttributes(fVisAttrOrange);

	FillFibersInMother(pbLog, 21.73*mm, cladLog_SFIL);
    }
    for (auto pbLog : logicPbBulk) {
    G4LogicalVolume* cladLog_Bulk = new G4LogicalVolume(fiberClad_Bulk, FindMaterial("PMMA"), "fiberClad_Bulk_Log");
G4LogicalVolume* coreLog_Bulk = new G4LogicalVolume(fiberCore_Bulk, FindMaterial("Polystyrene"), "fiberCore_Bulk_Log");
new G4PVPlacement(0, G4ThreeVector(0,0,0), coreLog_Bulk, "fiberCore_Phys", cladLog_Bulk, false, 0, false);
                cladLog_Bulk->SetVisAttributes(fVisAttrGray);
		coreLog_Bulk->SetVisAttributes(fVisAttrOrange);

	FillFibersInMother(pbLog, 20.889*mm, cladLog_Bulk); 
    }
}
