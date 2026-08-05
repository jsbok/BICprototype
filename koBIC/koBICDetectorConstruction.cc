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
  G4VisAttributes* glueVis = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3)); 
  glueVis->SetForceSolid(true);
  glueVis->SetVisibility(true);
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
ModuleLogical[0] = new G4LogicalVolume(moduleSolid, FindMaterial("G4_Galactic"), "ModuleLogical");

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
fiberClad_SFIL = new G4Tubs("fiberClad_SFIL", 0, clad_S_rMax, 700. / 2., 0 * deg, 360. * deg);
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
    G4SDManager *SDman = G4SDManager::GetSDMpointer();
  
if (doPMT) {
        if (fModuleProp.size() < 50) {
            G4cout << "Error: fModuleProp size is " << fModuleProp.size() << ", expected 50." << G4endl;
            return;
        }

        for (int i = 0; i < 50; i++) {
            // [중요] 이름을 "Module" (대문자 M)로 통일합니다. 
            // 시스템이 "ModuleC0"을 찾고 있다면 정확히 대소문자가 일치해야 합니다.
            
            G4String nameLeft  = "Module" + std::to_string(2 * i);
            G4String collLeft  = "ModuleC" + std::to_string(2 * i);
            G4String nameRight = "Module" + std::to_string(2 * i + 1);
            G4String collRight = "ModuleC" + std::to_string(2 * i + 1);

            // Left SD
            koBICSiPMSD *SiPMSD_left = new koBICSiPMSD(nameLeft, collLeft, 0, fModuleProp.at(i));
            SDman->AddNewDetector(SiPMSD_left);

            // Right SD
            koBICSiPMSD *SiPMSD_right = new koBICSiPMSD(nameRight, collRight, 1, fModuleProp.at(i));
            SDman->AddNewDetector(SiPMSD_right);

            // PMTcathLogical에 SD 할당
            if (PMTcathLogical[2 * i]) {
                PMTcathLogical[2 * i]->SetSensitiveDetector(SiPMSD_left);
            }
            if (PMTcathLogical[2 * i + 1]) {
                PMTcathLogical[2 * i + 1]->SetSensitiveDetector(SiPMSD_right);
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
    G4double sipmThick = PMTT; // 0.3mm
    
    // 전체 14개 레이어 순서대로 두께 배열 (Vac 포함)
    std::vector<G4double> layerThicks = {
        17.0*mm, 21.73*mm, 17.0*mm, 21.73*mm, 17.0*mm, 21.73*mm, 17.0*mm, // 0~6: SFIL 구간 (Vac-Pb 교차 + Vac_End)
        20.889*mm, 20.889*mm, 20.889*mm, 20.889*mm, 20.889*mm, 20.889*mm, 20.889*mm // 7~13: Bulk 구간
    };

    G4double runningZ = -pDz; // -500mm 부터 시작
    int pbLayerIdx = 0;       // 납 레이어 번호 (0 ~ 9)

    G4VSolid* sipmSolid = new G4Box("SiPMSolid", sipmSize/2., sipmThick/2., sipmSize/2.);

    for (int i = 0; i < 14; i++) {
        G4double thick = layerThicks[i];
        G4double zCenter = runningZ + thick / 2.0;

        // 납(Pb) 레이어인지 확인: 인덱스 1, 3, 5 (Pb_SFIL) 및 7 이상 (Pb_Bulk)
        bool isPbLayer = (i == 1 || i == 3 || i == 5 || i >= 7);

        if (isPbLayer) {
            // 1. Y축 기준: Bulk 구간은 150mm, SFIL 구간은 pDy 사용
            G4double currentDy = (i >= 7) ? (300.0 * mm / 2.0) : pDy;

            // 2. X축 기준: 사다리꼴(Trd)의 기울기를 반영하여 현재 Z 위치에서의 절반 너비(Dx) 계산
            G4double currentDx = pDx1 + (pDx2 - pDx1) * ((zCenter + pDz) / (2.0 * pDz));
            G4double stepX = (currentDx * 2.0) / 5.0; // 5등분 간격

            for (int xIdx = 0; xIdx < 5; xIdx++) {
                int s = pbLayerIdx * 5 + xIdx; // 0 ~ 49

                G4String nameL = "ModuleC" + std::to_string(2 * s);
                G4String nameR = "ModuleC" + std::to_string(2 * s + 1);

                PMTcathLogical[2 * s]   = new G4LogicalVolume(sipmSolid, FindMaterial("Silicon"), nameL);
                PMTcathLogical[2 * s + 1] = new G4LogicalVolume(sipmSolid, FindMaterial("Silicon"), nameR);

                // X 좌표: 사다리꼴 빗면에 맞춰 현재 currentDx 범위 내에서 중앙 정렬
                G4double xPos = -currentDx + (xIdx + 0.5) * stepX;
                
                // Y 좌표: 현재 레이어의 Y 높이에 맞춰 모듈 외벽에 부착
                G4double yPos_L = -(currentDy + sipmThick/2. + 0.01*mm);
                G4double yPos_R =  (currentDy + sipmThick/2. + 0.01*mm);

               
                new G4PVPlacement(0, G4ThreeVector(xPos, yPos_L, zCenter), PMTcathLogical[2 * s], 
                                  nameL, ModuleLogical_[0], false, 2*s, checkOverlaps);
                new G4PVPlacement(0, G4ThreeVector(xPos, yPos_R, zCenter), PMTcathLogical[2 * s + 1], 
                                  nameR, ModuleLogical_[0], false, 2*s+1, checkOverlaps);

                new G4LogicalSkinSurface(nameL + "_surf", PMTcathLogical[2 * s], FindSurface("SiPMSurf"));
                new G4LogicalSkinSurface(nameR + "_surf", PMTcathLogical[2 * s + 1], FindSurface("SiPMSurf"));
                PMTcathLogical[2 * s]->SetVisAttributes(fVisAttrGreen);
                PMTcathLogical[2 * s + 1]->SetVisAttributes(fVisAttrGreen);
            }
            pbLayerIdx++; // SiPM이 부착된 납 레이어 카운트 증가
        }
        runningZ += thick; // SiPM 부착 여부와 상관없이 Z축은 계속 전진 (Vacuum 포함)
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
    std::vector<G4LogicalVolume*> fiberUnitIntersection__[],
    std::vector<G4LogicalVolume*> fiberCladIntersection__[],
    std::vector<G4LogicalVolume*> fiberCoreIntersection__[]) {

    if (!doFiber) return;

    // 1. 치수 및 원본 솔리드 정의
    G4double pitch_horizontal = 1.35 * mm; 
    G4double pitch_vertical = 1.22 * mm;   
    G4int fiberId = 0;
    G4double longFiberLen = 1000.0 * mm;

    // 반지름 설정
    G4double clad2_rMax = 0.553 * mm; // 요청하신 수치
    G4double clad_rMax  = 0.500 * mm; 
    G4double core_rMax  = 0.485 * mm;

    // 원본 튜브 정의
    G4Tubs* fiberClad2_Long = new G4Tubs("fClad2_L", 0, clad2_rMax, longFiberLen / 2.0, 0, 360*deg);
    G4Tubs* fiberClad_Long  = new G4Tubs("fClad_L",  0, clad_rMax,  longFiberLen / 2.0, 0, 360*deg);
    G4Tubs* fiberCore_Long  = new G4Tubs("fCore_L",  0, core_rMax,  longFiberLen / 2.0, 0, 360*deg);

    G4RotationMatrix* fiberRot = new G4RotationMatrix();
    fiberRot->rotateX(90*deg);

    // 시각화 속성
    G4VisAttributes* cladVis = new G4VisAttributes(G4Colour(0.0, 0.0, 0.0, 0.1)); // 검은색
    cladVis->SetForceSolid(true);
    G4VisAttributes* coreVis = new G4VisAttributes(G4Colour(1.0, 0.5, 0.0,0.7));      // 주황색
    coreVis->SetForceSolid(true);
    coreVis->SetVisibility(true);
    G4VisAttributes* glueVis = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.3)); 
    glueVis->SetForceSolid(true);
    glueVis->SetVisibility(true);

    // 2. 람다 함수 (레이어별 배치)
    auto FillFibersInMother = [&](G4LogicalVolume* motherLog) {
        G4Trd* motherSolid = (G4Trd*)motherLog->GetSolid();
        G4double h1 = motherSolid->GetXHalfLength1();
        G4double h2 = motherSolid->GetXHalfLength2();
        G4double zHalf = motherSolid->GetZHalfLength();
        G4double thickness = zHalf * 2.0;

        int numRows = std::floor(thickness / pitch_vertical);
        G4double startV = -((numRows - 1) * pitch_vertical) / 2.0;

        for (int k = 0; k < numRows; k++) {
            G4double localV = startV + k * pitch_vertical; 
            G4double currentHalfWidth = h1 + (h2 - h1) * (localV + zHalf) / (2.0 * zHalf);
            
            int numCols = std::floor((currentHalfWidth * 2.0 - 1.0*mm) / pitch_horizontal); 
            G4double startH = -((numCols - 1) * pitch_horizontal) / 2.0;

            for (int j = 0; j < numCols; j++) {
                G4double localH = startH + j * pitch_horizontal;
                if (k % 2 != 0) localH += pitch_horizontal / 2.0;
                if (std::abs(localH) + clad2_rMax > currentHalfWidth - 0.5*mm) continue;

                G4ThreeVector fiberPos(localH, 0, localV);
                G4Transform3D transform(*fiberRot, fiberPos);

                // --- 1단계: Clad2 (납을 뚫는 Galactic 구멍) ---
                G4IntersectionSolid* intClad2 = new G4IntersectionSolid("Clad2_Int", motherSolid, fiberClad2_Long, transform);
                G4LogicalVolume* clad2Log = new G4LogicalVolume(intClad2, FindMaterial("G4_Galactic"), "Clad2_Log");
                new G4PVPlacement(0, G4ThreeVector(0,0,0), clad2Log, "Clad2_Phys", motherLog, false, fiberId, false);

                clad2Log->SetVisAttributes(glueVis);

                // --- 2단계: Clad (PMMA 클래딩) ---
                G4IntersectionSolid* intClad = new G4IntersectionSolid("Clad_Int", motherSolid, fiberClad_Long, transform);
                G4LogicalVolume* cladLog = new G4LogicalVolume(intClad, FindMaterial("PMMA"), "Clad_Log");
                cladLog->SetVisAttributes(cladVis);
                new G4PVPlacement(0, G4ThreeVector(0,0,0), cladLog, "Clad_Phys", clad2Log, false, fiberId, false);

                // --- 3단계: Core (Polystyrene 코어) ---
                G4IntersectionSolid* intCore = new G4IntersectionSolid("Core_Int", motherSolid, fiberCore_Long, transform);
                G4LogicalVolume* coreLog = new G4LogicalVolume(intCore, FindMaterial("Polystyrene"), "Core_Log");
                coreLog->SetVisAttributes(coreVis);
                new G4PVPlacement(0, G4ThreeVector(0,0,0), coreLog, "Core_Phys", cladLog, false, fiberId, false);

                // 리전 및 데이터 저장
 //               fScintRegion->AddRootLogicalVolume(cladLog);
 //               fScintRegion->AddRootLogicalVolume(coreLog);
 //               cladLog->SetRegion(fScintRegion);
 //               coreLog->SetRegion(fScintRegion);

                fiberCladIntersection__[i].push_back(cladLog);
                fiberCoreIntersection__[i].push_back(coreLog);

                fiberId++;
            }
        }
    };

    for (auto pbLog : logicPbSFIL) FillFibersInMother(pbLog);
    for (auto pbLog : logicPbBulk) FillFibersInMother(pbLog);
}

