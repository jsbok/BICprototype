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
  fModuleW = 700.0; // 가로 폭 (필요시 수정)
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

  doFiber = false;
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

G4cout << "Step 1" << G4endl;

  if (doPMT) {
    // fModuleProp의 첫 번째 데이터(인덱스 0)가 반드시 존재해야 합니다.
    if (fModuleProp.empty()) return; 

    // 왼쪽 SiPM (Index 0)
    koBICSiPMSD *SiPMSDmodule_left = new koBICSiPMSD(
        "Module0", "ModuleC0", 0, fModuleProp.at(0));
    
    // 오른쪽 SiPM (Index 1)
    koBICSiPMSD *SiPMSDmodule_right = new koBICSiPMSD(
        "Module1", "ModuleC1", 1, fModuleProp.at(0));

    SDman->AddNewDetector(SiPMSDmodule_left);
    SDman->AddNewDetector(SiPMSDmodule_right);

    // PMTcathLogical의 0번과 1번에 각각 할당
    if (PMTcathLogical[0]) PMTcathLogical[0]->SetSensitiveDetector(SiPMSDmodule_left);
    if (PMTcathLogical[1]) PMTcathLogical[1]->SetSensitiveDetector(SiPMSDmodule_right);
  }
G4cout << "Step 2" << G4endl;
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

    G4Trd* parentTrd = (G4Trd*)ModuleLogical_[0]->GetSolid();
    G4double pDx1 = parentTrd->GetXHalfLength1();
    G4double pDx2 = parentTrd->GetXHalfLength2();
    G4double pDy  = parentTrd->GetYHalfLength1();
    G4double pDz  = parentTrd->GetZHalfLength();

    G4RotationMatrix* rot180 = new G4RotationMatrix();
    rot180->rotateZ(180*deg); 

    G4double currentZ = -pDz;
    int totalSectionID = 0; // 0~49까지의 Readout ID

    // --- [Helper 1] 진공 레이어용 (기존 방식 유지, 5분할 안 함) ---
    auto PlaceVacLayer = [&](G4String name, G4double thick) {
        G4double zStart = currentZ + pDz;
        G4double zEnd = zStart + thick;
        G4double lDx1 = pDx1 + (pDx2 - pDx1) * (zStart / (2 * pDz));
        G4double lDx2 = pDx1 + (pDx2 - pDx1) * (zEnd / (2 * pDz));

        G4VSolid* vSol = new G4Trd(name+"_sol", lDx1, lDx2, pDy, pDy, thick/2.0);
        G4LogicalVolume* vLog = new G4LogicalVolume(vSol, vacMat, name+"_log");
        vLog->SetVisAttributes(G4VisAttributes::GetInvisible()); // 진공은 안 보이게

        new G4PVPlacement(0, G4ThreeVector(0, 0, currentZ + thick/2.0), 
                          vLog, name+"_phys", ModuleLogical_[0], false, 100, checkOverlaps);
        currentZ += thick;
    };

    // --- [Helper 2] 납 레이어용 (5분할 + 지그재그 + 빗각 일치) ---
auto PlaceDividedPbLayer = [&](G4String namePrefix, G4double thick, G4double customDy = -1.0) {
        G4double zStart = currentZ + pDz;
        G4double zEnd = zStart + thick;
        
        // 현재 레이어의 Z축 시작과 끝에서의 전체 너비(Half-width)
        G4double fullDx1 = pDx1 + (pDx2 - pDx1) * (zStart / (2 * pDz));
        G4double fullDx2 = pDx1 + (pDx2 - pDx1) * (zEnd / (2 * pDz));
        G4double currentDy = (customDy > 0) ? customDy : pDy;

        // [핵심 1] 사용자님의 정답: "큰변 - 짧은변 차이를 유지해라"
        G4double delta = fullDx2 - fullDx1; 
        
        // [핵심 2] 5분할된 조각의 평균 너비 (전체 너비의 1/5)
        G4double w_avg = (fullDx1 + fullDx2) / 10.0; // Half-width 기준이므로 10으로 나눔

        // [핵심 3] 작은 사다리꼴의 치수 계산 (기울기가 엄마 볼륨과 100% 동일해짐)
        G4double subDx1 = w_avg - delta / 2.0;
        G4double subDx2 = w_avg + delta / 2.0;

        for (int s = 0; s < 5; s++) {
            G4double finalDx1, finalDx2;
            
            // 교차 배치: 짝수번은 정방향, 홀수번은 180도 뒤집힌 상태(치수 교체)
            if (s % 2 == 0) { 
                finalDx1 = subDx1; 
                finalDx2 = subDx2;
            } else { 
                // 회전 행렬 없이 치수만 반대로 넣어주면 기하학적으로 완벽한 180도 뒤집힘입니다.
                finalDx1 = subDx2; 
                finalDx2 = subDx1;
            }

            G4String name = "Module" + std::to_string(totalSectionID);
            G4VSolid* subSol = new G4Trd(name+"_sol", finalDx1, finalDx2, currentDy, currentDy, thick/2.0);
            G4LogicalVolume* subLog = new G4LogicalVolume(subSol, pbMat, name+"_log");

            // 회색 납 시각화 유지
            G4VisAttributes* pbVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 1.0));
            pbVis->SetForceSolid(true);
            subLog->SetVisAttributes(pbVis);

            // [핵심 4] X축 배치 위치 (간격 없이 완벽하게 맞물리는 중심 좌표 계산)
            G4double posX = (2 * s - 4) * w_avg;

            // 회전 행렬(Rotation) 없이 그대로 밀어 넣습니다.
            new G4PVPlacement(0, G4ThreeVector(posX, 0, currentZ + thick/2.0),
                              subLog, name+"_phys", ModuleLogical_[0], false, totalSectionID, checkOverlaps);
            totalSectionID++;
        }
        currentZ += thick;
    };

    // ========================================================================
    // 3. 실제 적층 시작 (삭제되었던 반복문들 복구)
    // ========================================================================

    // [A] SFIL 구역: 3개 레이어 반복 (납은 5분할 배치)
    for(int r = 0; r < 3; r++) {
        PlaceVacLayer("Vac_SFIL_" + std::to_string(r), 17.0*mm);
        PlaceDividedPbLayer("Pb_SFIL", 21.73*mm, -1.0); // SFIL 높이(-1.0은 pDy 사용)
    }

    // SFIL 마지막 진공층
    PlaceVacLayer("Vac_SFIL_End", 17.0*mm);

    // [B] BULK 구역: 7개 레이어 반복 (납은 5분할 배치, 높이 30cm)
    G4double bulkHalfHeight = (300.0 * mm) / 2.0; 
    for(int r = 0; r < 7; r++) {
        // Bulk 구역은 진공층 없이 납 레이어만 7개 쌓는 것인지, 
        // 혹은 다른 물질이 섞이는지 확인이 필요하나 우선 요청하신 대로 7번 반복 배치합니다.
        PlaceDividedPbLayer("Pb_Bulk", 20.889*mm, bulkHalfHeight);
    }


    // 5. Fiber 배치 (Fiber도 Z축 방향으로 길게 뻗어야 함)
    FiberImplement(0, ModuleLogical_, fiberUnitIntersection_, fiberCladIntersection_, fiberCoreIntersection_);

    // =========================================================
    // 3. PMT 및 속성 설정
    // =========================================================
 /*   if (doPMT) {
        dimCalc->SetisModule(false);
        pmtg = new G4Box("PMTG", (fModuleH / 2.) * mm, (fModuleW / 2.) * mm, PMTT / 2. * mm);
        PMTGLogical_[0] = new G4LogicalVolume(pmtg, FindMaterial("G4_AIR"), moduleName);
        PMTGLogical_[1] = new G4LogicalVolume(pmtg, FindMaterial("G4_AIR"), moduleName);
        
        new G4PVPlacement(zRot, dimCalc->GetOrigin_PMTG(0), PMTGLogical_[0],
                          moduleName, worldLogical, false, 0, checkOverlaps);
        new G4PVPlacement(zRot, dimCalc->GetOrigin_PMTG(1), PMTGLogical_[1], 
                          moduleName, worldLogical, false, 0, checkOverlaps);
    }
*/
    koBICInterface::koBICModuleProperty ModulePropSingle;
    ModulePropSingle.towerXY = fTowerXY;
    ModulePropSingle.ModuleNum = 0;
    ModuleProp_.push_back(ModulePropSingle);
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
    G4float fX = -fModuleH * mm / 2 + k * fFiber_vert_dis * mm + fFiber_vert_dis / 2 * mm;

bool inVacuum = false;
    double relX = fX - (-fModuleH * mm / 2.0); 
    double fStart = relX - clad_S_rMax;
    double fEnd = relX + clad_S_rMax;

    double curr = 0;
    // (1) 앞선 3번의 반복: 진공(17) + SFIL(21.73)
    for (int r = 0; r < 3; r++) {
      if (fStart < (curr + 17.0) && fEnd > curr) { inVacuum = true; break; }
      curr += (17.0 + 21.730);
    }
    // (2) 마지막 4번째 진공층(17) 검사
    if (!inVacuum) {
      if (fStart < (curr + 17.0) && fEnd > curr) { inVacuum = true; }
    }

    if (inVacuum) {
      fWhich = !fWhich; 
      continue; 
    }
    // ==========================================

    G4double rMin = 904.485; // XML 기준 SFIL 시작점
    G4double limit_Y = (rMin + relX) * std::tan((360.*deg/48.0)/2.0);

    for (int j = 0; j < NofFiber; j++) {
      if (fWhich && j == NofFiber - 1) break;
      G4float fY = -fModuleW * mm / 2 + j * fFiber_hori_dis * mm + fFiber_hori_dis / 2 * mm;
      if (fWhich) fY += fFiber_hori_dis / 2 * mm;
      
      // 사다리꼴 범위를 벗어나면 생성 안 함
      if (std::abs(fY) > limit_Y) continue; 
      
      fFiberX.push_back(fX);
      fFiberY.push_back(fY);
      fFiberWhich.push_back(fWhich);
    }
    fWhich = !fWhich;
  }
  if (doFiber) {
    for (unsigned int fiberId = 0; fiberId < fFiberX.size(); fiberId++) {

      double currentRelX = fFiberX.at(fiberId) - (-fModuleH * mm / 2.0);

      G4VSolid* targetClad = (currentRelX <= 133.19) ? fiberClad_SFIL : fiberClad_Bulk;
      G4VSolid* targetCore = (currentRelX <= 133.19) ? fiberCore_SFIL : fiberCore_Bulk;
      G4VSolid* currentModuleSolid = ModuleLogical__[i]->GetSolid();

      tfiberCladIntersection = new G4IntersectionSolid(
          "fiberClad", targetClad, currentModuleSolid, 0,
          G4ThreeVector(-fFiberX.at(fiberId), -fFiberY.at(fiberId), 0.));
      fiberCladIntersection__[i].push_back(new G4LogicalVolume(
          tfiberCladIntersection, FindMaterial("PMMA"), name));
      new G4PVPlacement(
          0, G4ThreeVector(fFiberX.at(fiberId), fFiberY.at(fiberId), 0),
          fiberCladIntersection__[i].at(fiberId), name, ModuleLogical__[i],
          false, fiberId, checkOverlaps);

      tfiberCoreIntersection = new G4IntersectionSolid(
          "fiberCore", targetCore, currentModuleSolid, 0,
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
