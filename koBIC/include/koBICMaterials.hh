#ifndef koBICMaterials_h
#define koBICMaterials_h 1

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4OpticalSurface.hh"

class koBICMaterials {
public:
  virtual ~koBICMaterials();
  static koBICMaterials* GetInstance();
  G4Material* GetMaterial(const G4String);
  G4OpticalSurface* GetOpticalSurface(const G4String);

private:
  koBICMaterials();
  void CreateMaterials();

  static koBICMaterials* fInstance;
  G4NistManager* fNistMan;
  G4Material* fGlue;
  G4Material* fVacuum;
  G4Material* fAir;
  G4Material* fFluoPoly;
  G4Material* fPMMA;
  G4Material* fPS;
  G4Material* fCu;
  G4Material* fW;
  G4Material* fPb;
  G4Material* fFe;
  G4Material* fSi;
  G4Material* fAl;
  G4Material* fGlass;
  G4Material* fGelatin;
  G4Material* fMirror;
  G4Material* fBlackpaint;
  G4OpticalSurface* fSiPMSurf;
  G4OpticalSurface* fFilterSurf;
  G4OpticalSurface* fMirrorSurf;
};

#endif
