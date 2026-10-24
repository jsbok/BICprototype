#include "koBICMaterials.hh"
#include "G4SystemOfUnits.hh"

#include <stdio.h>
#include <algorithm>
#include <cmath>

koBICMaterials* koBICMaterials::fInstance = 0;

koBICMaterials::koBICMaterials() {
  fNistMan = G4NistManager::Instance();
  CreateMaterials();
}

koBICMaterials::~koBICMaterials() {}

koBICMaterials* koBICMaterials::GetInstance() {
  if (fInstance==0) fInstance = new koBICMaterials();

  return fInstance;
}

G4Material* koBICMaterials::GetMaterial(const G4String matName) {
  G4Material* mat = fNistMan->FindOrBuildMaterial(matName);

  if (!mat) mat = G4Material::GetMaterial(matName);
  if (!mat) {
    std::ostringstream o;
    o << "Material " << matName << " not found!";
    G4Exception("koBICMaterials::GetMaterial","",FatalException,o.str().c_str());
  }

  return mat;
}

G4OpticalSurface* koBICMaterials::GetOpticalSurface(const G4String surfName) {
  if (surfName=="SiPMSurf") return fSiPMSurf;
  else if (surfName=="FilterSurf") return fFilterSurf;
  else if (surfName=="MirrorSurf") return fMirrorSurf;
  else if (surfName=="AlSurf") return fAlSurf;
  else {
    std::ostringstream o;
    o << "OpticalSurface " << surfName << " not found!";
    G4Exception("koBICMaterials::GetOpticalSurface","",FatalException,o.str().c_str());
  }

  return nullptr;
}

void koBICMaterials::CreateMaterials() {
  fNistMan->FindOrBuildMaterial("G4_Galactic");
  fNistMan->FindOrBuildMaterial("G4_AIR");

  G4String symbol;
  G4double a, z, density;
  G4int ncomponents, natoms;
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z=1., a=1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z=6., a=12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z=7., a=14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  ,symbol="O" , z=8., a=16.00*g/mole);
  G4Element* F  = new G4Element("Fluorine",symbol="F" , z=9., a=18.9984*g/mole);
  G4Element* Si = new G4Element("Silicon", "Si", 14., 28.0855*g/mole);
  G4Element* Cl = new G4Element("Chlorine", "Cl", 17., 35.453*g / mole);

  fCu = new G4Material("Copper"  , z = 29., a = 63.546 * g/mole, density = 8.96  * g/cm3);
  fW  = new G4Material("Tungsten", z = 74., a = 183.84 * g/mole, density = 19.30 * g/cm3);
  fFe = new G4Material("Iron"    , z = 26., a = 55.845 * g/mole, density = 7.874 * g/cm3);
  fPb = new G4Material("Lead"    , z = 82., a = 207.2  * g/mole, density = 11.35 * g/cm3);
  fSi = new G4Material("Silicon" , z = 14., a = 28.09  * g/mole, density = 2.33  * g/cm3);
  fAl = new G4Material("Aluminum", z = 13., a = 26.98  * g/mole, density = 2.699 * g/cm3);

  fVacuum = G4Material::GetMaterial("G4_Galactic");
  fAir = G4Material::GetMaterial("G4_AIR");

  fFluoPoly = new G4Material("FluorinatedPolymer", density=1.43*g/cm3, ncomponents=2);
  fFluoPoly->AddElement(C, 2);
  fFluoPoly->AddElement(F, 2);

  fGlass = new G4Material("Glass", density=1.032*g/cm3, 2);
  fGlass->AddElement(C, 91.533*perCent);
  fGlass->AddElement(H, 8.467*perCent);

  fPS = new G4Material("Polystyrene", density=1.05*g/cm3, ncomponents=2);
  fPS->AddElement(C, natoms=8);
  fPS->AddElement(H, natoms=8);

  fPMMA = new G4Material("PMMA", density= 1.19*g/cm3, ncomponents=3);
  fPMMA->AddElement(C, natoms=5);
  fPMMA->AddElement(H, natoms=8);
  fPMMA->AddElement(O, natoms=2);

  //Acrylic at 20oC
  fAcrylic = new G4Material("Acrylic", 1.190 * g/cm3, 3);
  fAcrylic->AddElement(C, 5);
  fAcrylic->AddElement(H, 8);
  fAcrylic->AddElement(O, 2);

  //Optical cookie
  fEJ560 = new G4Material("EJ560", 1.03 * g/cm3, 2);
  fEJ560->AddElement(C, 2);
  fEJ560->AddElement(Si, 1);

  //MPPC Window Epoxy
  fEpoxyResin = new G4Material("EpoxyResin", 2.25 * g / cm3, 4);
  fEpoxyResin->AddElement(C, 21);
  fEpoxyResin->AddElement(H, 25);
  fEpoxyResin->AddElement(O, 5);
  fEpoxyResin->AddElement(Cl, 1);

  fGelatin = new G4Material("Gelatin", density=1.27*g/cm3, ncomponents=4);
  fGelatin->AddElement(C, natoms=102);
  fGelatin->AddElement(H, natoms=151);
  fGelatin->AddElement(N, natoms=31);
  fGelatin->AddElement(O, natoms=39);

  fGlue = new G4Material("Glue", density=1.18*g/cm3, ncomponents=4);
  fGlue->AddElement(C, 0.822);  // Carbon 82.2%
  fGlue->AddElement(H, 0.091);  // Hydrogen 9.1%
  fGlue->AddElement(O, 0.055);  // Oxygen 5.5%
  fGlue->AddElement(N, 0.032);  // Nitrogen 3.2%

  G4MaterialPropertiesTable* mpVacuum;
  G4MaterialPropertiesTable* mpAir;
  G4MaterialPropertiesTable* mpPS;
  G4MaterialPropertiesTable* mpPMMA;
  G4MaterialPropertiesTable* mpAcrylic;
  G4MaterialPropertiesTable* mpFluoPoly;
  G4MaterialPropertiesTable* mpGlass;
  G4MaterialPropertiesTable* mpAl;
  G4MaterialPropertiesTable* mpAlSurf;
  G4MaterialPropertiesTable* mpEJ560;
  G4MaterialPropertiesTable* mpGlue;
  G4MaterialPropertiesTable* mpEpoxyResin;
  G4MaterialPropertiesTable* mpSiPM;
  G4MaterialPropertiesTable* mpFilter;
  G4MaterialPropertiesTable* mpFilterSurf;
  G4MaterialPropertiesTable* mpMirror;
  G4MaterialPropertiesTable* mpMirrorSurf;

  G4double opEn[] = { // from 900nm to 300nm with 25nm step
    1.37760*eV, 1.41696*eV, 1.45864*eV, 1.50284*eV, 1.54980*eV, 1.59980*eV, 1.65312*eV, 1.71013*eV,
    1.77120*eV, 1.83680*eV, 1.90745*eV, 1.98375*eV, 2.06640*eV, 2.15625*eV, 2.25426*eV, 2.36160*eV,
    2.47968*eV, 2.61019*eV, 2.75520*eV, 2.91728*eV, 3.09960*eV, 3.30625*eV, 3.54241*eV, 3.81490*eV, 4.13281*eV
  };

  const G4int nEnt = sizeof(opEn) / sizeof(G4double);

  G4double RI_Air[nEnt]; std::fill_n(RI_Air,nEnt,1.0);
  mpAir = new G4MaterialPropertiesTable();
  mpAir->AddProperty("RINDEX",opEn,RI_Air,nEnt);
  fAir->SetMaterialPropertiesTable(mpAir);
/*
  G4double RI_Vacuum[nEnt]; std::fill_n(RI_Vacuum,nEnt,1.0);
  mpVacuum = new G4MaterialPropertiesTable();
  mpVacuum->AddProperty("RINDEX",opEn,RI_Vacuum,nEnt);
  fVacuum->SetMaterialPropertiesTable(mpVacuum);
*/
  G4double RI_PMMA[nEnt] = {
    1.48329, 1.48355, 1.48392, 1.48434, 1.48467, 1.48515, 1.48569, 1.48628,
    1.48677, 1.48749, 1.48831, 1.48899, 1.49000, 1.49119, 1.49219, 1.49372,
    1.49552, 1.49766, 1.49953, 1.50252, 1.50519, 1.51000, 1.51518, 1.52182, 1.53055
  };
  G4double AbsLen_PMMA[nEnt] = {
    0.414*m, 0.543*m, 0.965*m, 2.171*m, 2.171*m, 3.341*m, 4.343*m, 1.448*m,
    4.343*m, 14.48*m, 21.71*m, 8.686*m, 28.95*m, 54.29*m, 43.43*m, 48.25*m,
    54.29*m, 48.25*m, 43.43*m, 28.95*m, 21.71*m, 4.343*m, 2.171*m, 0.869*m, 0.434*m
  };
  mpPMMA = new G4MaterialPropertiesTable();
  mpPMMA->AddProperty("RINDEX",opEn,RI_PMMA,nEnt);
  mpPMMA->AddProperty("ABSLENGTH",opEn,AbsLen_PMMA,nEnt);
  fPMMA->SetMaterialPropertiesTable(mpPMMA);

  G4double RI_Acrylic[nEnt] = {
    1.48329, 1.48355, 1.48392, 1.48434, 1.48467, 1.48515, 1.48569, 1.48628,
    1.48677, 1.48749, 1.48831, 1.48899, 1.49000, 1.49119, 1.49219, 1.49372,
    1.49552, 1.49766, 1.49953, 1.50252, 1.50519, 1.51000, 1.51518, 1.52182, 1.53055
  };
  G4double AbsLen_Acrylic[nEnt] = {
    0.414*m, 0.543*m, 0.965*m, 2.171*m, 2.171*m, 3.341*m, 4.343*m, 1.448*m,
    4.343*m, 14.48*m, 21.71*m, 8.686*m, 28.95*m, 54.29*m, 43.43*m, 48.25*m,
    54.29*m, 48.25*m, 43.43*m, 28.95*m, 21.71*m, 4.343*m, 2.171*m, 0.869*m, 0.434*m
  };
  //G4double rindexAcrylic[nEnt]; std::fill_n(rindexAcrylic,nEnt,1.502);  // RI at 436nm
  //G4double AbsAcrylic[nEnt]; std::fill_n(AbsAcrylic,nEnt,2.4*m);
  mpAcrylic = new G4MaterialPropertiesTable();
  mpAcrylic->AddProperty("RINDEX", opEn, RI_Acrylic, nEnt);
  mpAcrylic->AddProperty("ABSLENGTH", opEn, AbsLen_Acrylic, nEnt);
  fAcrylic->SetMaterialPropertiesTable(mpAcrylic);

  G4double RI_FluoPoly[nEnt]; std::fill_n(RI_FluoPoly, nEnt, 1.42);
  mpFluoPoly = new G4MaterialPropertiesTable();
  mpFluoPoly->AddProperty("RINDEX",opEn,RI_FluoPoly,nEnt);
  fFluoPoly->SetMaterialPropertiesTable(mpFluoPoly);

  G4double RI_PS[nEnt] = {
    1.57483, 1.57568, 1.57644, 1.57726, 1.57817, 1.57916, 1.58026, 1.58148,
    1.58284, 1.58435, 1.58605, 1.58796, 1.59013, 1.59328, 1.59621, 1.59960,
    1.60251, 1.60824, 1.61229, 1.62032, 1.62858, 1.63886, 1.65191, 1.66888, 1.69165
  };
  G4double AbsLen_PS[nEnt] = {
    2.714*m, 3.102*m, 3.619*m, 4.343*m, 5.791*m, 7.896*m, 4.343*m, 7.896*m,
    5.429*m, 36.19*m, 17.37*m, 36.19*m, 5.429*m, 28.95*m, 21.71*m, 14.48*m,
    12.41*m, 8.686*m, 7.238*m, 1.200*m, 0.200*m, 0.500*m, 0.200*m, 0.100*m, 0.100*m
  };
  G4double scintFast_PS[nEnt] = {
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.00, 0.03, 0.07, 0.13,
    0.33, 0.63, 1.00, 0.50, 0.00, 0.00, 0.00, 0.00, 0.00
  };
  mpPS = new G4MaterialPropertiesTable();
  mpPS->AddProperty("RINDEX",opEn,RI_PS,nEnt);
  mpPS->AddProperty("ABSLENGTH",opEn,AbsLen_PS,nEnt);
  mpPS->AddConstProperty("SCINTILLATIONYIELD",10./keV);
  mpPS->AddConstProperty("RESOLUTIONSCALE",1);
  
#ifdef G4_LATEST
//  mpPS->AddProperty("FASTCOMPONENT",opEn,scintFast_PS,nEnt,true);
//  mpPS->AddConstProperty("FASTTIMECONSTANT",2.8*ns,true);

  mpPS->AddProperty("SCINTILLATIONCOMPONENT1", opEn, scintFast_PS, nEnt, true);
  mpPS->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.8*ns, true);
  mpPS->AddConstProperty("SCINTILLATIONBYPARTICLETYPE", false, true);
#else
  mpPS->AddProperty("FASTCOMPONENT",opEn,scintFast_PS,nEnt);
  mpPS->AddConstProperty("FASTTIMECONSTANT",2.8*ns);
  mpPS->AddConstProperty("SCINTILLATIONBYPARTICLETYPE", false);
#endif
  fPS->SetMaterialPropertiesTable(mpPS);
//  fPS->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
       
  G4double RI_Glass[nEnt]; std::fill_n(RI_Glass, nEnt, 1.52);
  G4double Abslength_Glass[nEnt]; std::fill_n(Abslength_Glass, nEnt, 420.*cm);
  mpGlass = new G4MaterialPropertiesTable();
  mpGlass->AddProperty("RINDEX",opEn,RI_Glass,nEnt);
  mpGlass->AddProperty("ABSLENGTH",opEn,Abslength_Glass,nEnt);
  fGlass->SetMaterialPropertiesTable(mpGlass);

  G4double refl_SiPM[nEnt]; std::fill_n(refl_SiPM, nEnt, 0.);
  G4double RI_EJ560[nEnt]; std::fill_n(RI_EJ560,nEnt,1.43); // at 436nm
  G4double AbsEJ560[nEnt]; std::fill_n(AbsEJ560,nEnt,2.4*m);
  mpEJ560 = new G4MaterialPropertiesTable();
  mpEJ560->AddProperty("RINDEX", opEn, RI_EJ560, nEnt);
  mpEJ560->AddProperty("ABSLENGTH", opEn, AbsEJ560, nEnt);
  fEJ560->SetMaterialPropertiesTable(mpEJ560);

  //G4double reflectivityEpoxyResin[2] = {0.99, 0.99};
  //G4double Abs_EpoxyResin[2] = {5.*m, 5.*m};
  G4double RI_EpoxyResin[nEnt]; std::fill_n(RI_EpoxyResin,nEnt,1.57);
  mpEpoxyResin = new G4MaterialPropertiesTable();
  mpEpoxyResin->AddProperty("RINDEX", opEn, RI_EpoxyResin, nEnt);
  //EpoxyResin_mpt->AddProperty("ABSLENGTH", opEn, Abs_EpoxyResin, 2);
  fEpoxyResin->SetMaterialPropertiesTable(mpEpoxyResin);

/*
G4double eff_SiPM[nEnt] = {
    1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
    1.00, 1.00, 1.00, 1.00, 1.00, 1.0, 1.00, 1.0,
    1.00, 1.00, 1.00, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
  };

G4double eff_SiPM[nEnt] = {
    0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
    0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001,
    0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001
  };
*/

// Bad
/*
G4double eff_SiPM[nEnt] = {
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.01, 0.01, 0.02, 0.03,
    0.05, 0.06, 0.08, 0.09, 0.10, 0.10, 0.10, 0.10, 0.09
  };*/
/*
  G4double eff_SiPM[nEnt] = {  // from 900nm to 300nm with 25nm step
    0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.01, 0.02, 0.03, 0.07,
    0.11, 0.13, 0.17, 0.19, 0.21, 0.21, 0.21, 0.20, 0.19
  };
*/
/*
//real PMT
   G4double eff_SiPM[nEnt] = {
0,0,0,0,0,0,0,0,0,0,0,0.01,0.02,0.04,0.06,0.12,0.20,0.23,0.28,0.32, 0.35, 0.35, 0.35, 0.32, 0.3
};
*/
// =========================================================
  // New SiPM PDE Data (from 900nm to 300nm)
  // Matching the existing opEn array (1.3776 eV to 4.13281 eV)
  // =========================================================
  G4double eff_SiPM[nEnt] = {
    0.040435, 0.050870, 0.063913, 0.078261, 0.095217, 
    0.109565, 0.129130, 0.147391, 0.170870, 0.195652, 
    0.226957, 0.259565, 0.296087, 0.335217, 0.371739, 
    0.414783, 0.457826, 0.498261, 0.508696, 0.499565, 
    0.473478, 0.420000, 0.392609, 0.354783, 0.238696
  };

/*
  //S14 QE
  G4double eff_SiPM[nEnt] = {
    0.04, 0.05, 0.07, 0.10, 0.12, 0.14, 0.15, 0.16,
    0.17, 0.20, 0.22, 0.27, 0.30, 0.32, 0.37, 0.42,
    0.45, 0.50, 0.51, 0.50, 0.46, 0.43, 0.40, 0.35, 0.20
  };

  //S13 QE
  G4double eff_SiPM[nEnt] = {
    0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.11, 0.14,
    0.15, 0.17, 0.20, 0.23, 0.26, 0.30, 0.33, 0.35,
    0.38, 0.39, 0.40, 0.37, 0.34, 0.30, 0.20, 0.05, 0.00
  };
*/

  mpSiPM = new G4MaterialPropertiesTable();
  mpSiPM->AddProperty("REFLECTIVITY",opEn,refl_SiPM,nEnt);
  mpSiPM->AddProperty("EFFICIENCY",opEn,eff_SiPM,nEnt);
  fSiPMSurf = new G4OpticalSurface("SiPMSurf",glisur,polished,dielectric_metal);
  fSiPMSurf->SetMaterialPropertiesTable(mpSiPM);
/*
  G4double filterEff[nEnt] = {
    1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
    1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.00, 1.,
    1., 1., 1., 1.00, 1.000, 1.000, 1.000, 1.000, 1.000
  };
*/

  G4double filterEff[nEnt] = {
    1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000,
    1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.900, 0.734,
    0.568, 0.402, 0.296, 0.070, 0.000, 0.000, 0.000, 0.000, 0.000
  };

  G4double filterRef[nEnt]; std::fill_n(filterRef,nEnt,0.);
  G4double RI_gel[nEnt]; std::fill_n(RI_gel,nEnt,1.52);
  mpFilter = new G4MaterialPropertiesTable();
  mpFilter->AddProperty("RINDEX",opEn,RI_gel,nEnt);
  fGelatin->SetMaterialPropertiesTable(mpFilter);
  mpFilterSurf = new G4MaterialPropertiesTable();
  mpFilterSurf->AddProperty("TRANSMITTANCE",opEn,filterEff,nEnt);
  mpFilterSurf->AddProperty("REFLECTIVITY",opEn,filterRef,nEnt);
  fFilterSurf = new G4OpticalSurface("FilterSurf",glisur,polished,dielectric_dielectric);
  fFilterSurf->SetMaterialPropertiesTable(mpFilterSurf);

  G4double MirrorRef[nEnt]; std::fill_n(MirrorRef, nEnt, 0.9);
  G4double MirrorEff[nEnt]; std::fill_n(MirrorEff, nEnt, 0.);

  mpMirrorSurf = new G4MaterialPropertiesTable();
  mpMirrorSurf->AddProperty("TRANSMITTANCE",opEn,MirrorEff,nEnt);
  mpMirrorSurf->AddProperty("REFLECTIVITY",opEn,MirrorRef,nEnt);
  fMirrorSurf = new G4OpticalSurface("MirrorSurf",glisur,polished,dielectric_metal);
  fMirrorSurf->SetMaterialPropertiesTable(mpMirrorSurf);

  G4double RI_glue[nEnt]; std::fill_n(RI_glue,nEnt,1.5);
  G4double Ab_glue[nEnt]; std::fill_n(Ab_glue,nEnt,5.0*m);

  mpGlue = new G4MaterialPropertiesTable();
  mpGlue->AddProperty("RINDEX",    opEn, RI_glue,  nEnt);
  mpGlue->AddProperty("ABSLENGTH", opEn, Ab_glue, nEnt);
  fGlue->SetMaterialPropertiesTable(mpGlue);

  //Alumiminium Mylar
  G4double RI_Al[nEnt]; std::fill_n(RI_Al,nEnt,1.2);
  G4double AbsAl[nEnt]; std::fill_n(AbsAl,nEnt, 64.8 * cm);

  G4double ReflectivityAl[nEnt]; std::fill_n(ReflectivityAl,nEnt,0.96); //Al
  //G4double ReflectivityAl[nEnt]; std::fill_n(ReflectivityAl,nEnt,0.0); //Perfect Blacksheet

  G4double EfficiencyAl[nEnt]; std::fill_n(EfficiencyAl,nEnt, 0.04);
  mpAl = new G4MaterialPropertiesTable();
  mpAl->AddProperty("RINDEX", opEn, RI_Al, nEnt);
  //mpAl->AddProperty("ABSLENGTH", opEn, AbsAl, 2);
  mpAl->AddProperty("REFLECTIVITY", opEn, ReflectivityAl, nEnt);
  mpAl->AddProperty("EFFICIENCY", opEn, EfficiencyAl, nEnt);
  fAl->SetMaterialPropertiesTable(mpAl);
  fAlSurf = new G4OpticalSurface("AlSurf", unified, polishedbackpainted, dielectric_dielectric);
  fAlSurf->SetMaterialPropertiesTable(mpAl);

}
