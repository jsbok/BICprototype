#include "koBICLightguide.hh"
#include "G4SystemOfUnits.hh"
#include "G4Trd.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"

#include "G4Trap.hh"
#include "G4IntersectionSolid.hh"
#include "G4GenericTrap.hh"
#include "G4TwoVector.hh"

#include "G4VisAttributes.hh"
#include <cmath>

LG::LG(G4double length, G4Material* acrylic, G4Material* al, G4OpticalSurface* alSurface)
: LGLength(length), Acrylic(acrylic), Al(al), AlSurface(alSurface)
{
}

void LG::TrapParamAdjustment(const double (&in)[11], double (&out)[11])
{
  // get input
  double dy1 = in[3], dx1 = in[4], dx2 = in[5], alp1 = in[6];
  double dy2 = in[7], dx3 = in[8], dx4 = in[9], alp2 = in[10];

  // adjust alpha
  //alp1 = alp2 = (alp1 + alp2)/2.;
  alp1 = alp2 = std::atan((std::tan(alp1) + std::tan(alp2))/2.);

  // adjust dx
  double k1 = (dx1 - dx2)/(dy1 + dy1);
  double k2 = (dx3 - dx4)/(dy2 + dy2) ;
  double k  = (k1 + k2)/2.;
  double m1 = (dx1 + dx2)/2.;
  double m2 = (dx3 + dx4)/2.;
  dx1 = m1 + k*dy1; dx2 = m1 - k*dy1;
  dx3 = m2 + k*dy2; dx4 = m2 - k*dy2;

  // set output
  out[0] = in[0]; out[1] = in[1]; out[2] = in[2];
  out[3] = in[3]; out[4] = dx1; out[5] = dx2; out[6]  = alp1;
  out[7] = in[7]; out[8] = dx3; out[9] = dx4; out[10] = alp2;
}



void LG::Construct()
{
  //-------------Square-square Trapezoid LG----------------------
  G4double airGap = 0.1*mm;
  G4double mylarThickness = 0.008*mm;
  G4double sideThickness = std::cos(std::atan((8.5*mm)/LGLength));

  solidAcrylicTrd = new G4Trd("solidAcrylicTrd", 15.0*mm, 6.5*mm, 15.0*mm, 6.5*mm, LGLength/2.0);
  logicAcrylicTrd = new G4LogicalVolume(solidAcrylicTrd, Acrylic, "logicAcrylicTrd");

  G4double outerX1 = 15.1*mm + mylarThickness/sideThickness;
  G4double outerX2 = 6.6*mm + mylarThickness/sideThickness;
  G4double outerY1 = 15.1*mm + mylarThickness/sideThickness;
  G4double outerY2 = 6.6*mm + mylarThickness/sideThickness;

  solidAlTrd = new G4Trd("solidAlTrd", outerX1, outerX2, outerY1, outerY2, LGLength/2.0);
  solidAirTrd = new G4Trd("solidAirTrd", 15.1*mm, 6.6*mm, 15.1*mm, 6.6*mm, LGLength/2.0);
  solidAlMylarRaw = new G4SubtractionSolid("solidAlMylarRaw", solidAlTrd, solidAirTrd);

  solidAlMylar = new G4IntersectionSolid("solidAlMylar", solidAlMylarRaw, solidPMTGClip);

  logicAlMylar = new G4LogicalVolume(solidAlMylar, Al, "logicAlMylar");
  new G4LogicalSkinSurface("LGAlSkin", logicAlMylar, AlSurface);


  // Nominal, Hexagonal Layer 1 - 12 Lightguide =========================== 
  double in[11];
  double out[11];  


  if (lightguideShape == "b1" || lightguideShape == "bh1") {
    lightguideSizeAbottom = 20.937 * mm;
    lightguideSizeAtop = 23.786 * mm;
    lightguideSizeAheight = 21.730 * mm;
    lightguideSizeBbottom = 13.0 * mm;
    lightguideSizeBtop = 13.0 * mm;
    lightguideSizeBheight = 13.0 * mm;
  }
  else if (lightguideShape == "b2" || lightguideShape == "bh2") {
    lightguideSizeAbottom = 21.953 * mm;
    lightguideSizeAtop = 24.801 * mm;
    lightguideSizeAheight = 21.730 * mm;
    lightguideSizeBbottom = 13.0 * mm;
    lightguideSizeBtop = 13.0 * mm;
    lightguideSizeBheight = 13.0 * mm;
  }
  else if (lightguideShape == "b3" || lightguideShape == "bh3") {
    lightguideSizeAbottom = 22.968 * mm;
    lightguideSizeAtop = 25.817 * mm;
    lightguideSizeAheight = 21.730 * mm;
    lightguideSizeBbottom = 13.0 * mm;
    lightguideSizeBtop = 13.0 * mm;
    lightguideSizeBheight = 13.0 * mm;
  }
  else if (lightguideShape == "b4" || lightguideShape == "bh4") {
    lightguideSizeAbottom = 23.983 * mm;
    lightguideSizeAtop = 26.831 * mm;
    lightguideSizeAheight = 21.730 * mm;
    lightguideSizeBbottom = 13.0 * mm;
    lightguideSizeBtop = 13.0 * mm;
    lightguideSizeBheight = 13.0 * mm;
  }
  else if (lightguideShape == "b5" || lightguideShape == "bh5") {
    lightguideSizeAbottom = 24.998 * mm;
    lightguideSizeAtop = 27.847 * mm;
    lightguideSizeAheight = 21.730 * mm;
    lightguideSizeBbottom = 13.0 * mm;
    lightguideSizeBtop = 13.0 * mm;
    lightguideSizeBheight = 13.0 * mm;
  }
  else if (lightguideShape == "b6" || lightguideShape == "bh6") {
    lightguideSizeAbottom = 26.058 * mm;
    lightguideSizeAtop = 28.795 * mm;
    lightguideSizeAheight = 20.879 * mm;
    lightguideSizeBbottom = 13.0 * mm;
    lightguideSizeBtop = 13.0 * mm;
    lightguideSizeBheight = 13.0 * mm;
  }
  else if (lightguideShape == "b7" || lightguideShape == "bh7") {
    lightguideSizeAbottom = 26.606 * mm;
    lightguideSizeAtop = 29.342 * mm;
    lightguideSizeAheight = 20.879 * mm;
    lightguideSizeBbottom = 13.0 * mm;
    lightguideSizeBtop = 13.0 * mm;
    lightguideSizeBheight = 13.0 * mm;
  }
  else if (lightguideShape == "b8" || lightguideShape == "bh8") {
    lightguideSizeAbottom = 27.153 * mm;
    lightguideSizeAtop = 29.890 * mm;
    lightguideSizeAheight = 20.879 * mm;
    lightguideSizeBbottom = 13.0 * mm;
    lightguideSizeBtop = 13.0 * mm;
    lightguideSizeBheight = 13.0 * mm;
  }
  else if (lightguideShape == "b9" || lightguideShape == "bh9") {
    lightguideSizeAbottom = 27.700 * mm;
    lightguideSizeAtop = 30.437 * mm;
    lightguideSizeAheight = 20.879 * mm;
    lightguideSizeBbottom = 13.0 * mm;
    lightguideSizeBtop = 13.0 * mm;
    lightguideSizeBheight = 13.0 * mm;
  }
  else if (lightguideShape == "b10" || lightguideShape == "bh10") {
    lightguideSizeAbottom = 28.248 * mm;
    lightguideSizeAtop = 30.985 * mm;
    lightguideSizeAheight = 20.879 * mm;
    lightguideSizeBbottom = 13.0 * mm;
    lightguideSizeBtop = 13.0 * mm;
    lightguideSizeBheight = 13.0 * mm;
  }
  else if (lightguideShape == "b11" || lightguideShape == "bh11") {
    lightguideSizeAbottom = 28.796 * mm;
    lightguideSizeAtop = 31.532 * mm;
    lightguideSizeAheight = 20.879 * mm;
    lightguideSizeBbottom = 13.0 * mm;
    lightguideSizeBtop = 13.0 * mm;
    lightguideSizeBheight = 13.0 * mm;
  }
  else if (lightguideShape == "b12" || lightguideShape == "bh12") {
    lightguideSizeAbottom = 29.342 * mm;
    lightguideSizeAtop = 32.079 * mm;
    lightguideSizeAheight = 20.879 * mm;
    lightguideSizeBbottom = 13.0 * mm;
    lightguideSizeBtop = 13.0 * mm;
    lightguideSizeBheight = 13.0 * mm;
  }
  //Nominal
  if (lightguideShape == "b1" || lightguideShape == "b2" || lightguideShape == "b3" || lightguideShape == "b4" || lightguideShape == "b5" || lightguideShape == "b6" || lightguideShape == "b7" || lightguideShape == "b8" || lightguideShape == "b9" || lightguideShape == "b10" || lightguideShape == "b11" || lightguideShape == "b12" ) 
  {
    std::vector<G4TwoVector> vertices;
    vertices.push_back(G4TwoVector(-lightguideSizeAtop/2.,lightguideSizeAheight/2.));
    vertices.push_back(G4TwoVector(lightguideSizeAtop/2.,lightguideSizeAheight/2.));
    vertices.push_back(G4TwoVector(lightguideSizeAbottom/2.,-lightguideSizeAheight/2.));
    vertices.push_back(G4TwoVector(-lightguideSizeAbottom/2.,-lightguideSizeAheight/2.));
    vertices.push_back(G4TwoVector(-lightguideSizeBtop/2.,lightguideSizeBheight/2.));
    vertices.push_back(G4TwoVector(lightguideSizeBtop/2.,lightguideSizeBheight/2.));
    vertices.push_back(G4TwoVector(lightguideSizeBbottom/2.,-lightguideSizeBheight/2.));
    vertices.push_back(G4TwoVector(-lightguideSizeBbottom/2.,-lightguideSizeBheight/2.));
    G4GenericTrap *lightGuideSolid = new G4GenericTrap("lightGuideSolid", LGLength / 2., vertices);
  }

  //Hexagonal light guides
  else if (lightguideShape == "bh1" || lightguideShape == "bh2" || lightguideShape == "bh3" || lightguideShape == "bh4" || lightguideShape == "bh5" || lightguideShape == "bh6" || lightguideShape == "bh7" || lightguideShape == "bh8" || lightguideShape == "bh9" || lightguideShape == "bh10" || lightguideShape == "bh11" || lightguideShape == "bh12") 
  { 
    in[0] = LGLength / 2.;
    in[1] = 0;
    in[2] = 0;
    in[3] = lightguideSizeAheight / 2.;
    in[4] = lightguideSizeAbottom / 2.;
    in[5] = lightguideSizeAtop / 2.;
    in[6] = 0; 
    in[7] = lightguideSizeBheight / 2.;
    in[8] = lightguideSizeBbottom / 2.;
    in[9] = ((lightguideSizeAtop-lightguideSizeAbottom)/lightguideSizeAheight*lightguideSizeBheight + lightguideSizeBbottom) / 2.;
    in[10] = 0;

    TrapParamAdjustment(in, out);
    lightguideSizeAbottom = out[4]*2.;
    lightguideSizeAtop = out[5]*2.;
    lightguideSizeBbottom = out[8]*2.;
    lightguideSizeBtop = out[8]*2.; // Update to adjusted dimensions for particle gun generation

    G4double motherLength = LGLength + 3.0*mm + 0.15*mm;
    solidOptMother = new G4Trap("OptMotherSolid", motherLength/2.0, 0.0*deg, 0.0*deg, out[3], out[4], out[5], out[6], out[3], out[4], out[5], out[6]);

    G4Trap *lightGuideASolid = new G4Trap("lightGuideASolid", out[0], out[1], out[2], out[3], out[4], out[5], out[6], out[7], out[8], out[9], out[10]);
    G4Trap *lightGuideBSolid = new G4Trap("lightGuideBSolid", out[0], out[1], out[2], out[3], out[5], out[5], out[6], out[7], out[8], out[8], out[10]);

    G4RotationMatrix *noRotation = new G4RotationMatrix();
    solidLightGuide = new G4IntersectionSolid("lightGuideSolid", lightGuideASolid, lightGuideBSolid, noRotation, G4ThreeVector(0, 0, 0));

    logicLightGuide = new G4LogicalVolume(solidLightGuide, Acrylic, "logicLightGuide");

    // Set LG volume (for Stepping Action)
    //fLGVolume = lightGuideLogical;


    in[0] = LGLength / 2.;
    in[1] = 0;
    in[2] = 0;
    in[3] = lightguideSizeAheight / 2. + wrapThickness + 1 * airGap;
    in[4] = lightguideSizeAbottom / 2. + wrapThickness + 1 * airGap;
    in[5] = lightguideSizeAtop / 2. + wrapThickness + 1 * airGap;
    in[6] = 0;
    in[7] = lightguideSizeBheight / 2. + wrapThickness + 1 * airGap;
    in[8] = lightguideSizeBbottom / 2. + wrapThickness + 1 * airGap;
    in[9] = ((lightguideSizeAtop-lightguideSizeAbottom)/lightguideSizeAheight*lightguideSizeBheight + lightguideSizeBbottom) / 2. + wrapThickness + 1 * airGap;
    in[10] = 0;

    TrapParamAdjustment(in, out);
    G4Trap *lightGuideAlASolid = new G4Trap("lightGuideAlASolid", out[0], out[1], out[2], out[3], out[4], out[5], out[6], out[7], out[8], out[9], out[10]);
    G4Trap *lightGuideAlBSolid = new G4Trap("lightGuideAlBSolid", out[0], out[1], out[2], out[3], out[5], out[5], out[6], out[7], out[8], out[8], out[10]);
    solidLightGuideAl = new G4IntersectionSolid("lightGuideSAlolid", lightGuideAlASolid, lightGuideAlBSolid, noRotation, G4ThreeVector(0, 0, 0));

    in[0] = LGLength / 2.;
    in[1] = 0;
    in[2] = 0;
    in[3] = lightguideSizeAheight / 2. + 1 * airGap;
    in[4] = lightguideSizeAbottom / 2. + 1 * airGap;
    in[5] = lightguideSizeAtop / 2.    + 1 * airGap;
    in[6] = 0; 
    in[7] = lightguideSizeBheight / 2. + 1 * airGap;
    in[8] = lightguideSizeBbottom / 2. + 1 * airGap;
    in[9] = ((lightguideSizeAtop-lightguideSizeAbottom)/lightguideSizeAheight*lightguideSizeBheight + lightguideSizeBbottom) / 2. + 1 * airGap;
    in[10] = 0;
    TrapParamAdjustment(in, out);
    G4Trap *lightGuideAirASolid = new G4Trap("lightGuideAirASolid", out[0], out[1], out[2], out[3], out[4], out[5], out[6], out[7], out[8], out[9], out[10]);                                                   
    G4Trap *lightGuideAirBSolid = new G4Trap("lightGuideAirBSolid", out[0], out[1], out[2], out[3], out[5], out[5], out[6], out[7], out[8], out[8], out[10]);
    solidLightGuideAir = new G4IntersectionSolid("lightGuideAirSolid", lightGuideAirASolid, lightGuideAirBSolid, noRotation, G4ThreeVector(0, 0, 0));


    solidLightGuideWrap = new G4SubtractionSolid("lightGuideWrapSolid", solidLightGuideAl, solidLightGuideAir, noRotation, G4ThreeVector(0, 0, 0));
    logicLightGuideWrap = new G4LogicalVolume(solidLightGuideWrap, Al, "logicLightGuideWrap");
    
    new G4LogicalSkinSurface("BH_LGAlSkin", logicLightGuideWrap, AlSurface);

    G4VisAttributes* wrapLogicalVisAttr = new G4VisAttributes();
    logicLightGuideWrap->SetVisAttributes(wrapLogicalVisAttr);
    wrapLogicalVisAttr->SetForceWireframe();

  }
}
