#ifndef koBICLightguide_h
#define koBICLightguide_h 1

#include "globals.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SystemOfUnits.hh"

class G4Material;
//class G4LogicalVolume;
//class G4VSolid;
class G4OpticalSurface;


class LG {
  public:
    LG(G4double length, G4Material* acrylic, G4Material* al, G4OpticalSurface* alSurface);
    ~LG() = default;

    void SetShape(const G4String& shape) { lightguideShape = shape; }
    void Construct();

    G4double GetLength() const { return LGLength; }

    // Square-square trapezoid LG
    void SetClipSolid(G4VSolid* solid) { solidPMTGClip = solid; }
    G4LogicalVolume* GetSquareAcrylicLogical() const { return logicAcrylicTrd; }
    G4LogicalVolume* GetSquareMylarLogical() const { return logicAlMylar; }

    // b/bh LG
    G4VSolid* GetMotherSolid() const { return solidOptMother; }
    G4LogicalVolume* GetBHAcrylicLogical() const { return logicLightGuide; }
    G4LogicalVolume* GetBHMylarLogical() const { return logicLightGuideWrap; }

    G4VSolid* GetMotherSolid(G4double motherLength);


  private:

    void TrapParamAdjustment(const double (&in)[11], double (&out)[11]);

    G4double LGLength;

    G4Material* Acrylic = nullptr;
    G4Material* Al = nullptr;
    G4OpticalSurface* AlSurface = nullptr;

    G4String lightguideShape;

    // b1-b12 / bh1-bh12 dimensions
    G4double lightguideSizeAbottom = 0.;
    G4double lightguideSizeAtop = 0.;
    G4double lightguideSizeAheight = 0.;
    G4double lightguideSizeBbottom = 0.;
    G4double lightguideSizeBtop = 0.;
    G4double lightguideSizeBheight = 0.;
    G4double wrapThickness = 0.008*mm;


    G4VSolid* solidPMTGClip = nullptr;
    G4VSolid* solidAcrylicTrd = nullptr;
    G4VSolid* solidAlTrd = nullptr;
    G4VSolid* solidAirTrd = nullptr;
    G4VSolid* solidLightGuideAl = nullptr;
    G4VSolid* solidLightGuideAir = nullptr;

    G4SubtractionSolid* solidAlMylarRaw = nullptr;
    G4IntersectionSolid* solidAlMylar = nullptr;

    G4LogicalVolume* logicAcrylicTrd = nullptr;
    G4LogicalVolume* logicAlMylar = nullptr;

    G4VSolid* solidOptMother = nullptr;
    G4VSolid* solidLightGuide = nullptr;
    G4VSolid* solidLightGuideWrap = nullptr;

    G4LogicalVolume* logicLightGuide = nullptr;
    G4LogicalVolume* logicLightGuideWrap = nullptr;

};

#endif
