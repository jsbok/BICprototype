#include "koBICMirrorParameterisation.hh"
#include "koBICCellParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

koBICMirrorParameterisation::koBICMirrorParameterisation(const G4int numx, const G4int numy)
: G4VPVParameterisation()
{
  for ( G4int copyNo = 0; copyNo < numx*numy; copyNo++ ) {

    G4int column = copyNo / numy;
    G4int row = copyNo % numy;
  }
  fNumx = numx;
  fNumy = numy;
}

koBICMirrorParameterisation::~koBICMirrorParameterisation() {}

void koBICMirrorParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const {
  physVol->SetTranslation(G4ThreeVector(fXMirror[copyNo],fYMirror[copyNo],0.));
}
