#include "DRsimCellParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

DRsimCellParameterisation::DRsimCellParameterisation(const G4int numx, const G4int numy, const G4double moduleH, const G4double moduleW)
: G4VPVParameterisation()
{
  for (G4int k = 0; k < numx; k++ ) {
    for (G4int j = 0; j < numy; j++ ) {

      if (k%2!=0 && j==numy-1) break;
      fXCell.push_back( -moduleH*mm/2 + k*1.22*mm + 0.61*mm );
      if (k%2==0) {fYCell.push_back( -moduleW*mm/2 + j*1.35*mm + 0.675*mm );}
      else {fYCell.push_back( -moduleW*mm/2 + j*1.35*mm + 1.35*mm );}
    }
  }
}

DRsimCellParameterisation::~DRsimCellParameterisation()
{}

void DRsimCellParameterisation::ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const {
  physVol->SetTranslation(G4ThreeVector(fXCell[copyNo],fYCell[copyNo],0.));
}
