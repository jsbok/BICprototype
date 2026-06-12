#ifndef koBICFilterParameterisation_h
#define koBICFilterParameterisation_h 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4VisAttributes.hh"
#include <vector>

class G4VPhysicalVolume;

class koBICFilterParameterisation : public G4VPVParameterisation {
public:
  koBICFilterParameterisation(const G4int numx, const G4int numy, const G4double moduleH, const G4double moduleW);
  virtual ~koBICFilterParameterisation();

  virtual void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;

private:
  std::vector<G4double> fXFilter;
  std::vector<G4double> fYFilter;
};

#endif
