#ifndef koBICMirrorParameterisation_h
#define koBICMirrorParameterisation_h 1

#include "globals.hh"
#include "G4VPVParameterisation.hh"
#include "G4VisAttributes.hh"
#include <vector>

class G4VPhysicalVolume;

class koBICMirrorParameterisation : public G4VPVParameterisation {
public:
  koBICMirrorParameterisation(const G4int numx, const G4int numy);
  virtual ~koBICMirrorParameterisation();

  virtual void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;

private:
  std::vector<G4double> fXMirror;
  std::vector<G4double> fYMirror;
  G4int fNumx;
  G4int fNumy;
};

#endif
