#include "dimensionCalc.hh"

#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "geomdefs.hh"

#include <cmath>
#include <stdio.h>
#include <float.h>

using namespace std;

dimensionCalc::dimensionCalc() {

  fFrontL       = 0;
  fNofModules   = 0;
  fNofRow       = 0;
  fNofCol       = 0;
  fmodule_height= 0;
  fmodule_width = 0;
  ftower_height = 0;
  fPMTT         = 0;
  fisModule     = false;

}

dimensionCalc::~dimensionCalc() {}

G4ThreeVector dimensionCalc::GetOrigin(G4int i) {

  int row = i/fNofCol;
  int col = i%fNofCol;

  return G4ThreeVector( -fmodule_height * (double)fNofRow/2. + row * fmodule_height + fmodule_height/2., -fmodule_width * (double)fNofCol/2. + col * fmodule_width + fmodule_width/2. + fFrontL + fNofCol * fmodule_width/2., 0.);
}

G4ThreeVector dimensionCalc::GetOrigin_PMTG(G4int i) {

  int row = (i/2)/fNofCol;
  int col = (i/2)%fNofCol;

  G4ThreeVector returnVector;

  if (i%2==0) {
    returnVector = G4ThreeVector( -fmodule_height * (double)fNofRow/2. + row * fmodule_height + fmodule_height/2., -fmodule_width * (double)fNofCol/2. + col * fmodule_width + fmodule_width/2. + fFrontL + fNofCol * fmodule_width/2., ftower_height/2. + fPMTT/2.);
  } else {
    returnVector = G4ThreeVector( -fmodule_height * (double)fNofRow/2. + row * fmodule_height + fmodule_height/2., -fmodule_width * (double)fNofCol/2. + col * fmodule_width + fmodule_width/2. + fFrontL + fNofCol * fmodule_width/2., -ftower_height/2. - fPMTT/2.);
  }

  return returnVector;
}
