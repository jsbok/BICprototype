#include "koBICMagneticField.hh"

#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

koBICMagneticField::koBICMagneticField(G4int Bnum)
: G4MagneticField(), fMessenger(0), fBy(0.5*tesla)
{
  // define commands for this class
  DefineCommands(Bnum);
}

koBICMagneticField::~koBICMagneticField() {
  delete fMessenger;
}

void koBICMagneticField::GetFieldValue(const G4double [4],double *bField) const {
  bField[0] = 0.;
  bField[1] = fBy;
  bField[2] = 0.;
}

void koBICMagneticField::DefineCommands(G4int Bnum) {
  // Define /koBIC/field command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/koBIC/Magneticfield/", "Field control");

  // fieldValue command
  G4GenericMessenger::Command& valueCmd = fMessenger->DeclareMethodWithUnit("value","tesla", &koBICMagneticField::SetField, "Set field strength.");
  valueCmd.SetParameterName("field", true);
  valueCmd.SetDefaultValue("1.");
}
