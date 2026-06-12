#ifndef koBICMagneticField_h
#define koBICMagneticField_h 1

#include "globals.hh"
#include "G4MagneticField.hh"

class G4GenericMessenger;

class koBICMagneticField : public G4MagneticField {
public:
  koBICMagneticField(G4int);
  virtual ~koBICMagneticField();

  virtual void GetFieldValue(const G4double point[4],double* bField ) const;

  void SetField(G4double val) { fBy = val; }
  G4double GetField() const { return fBy; }

private:
  void DefineCommands(G4int);

  G4GenericMessenger* fMessenger;
  G4double fBy;
};

#endif
