#include "koBICSiPMHit.hh"

G4ThreadLocal G4Allocator<koBICSiPMHit>* koBICSiPMHitAllocator = 0;

koBICSiPMHit::koBICSiPMHit(G4int wavBin, G4int timeBin)
: G4VHit(),
  fSiPMnum(0),
  fPhotons(0),
  fSiPMpos(G4ThreeVector(0.,0.,0.)),
  fMuduleNum(999),
  fTowerXY(std::make_pair(-1,-1)),
  fInnerR(0.),
  fTowerH(0.),
  fSiPMXY(std::make_pair(-1,-1)),
  fWavBin(wavBin),
  fTimeBin(timeBin)
{}

koBICSiPMHit::~koBICSiPMHit() {}

koBICSiPMHit::koBICSiPMHit(const koBICSiPMHit &right)
: G4VHit() {
  fSiPMnum = right.fSiPMnum;
  fPhotons = right.fPhotons;
  fSiPMpos = right.fSiPMpos;
  fMuduleNum = right.fMuduleNum;
  fisLeft = right.fisLeft;
  fTowerXY = right.fTowerXY;
  fInnerR = right.fInnerR;
  fTowerH = right.fTowerH;
  fSiPMXY = right.fSiPMXY;
  fWavlenSpectrum = right.fWavlenSpectrum;
  fTimeStruct = right.fTimeStruct;
}

const koBICSiPMHit& koBICSiPMHit::operator=(const koBICSiPMHit &right) {
  fSiPMnum = right.fSiPMnum;
  fPhotons = right.fPhotons;
  fSiPMpos = right.fSiPMpos;
  fMuduleNum = right.fMuduleNum;
  fisLeft = right.fisLeft;
  fTowerXY = right.fTowerXY;
  fInnerR = right.fInnerR;
  fTowerH = right.fTowerH;
  fSiPMXY = right.fSiPMXY;
  fWavlenSpectrum = right.fWavlenSpectrum;
  fTimeStruct = right.fTimeStruct;
  return *this;
}

G4bool koBICSiPMHit::operator==(const koBICSiPMHit &right) const {
  return (fSiPMnum==right.fSiPMnum && fMuduleNum==right.fMuduleNum && fSiPMXY==right.fSiPMXY);
}

void koBICSiPMHit::Draw() {}

void koBICSiPMHit::Print() {}

void koBICSiPMHit::CountWavlenSpectrum(koBICInterface::hitRange range) {
  auto it = fWavlenSpectrum.find(range);
  if (it==fWavlenSpectrum.end()) fWavlenSpectrum.insert(std::make_pair(range,1));
  else it->second++;
}

void koBICSiPMHit::CountTimeStruct(koBICInterface::hitRange range) {
  auto it = fTimeStruct.find(range);
  if (it==fTimeStruct.end()) fTimeStruct.insert(std::make_pair(range,1));
  else it->second++;
}
