#ifndef koBICSiPMHit_h
#define koBICSiPMHit_h 1

#include "koBICInterface.h"

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class koBICSiPMHit : public G4VHit {
public:

  koBICSiPMHit(G4int wavBin, G4int timeBin);
  koBICSiPMHit(const koBICSiPMHit &right);
  virtual ~koBICSiPMHit();

  const koBICSiPMHit& operator=(const koBICSiPMHit &right);
  G4bool operator==(const koBICSiPMHit &right) const;

  inline void *operator new(size_t);
  inline void operator delete(void* aHit);

  void Draw();
  void Print();

  void photonCount() { fPhotons++; }
  G4int GetPhotonCount() const { return fPhotons; }

  void SetSiPMnum(G4int n) { fSiPMnum = n; }
  G4int GetSiPMnum() const { return fSiPMnum; }

  void SetSiPMpos(G4ThreeVector pos) { fSiPMpos = pos; }
  G4ThreeVector GetSiPMpos() const { return fSiPMpos; }

  void SetModuleNum(G4int MuduleNum) { fMuduleNum = MuduleNum; }
  G4int GetModuleNum() const { return fMuduleNum; }

  void SetisLeft(G4int isLeft) { fisLeft = isLeft; }
  G4int GetisLeft() const { return fisLeft; }

  void SetTowerXY(koBICInterface::hitXY xy) { fTowerXY = xy; }
  koBICInterface::hitXY GetTowerXY() const { return fTowerXY; }

  // void SetTowerInnerR(G4float innerR) { fInnerR = innerR; }
  // G4float GetTowerInnerR() const { return fInnerR; }

  // void SetTowerH(G4float towerH) { fTowerH = towerH; }
  // G4float GetTowerH() const { return fTowerH; }

  void SetSiPMXY(koBICInterface::hitXY xy) { fSiPMXY = xy; }
  koBICInterface::hitXY GetSiPMXY() const { return fSiPMXY; }

  void CountWavlenSpectrum(koBICInterface::hitRange range);
  koBICInterface::koBICWavlenSpectrum GetWavlenSpectrum() const { return fWavlenSpectrum; }

  void CountTimeStruct(koBICInterface::hitRange range);
  koBICInterface::koBICTimeStruct GetTimeStruct() const { return fTimeStruct; }

private:
  G4int fSiPMnum;
  G4int fPhotons;
  G4ThreeVector fSiPMpos;
  G4int fMuduleNum;
  G4int fisLeft;
  koBICInterface::hitXY fTowerXY;
  G4float fInnerR;
  G4float fTowerH;
  koBICInterface::hitXY fSiPMXY;
  koBICInterface::koBICWavlenSpectrum fWavlenSpectrum;
  koBICInterface::koBICTimeStruct fTimeStruct;

  G4int fWavBin;
  G4int fTimeBin;
};

typedef G4THitsCollection<koBICSiPMHit> koBICSiPMHitsCollection;
extern G4ThreadLocal G4Allocator<koBICSiPMHit>* koBICSiPMHitAllocator;

inline void* koBICSiPMHit::operator new(size_t) {
  if (!koBICSiPMHitAllocator) koBICSiPMHitAllocator = new G4Allocator<koBICSiPMHit>;
  return (void*)koBICSiPMHitAllocator->MallocSingle();
}

inline void koBICSiPMHit::operator delete(void*aHit) {
  koBICSiPMHitAllocator->FreeSingle((koBICSiPMHit*) aHit);
}

#endif
