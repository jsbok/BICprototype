#ifndef koBICInterface_h
#define koBICInterface_h 1

#include <vector>
#include <utility>
#include <map>
#include <tuple>

class koBICInterface {
public:
  koBICInterface();
  ~koBICInterface();

  typedef std::pair<float,float> hitRange;
  typedef std::pair<int,int> hitXY;
  typedef std::map<hitRange, int> koBICTimeStruct;
  typedef std::map<hitRange, int> koBICWavlenSpectrum;
  typedef std::tuple<float,float,float> threeVector;

  struct koBICModuleProperty {
    koBICModuleProperty() {};
    virtual ~koBICModuleProperty() {};

    int ModuleNum;
    koBICInterface::hitXY towerXY;
  };

  struct koBICSiPMData {
    koBICSiPMData() {};
    virtual ~koBICSiPMData() {};

    int count;
    int SiPMnum;
    int isleft;
    int x;    // plate num
    int y;    // fiber num on the plate
    threeVector pos;
    koBICTimeStruct timeStruct;
    koBICWavlenSpectrum wavlenSpectrum;
    double photonAngles;
  };

  struct koBICTowerData {
    koBICTowerData() {};
    virtual ~koBICTowerData() {};

    int ModuleNum;
    int numx;
    int numy;
    std::vector<koBICSiPMData> SiPMs;
  };

  struct koBICEdepData {
    koBICEdepData() {};
    virtual ~koBICEdepData() {};

    float EdepCore;
    float Edep;
    float EdepEle;
    float EdepGamma;
    float EdepCharged;
    int ModuleNum;
  };

  struct koBICLeakageData {
    koBICLeakageData() {};
    virtual ~koBICLeakageData() {};

    int ModuleNum;
    float EdepCore;
    float E;
    float px;
    float py;
    float pz;
    float vx;
    float vy;
    float vz;
    float vt;
    int pdgId;
    float kE;
  };

  struct koBICGenData {
    koBICGenData() {};
    virtual ~koBICGenData() {};

    float E;
    float px;
    float py;
    float pz;
    float vx;
    float vy;
    float vz;
    float vt;
    int pdgId;
  };

  struct koBICEventData {
    koBICEventData() {};
    virtual ~koBICEventData() {};

    int event_number;
    std::vector<koBICTowerData> towers;
    std::vector<koBICEdepData> Edeps;
    std::vector<koBICLeakageData> leaks;
    std::vector<koBICGenData> GenPtcs;
  };

};

#endif
