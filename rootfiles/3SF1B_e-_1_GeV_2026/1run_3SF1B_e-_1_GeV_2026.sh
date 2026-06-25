#! /bin/sh
export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.26.2/Linux-x86_64/bin:$PATH
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/11/x86_64-el9-gcc11-opt/setup.sh
source /u/user/changhui/software/root/bin/thisroot.sh
source /cvmfs/geant4.cern.ch/geant4/11.2/x86_64-el9-gcc11-optdeb/CMake-setup.sh
source /u/user/changhui/GEANT4/geant4-v10.5.1/install/bin/geant4.sh
export HEPMC_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/hepmc3/3.2.7/x86_64-el9-gcc11-opt
export FASTJET_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/fastjet/3.4.1/x86_64-el9-gcc11-opt
export PYTHIA_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/MCGenerators/pythia8/310/x86_64-el9-gcc11-opt
export PYTHIA8=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/MCGenerators/pythia8/310/x86_64-el9-gcc11-opt
export PYTHIA8DATA=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/MCGenerators/pythia8/310/x86_64-el9-gcc11-opt/share/Pythia8/xmldoc
export ROOT_INCLUDE_PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/hepmc3/3.2.7/x86_64-el9-gcc11-opt/include:$ROOT_INCLUDE_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HEPMC_DIR/lib64:$FASTJET_DIR/lib:$PYTHIA_DIR/lib:$PWD/lib
./bin/koBIC 1run_3SF1B_e-_1_GeV_2026.mac $1 /u/user/changhui/koBIC2026/BICprototype/rootfiles/3SF1B_e-_1_GeV_2026//root/3SF1B_e-_1_GeV_2026
