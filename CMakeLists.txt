cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(dual-readout)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_PREFIX_PATH "/Users/skku_server/geant4dir/geant4-install-10.6.0")

# set(HEPMC_DIR "$ENV{HEPMC_DIR}")
# set(PYTHIA_DIR "$ENV{PYTHIA_DIR}")
# set(FASTJET_DIR "$ENV{FASTJET_DIR}")

add_subdirectory(rootIO)
add_subdirectory(DRsim)
# add_subdirectory(Gen)
# add_subdirectory(Reco)
add_subdirectory(analysis)
