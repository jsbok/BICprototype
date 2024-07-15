cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(BIC)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_PREFIX_PATH "$ENV{Geant4_DIR}")

set(CMAKE_MACOSX_RPATH 1)
list(APPEND CMAKE_INSTALL_RPATH ${CMAKE_PREFIX_PATH}/lib)
list(APPEND CMAKE_INSTALL_RPATH ${PROJECT_SOURCE_DIR}/install/lib)
list(APPEND CMAKE_INSTALL_RPATH ${ROOT_DIR})

# set(HEPMC_DIR "$ENV{HEPMC_DIR}")
# set(PYTHIA_DIR "$ENV{PYTHIA_DIR}")
# set(FASTJET_DIR "$ENV{FASTJET_DIR}")

add_subdirectory(rootIO)
add_subdirectory(DRsim)
# add_subdirectory(Gen)
# add_subdirectory(Reco)
add_subdirectory(analysis)
