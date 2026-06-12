#----------------------------------------------------------------------------
# Setup the project
project(koBIC)

set(CMAKE_MACOSX_RPATH 1)
list(APPEND CMAKE_INSTALL_RPATH ${CMAKE_PREFIX_PATH}/lib)
list(APPEND CMAKE_INSTALL_RPATH ${ROOT_DIR})
list(APPEND CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})

message("Detected Geant4 Version: ${Geant4_VERSION}")
if("${Geant4_VERSION}" VERSION_GREATER_EQUAL "11.0.0")
    add_definitions(-DG4_LATEST)
    message("Status: Geant4 11+ detected. DG4_LATEST flag enabled.")
endif()
add_definitions(-DG4VIS_USE -DG4UI_USE)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(
  # ${HEPMC_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${Geant4_INCLUDE_DIR}
)
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(koBIC koBIC.cc ${sources} ${headers})
target_link_libraries(
  koBIC ${Geant4_LIBRARIES}
  # ${HEPMC_DIR}/lib64/libHepMC3.so
  # ${HEPMC_DIR}/lib64/libHepMC3rootIO.so
  rootIO
)

file(GLOB koBIC_MACROS ${PROJECT_SOURCE_DIR}/*.mac)
file(COPY ${koBIC_MACROS} DESTINATION ${PROJECT_BINARY_DIR})

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
#add_custom_target(koBIC)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS koBIC DESTINATION bin)
install(FILES ${koBIC_MACROS} DESTINATION .)
