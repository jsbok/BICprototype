# BIC TB Simulation
Repository for GEANT4 simulation &amp; analysis of BIC TB simulation.

## How-to
### Prerequisite
    
#### Geant4
Geant4 10.5.0 or 10.6.0 work on MacOS 14.4.1\
Simulation on Geant4 with higher version produce invalid simulation output.\
Recommend to use recent version only for geometry visualization.

##### For local installation
Geant4 must be built with cmake option `-DGEANT4_INSTALL_DATA=ON` for simulation.\
For more information, Refer to https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/

#### CERN ROOT
any recent version may work.

### Compile

    export Geant4_DIR=/dir/to/geant4-install  # Geant4 install directory
    source /dir/to/root/bin/thisroot.sh
    
    source setup.sh
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=../install ..
    make [-jN]
    make install

### Simulation
    
    cd install
    ./bin/DRsim [geant4-macro]
This will launch interactive console and visualizatino window if enabled.\
You can manually control geant4 with the console or use predefined geant4 macro like 'run_ele.mac'.

#### Visualization
You can use below Geant4 commands in the console to control visualization.\
`/vis/viewer/zoom N` N times magnification. (N can be less than 1)\
`/vis/viewer/set/targetPoint X Y Z unit` set center coordinate of rotation and zoom. (unit: m, mm, ...)\
`/vis/viewer/set/viewpointThetaPhi T P [unit]` rotate around target point.\
`/control/execute run_ele.mac` execute external Geant4 macro.\

### Analysis
analysis_full.cc: Analysis with full option.\
analysis.cc: Analysis with simplified options.\
When running analysis.cc, import it from analysis_full.cc or use it directly to perform the corresponding task.\

    ./bin/analysis <beam_energy of 5th file> <low_edge_of_hist> <high_edge>

e.g.)

    ./bin/analysis 1000 0 20
