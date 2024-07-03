# BIC TB Simulation
Repository for GEANT4 simulation &amp; analysis of BIC TB simulation.

## How-to
### Prerequisite
    
Geant4 10.5.0 or 10.6.0 confirmed on MacOS 14.4.1\
Simulation on Geant4 with higher version produce invalid output.\
Recommend to use it only for geometry visualization.

#### For local installation
Geant4 must be built with cmake option `-DGEANT4_INSTALL_DATA=ON` for simulation.
\
For more information, Refer to https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/
\
ROOT 6.30 confirmed on MacOS 14.4.1
Other versions may work.
\
Other prerequisites not listed here can exists.

### Compile

    export Geant4_DIR=/dir/to/geant4-install
    export ROOT_DIR=/dir/to/ROOT
    
    source setup.sh
    mkdir build && cd build
    cmake -DCMAKE_INSTALL_PREFIX=../install ..
    make [-jN]
    make install

### Simulation
    
    cd install
    ./bin/DRsim 
This will launch interactive console and visualizatino window if enabled.

#### Visualization
You can use below Geant4 commands to control visualization.
\
`/vis/viewer/zoom N` N times magnification. (N can be less than 1)\
`/vis/viewer/set/targetPoint X Y Z unit` set center coordinate of rotation and zoom. (unit: m, mm, ...)\
`/vis/viewer/set/viewpointThetaPhi T P [unit]` rotate around target point.\
`/control/execute run_ele.mac` execute external Geant4 macro.\

### Analysis

    ./bin/analysis <path_to_root_file> <low_edge_of_hist> <high_edge>

e.g.)

    ./bin/analysis /home/USER/20GeV_ele_data 0 20
