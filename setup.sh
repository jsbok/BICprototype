#!/usr/bin/env bash
INSTALLDIR=$(pwd)
if [ -z ${Geant4_DIR} ]; then
  echo "$0: error: 'Geant4_DIR' should be set."
  return
fi
if [ -z ${ROOTSYS} ]; then
  echo "$0: error: ROOT environment should be set. run 'source dir_to_root/bin/thisroot.sh'"
  return
fi

export G4_VERSION=$(${Geant4_DIR}/bin/geant4-config --version)
echo "Detected Geant4 version: ${G4_VERSION}"
if [[ "${G4_VERSION}" > "11.0.0" || "${G4_VERSION}" == "11.0.0" ]]; then
  echo "This Geant4 version is not compatible with analysis provided here and simulation may not work as expected."
fi

cd ${Geant4_DIR}/share/Geant4*/geant4make
source geant4make.sh
cd $INSTALLDIR
