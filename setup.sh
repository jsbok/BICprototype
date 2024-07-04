INSTALLDIR=$(pwd)
if [ -z ${Geant4_DIR} ]; then
  echo "$0: error: 'Geant4_DIR' should be set."
  return
fi
if [ -z ${ROOT_DIR} ]; then
  echo "$0: error: 'ROOT_DIR' should be set."
  return
fi
cd ${Geant4_DIR}/share/Geant4*/geant4make
source geant4make.sh
cd $INSTALLDIR
