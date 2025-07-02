# source init_lcg.sh does not work by script. you should do it manually

#run the bellow command
# nohup ./1run_simulation.sh > 1_nohup.out & tail -f 1_nohup.out


#for i in "e-" "pi-" "pi0" "proton"
for i in "e-"
# particle names(e-, pi-)
do
    run_beamOn="20"

    tower_material="3by8"


    for j in 500
    # energy with the unit of GeV.
    do
        #------------------------------------------------------------#
        # events number = beamon*queue. so change "/run/beamOn 10" when you simulate different events number.
        #------------------------------------------------------------#

        macroname="${tower_material}_${i}_${j}_MeV_2025_fast"
                
        gun_particle=$i
        gun_energy="$j GeV"
        root_name="${macroname}"

        ########################################################
        results="/u/user/changhui/BICprototype2025/rootfiles/$root_name/"
        ################################^^^^^^^^################

        echo "starting submit for $gun_energy $gun_particle and the output will be $root_name.root"   
        
        # ctrl+shift+L >> select same words    

        # Length from beam generator to tower front : 1000 mm
        # x = 1000*tan(theta)

        # ---------------------------
        # for Center of Module (theta 0)
        # theta : 0.
        # phi   : -90.0
        # x(cm) : 0.
        # z(cm) : 0.
        # ---------------------------
        # ---------------------------
        # for lefttop center of Module (theta 0)
        # theta : 0.
        # phi   : -90.0
        # x(cm) : 3.
        # z(cm) : 8.
        # ---------------------------
        # ---------------------------
        # for righttop center of Module (theta 0)
        # theta : 0.
        # phi   : -90.0
        # x(cm) : 3.
        # z(cm) : -8.
        # ---------------------------
        # ---------------------------
        # for leftbottom center of Module (theta 0)
        # theta : 0.
        # phi   : -90.0
        #x(cm) : -3.
        # z(cm) : 8.
        # ---------------------------
        # ---------------------------
        # for rightbottom center of Module (theta 0)
        # theta : 0.
        # phi   : -90.0
        # x(cm) : -3.
        # z(cm) : -8.
        # ---------------------------

        echo "/DRsim/action/useHepMC False" >> 1run_$macroname.mac
        echo "/DRsim/action/useCalib True" >> 1run_$macroname.mac
        echo "/vis/disable" >> 1run_$macroname.mac
        echo "/run/numberOfThreads 1" >> 1run_$macroname.mac
        echo "/run/initialize" >> 1run_$macroname.mac
        echo "/run/verbose 1" >> 1run_$macroname.mac
#	echo "/event/verbose 2" >> 1run_$macroname.mac
#	echo "/tracking/verbose 4" >> 1run_$macroname.mac

#        echo "/DRsim/generator/randx 2" >> 1run_$macroname.mac
#        echo "/DRsim/generator/randz 2" >> 1run_$macroname.mac
#        echo "/DRsim/generator/theta 0." >> 1run_$macroname.mac
#        echo "/DRsim/generator/phi -90." >> 1run_$macroname.mac
#        echo "/DRsim/generator/x0  0." >> 1run_$macroname.mac
#        echo "/DRsim/generator/y0  0" >> 1run_$macroname.mac
#        echo "/DRsim/generator/z0  0." >> 1run_$macroname.mac
        echo "/gun/particle $i" >> 1run_$macroname.mac
        echo "/gun/energy $j GeV" >> 1run_$macroname.mac
        echo "/run/beamOn $run_beamOn" >> 1run_$macroname.mac

        echo "#! /bin/sh" > 1run_$macroname.sh
	echo "export LD_LIBRARY_PATH=/usr/lib64:\$LD_LIBRARY_PATH" >> 1run_$macroname.sh
        echo "export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.26.2/Linux-x86_64/bin:\$PATH" >> 1run_$macroname.sh
        echo "source /cvmfs/sft.cern.ch/lcg/contrib/gcc/11/x86_64-el9-gcc11-opt/setup.sh" >> 1run_$macroname.sh
        echo "source /u/user/changhui/software/root/bin/thisroot.sh" >> 1run_$macroname.sh
        echo "source /cvmfs/geant4.cern.ch/geant4/11.2/x86_64-el9-gcc11-optdeb/CMake-setup.sh" >> 1run_$macroname.sh
 
	echo "source /u/user/changhui/GEANT4/geant4-v10.5.1/install/bin/geant4.sh" >> 1run_$macroname.sh
       
        echo "export HEPMC_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/hepmc3/3.2.7/x86_64-el9-gcc11-opt" >> 1run_$macroname.sh
        echo "export FASTJET_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/fastjet/3.4.1/x86_64-el9-gcc11-opt" >> 1run_$macroname.sh
        echo "export PYTHIA_DIR=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/MCGenerators/pythia8/310/x86_64-el9-gcc11-opt" >> 1run_$macroname.sh
        echo "export PYTHIA8=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/MCGenerators/pythia8/310/x86_64-el9-gcc11-opt" >> 1run_$macroname.sh
        echo "export PYTHIA8DATA=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/MCGenerators/pythia8/310/x86_64-el9-gcc11-opt/share/Pythia8/xmldoc" >> 1run_$macroname.sh
        echo "export ROOT_INCLUDE_PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_105/hepmc3/3.2.7/x86_64-el9-gcc11-opt/include:\$ROOT_INCLUDE_PATH" >> 1run_$macroname.sh
	echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$HEPMC_DIR/lib64:\$FASTJET_DIR/lib:\$PYTHIA_DIR/lib:\$PWD/lib" >> 1run_$macroname.sh
 #     	echo "echo \"[DEBUG] LD_LIBRARY_PATH = \$LD_LIBRARY_PATH\"" >> 1run_$macroname.sh
 #        echo "      ldd ./bin/DRsim" >> 1run_$macroname.sh
        echo "./bin/DRsim 1run_$macroname.mac \$1 $results/root/$root_name" >> 1run_$macroname.sh
     
        echo "getenv = True" >> 1run_$macroname.sub
        echo "universe = vanilla" >> 1run_$macroname.sub
        echo "executable = 1run_$macroname.sh" >> 1run_$macroname.sub
        echo "arguments = \$(ProcId)" >> 1run_$macroname.sub
        echo "core_size = 0" >> 1run_$macroname.sub
        echo "output = $results/log/\$(ProcId).out" >> 1run_$macroname.sub
        echo "error = $results/log/\$(ProcId).err" >> 1run_$macroname.sub
        echo "log = $results/log/\$(ProcId).log" >> 1run_$macroname.sub
        echo "request_memory = 1 GB" >> 1run_$macroname.sub
        echo "should_transfer_files = YES" >> 1run_$macroname.sub
        echo "when_to_transfer_output = ON_EXIT" >> 1run_$macroname.sub
        echo "transfer_input_files =./bin, ./lib, ./init.mac, ./1run_$macroname.mac" >> 1run_$macroname.sub
        echo "queue 500" >> 1run_$macroname.sub

        #cat 0_SimulationShFile.sh
        mkdir "$results"
        mkdir "$results/log"
        mkdir "$results/root"
        condor_submit 1run_$macroname.sub

        echo "Simulation Uploaded"
        echo "$i $j GeV $tower_material"

        cp 1run_$macroname.* $results
    done
    
done

