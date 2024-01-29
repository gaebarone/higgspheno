###### ON BRUX ######

#source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

#here=${PWD}
#cd /isilon/data/users/sellis9/mg_cvmfs/MG5_aMC_v3_5_2/Delphes
#source /isilon/data/users/sellis9/mg_cvmfs/MG5_aMC_v3_5_2/Delphes/DelphesEnv.sh
#cd $here

###### ON MAC ######

export PATH=/opt/homebrew/Cellar/madgraph5_amcatnlo/3.4.2_2/HEPTools/bin/:$PATH
echo $PATH

here=${PWD}
cd /opt/homebrew/Cellar/madgraph5_amcatnlo/3.4.2_2/Delphes
source DelphesEnv.sh
cd $here

