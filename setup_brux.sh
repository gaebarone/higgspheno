#source /cvmfs/sft.cern.ch/lcg/views/LCG_94/x86_64-slc6-gcc7-opt/setup.sh
#source /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-centos7-gcc8-opt/setup.sh 
#source /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-slc6-gcc8-opt/setup.sh 
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

#delphes setyp 
here=${PWD}
#change with your delphes dir
cd /isilon/data/users/gbarone1/higgsandmore/code/mg/MG5_aMC_v3_5_1/Delphes
source /isilon/data/users/gbarone1/higgsandmore/code/mg/MG5_aMC_v3_5_1/Delphes/DelphesEnv.sh
cd $here

