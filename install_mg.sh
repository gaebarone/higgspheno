source /cvmfs/sft.cern.ch/lcg/views/LCG_94/x86_64-slc6-gcc7-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/releases/ROOT/6.14.04-0d8dc/x86_64-slc6-gcc7-opt/bin/thisroot.sh 
#delphes setyp 
here=${PWD}
#change with your delphes dir
cd ../
mkdir -p mg 
cd mg
wget http://launchpad.net/madgraph5/3.0/3.5.x/+download/MG5_aMC_v3.5.1.tar.gz .
tar -xvf MG5_aMC_v3.5.1.tar.gz
cd MG5_aMC_v3.5.1 
lhapdfDIR=`lhapdf-config --prefix`
echo "set lhapdf $lhapdfDIR/bin/lhapdf-config" > startup.script
echo "install Delphes" >> startup.scrpit 
echo "install pythia8" >> startup.script 
./bin/mg5_aMC startup.script 
cd $here



