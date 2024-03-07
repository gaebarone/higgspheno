#!/bin/bash

#outdir="/isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs"

source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

#LHAPDFDIR=/isilon/data/users/sellis9/mg_cvmfs/MG5_aMC_v3_5_2/HEPTools/lhapdf6_py3
#export PATH=${LHAPDFDIR}/bin:$PATH
#export LD_LIBRARY_PATH=${LHAPDFDIR}/lib:$LD_LIBRARY_PATH
#export PYTHONPATH=${LHAPDFDIR}/lib/python3.10/site-packages/:$PYTHONPATH
#export PYTHIA8DATA=/isilon/data/users/sellis9/mg_cvmfs/MG5_aMC_v3_5_2/HEPTools/pythia8/share/Pythia8/xmldoc

here=${PWD}
cd /isilon/data/users/sellis9/mg_cvmfs/MG5_aMC_v3_5_2/Delphes
source /isilon/data/users/sellis9/mg_cvmfs/MG5_aMC_v3_5_2/Delphes/DelphesEnv.sh
cd $here

mkdir -p execute/job_$2
cd execute/job_$2
#cp -r /isilon/data/users/sellis9/mg_cvmfs/MG5_aMC_v3_5_2/Delphes/* .

ulimit -s unlimited
cd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj

./zAnalyzer $1 $2 $3 $4

#cp -p $2 $outdir/$2

cd $here/execute/
rm -r job_$2
