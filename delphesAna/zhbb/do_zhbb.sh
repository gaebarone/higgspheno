#!/bin/bash

outdir="/home/trussel1/MG2/Madgraph/MG5_aMC_v2_6_7/Delphes/examples/outputs/"

MGPATH=/home/trussel1/MG2/Madgraph
source /cvmfs/sft.cern.ch/lcg/views/LCG_94/x86_64-slc6-gcc7-opt/setup.sh
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/trussel1/MG2/Madgraph/lhapdf_build/lib/
PATH=/home/trussel1/MG2/Madgraph/lhapdf_build/bin:$PATH
LD_LIBRARY_PATH=/home/trussel1/MG2/Madgraph/lhapdf_build/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/home/trussel1/MG2/Madgraph/lhapdf_build/lib/python2.7/site-packages:$PYTHONPATH
export PYTHIA8DATA=/home/trussel1/MG2/Madgraph/MG5_aMC_v2_6_7/HEPTools/pythia8/share/Pythia8/xmldoc/

mkdir -p execute/job_$2
cd execute/job_$2
pwd
cp -r /home/trussel1/MG2/Madgraph/MG5_aMC_v2_6_7/Delphes/* .

ulimit -s unlimited
cd examples/higgsandmore/delphesAna/zhbb

./zhbb_analyze $1 $2 $3

cp -p $2 $outdir/$2

cd ../../../../..
rm -r job_$2