#!/bin/bash

cd outputs 

./do_hadding.sh

cd ../histograms 

./make_plots.sh
./copy_plots.sh

cd ..
