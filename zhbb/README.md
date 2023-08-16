# Quick steps
 - make
 - python generate_inputs.py
 - condor_submit submit_zhbb.sub
 - bash do_hadding.sh
 - root -l -b
 - .x draw_zhbb.C("../outputs/histograms/signal.root", "../outputs/histograms/all_bkg.root", "../outputs/histograms/ttbar.root", "../outputs/histograms/ttHbb.root", "../outputs/histograms/diboson.root", "../outputs/histograms/drellyan.root", "../outputs/plots")


# Extended steps for running zhbb analyzer
 - open zhbb_analyze.cpp and institute any cuts at the top of the main function
 - open do_zhbb.sh and change the paths at the top to suit your file system
 - run 'make' to compile zhbb_analyze.cpp into an executable
 - open generate_inputs and edit the process and directory lists to suit your needs (make sure all the processes have different names)
 - add each process and its corresponding cross section to get_cross_section.h
 - make a corresponding directory within the histograms directory for each sample directory being used (ensure they are all different)
 - run 'python generate_inputs.py' to automatically generate file lists for each process as well as argument lists for the condor submission
 - run 'condor_submit test_submit.sub' and check that the submission worked by viewing the job output in the condor directory
 - once the test submission runs correctly, proceed by running 'condor_submit submit_zhbb.sub'
 - open 'do_hadding.sh' and edit the hadding to fit your file system, sample directories, and background categories
 - hadd the histograms by running 'bash do_hadding.sh'
 - open 'draw_zhbb.C' and edit the background files accordingly and checking that the plots will be output where you want them
 - run draw_zhbb.C by running 'root -l -b' followed by '.x draw_zhbb.C(*insert arguments*)'
 - check output to change the order in which background categories are added to the stack if needed (smallest to largest)