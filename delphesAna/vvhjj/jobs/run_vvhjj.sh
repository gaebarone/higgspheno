cd ..
make
cd jobs

# clean
rm inputs/*.txt
rm logs/*.err
rm logs/*.log
rm logs/*.out

# generate inputs
python generate_sig_inputs.py
python generate_bckg_inputs.py

# submit 
condor_submit submit_sig_vvhjj.sub
condor_submit submit_bckg_vvhjj.sub
