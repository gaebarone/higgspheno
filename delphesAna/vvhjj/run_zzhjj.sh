make
bash do_inputs.sh
rm condor/*.err
rm condor/*.log
rm condor/*.out
rm condor/*
condor_submit submit_sig_zzhjj.sub
condor_submit submit_bckg_zzhjj.sub
