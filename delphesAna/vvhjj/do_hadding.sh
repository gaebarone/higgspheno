cd outputs/histograms
rm *.root
rm wwjj_j/wwjj_j.root
rm wzjj_j/wzjj_j.root
rm wz_wjj_123j/wz_wjj_123j.root
rm zzjj_j/zzjj_j.root
rm zz_zjj_123j/zz_zjj_123j.root
hadd ttbar.root ttbar012j/del*.root
hadd ttHbb.root ttHbb/del*.root
hadd wwjj_j/wwjj_j.root wwjj_j/del*.root
hadd wzjj_j/wzjj_j.root wzjj_j/del*.root
hadd wz_wjj_123j/wz_wjj_123j.root wz_wjj_123j/del*.root
hadd zzjj_j/zzjj_j.root zzjj_j/del*.root
hadd zz_zjj_123j/zz_zjj_123j.root zz_zjj_123j/del*.root
hadd diboson.root wwjj_j/wwjj_j.root wzjj_j/wzjj_j.root wz_wjj_123j/wz_wjj_123j.root zzjj_j/zzjj_j.root zz_zjj_123j/zz_zjj_123j.root
hadd drellyan.root zll_123j/del*.root
hadd signal.root zh_zll_hbb_012j/del*.root
hadd all_bkg.root ttbar.root ttHbb.root diboson.root drellyan.root
cd ../..
