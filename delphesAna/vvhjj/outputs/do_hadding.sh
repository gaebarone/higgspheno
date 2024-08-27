cd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs
rm *.root

cd /isilon/data/common/sellis9/vvhjj_outputs
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/ttbar.root ttbar012j/del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/ttHbb.root ttHbb/del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/wwjj_j.root wwjj_j/del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/wzjj_j.root wzjj_j/del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/wz_wjj_123j.root wz_wjj_123j/del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/zzjj_j.root zzjj_j/del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/zz_zjj_123j.root zz_zjj_123j/del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/drellyan.root DY2j3j/del*.root

hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/diboson.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/wwjj_j.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/wzjj_j.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/wz_wjj_123j.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/zzjj_j.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/zz_zjj_123j.root

hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/all_bkg.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/ttbar.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/ttHbb.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/diboson.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/drellyan.root


# vvhqq
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/zzhjj.root zzhjj/del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/wpwmhjj.root wpwmhjj/del*.root

# vvqq
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/zzjj.root zzjj/del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/wpwmjj.root wpwmjj/del*.root

# hvvqq ( h > vv )
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/hzzjj.root hzzjj/del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/hwpwmjj.root hwpwmjj/del*.root

cd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs
