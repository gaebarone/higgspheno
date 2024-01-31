selection='HZZJJ'

cd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs
rm *.root

cd /isilon/data/common/sellis9/zAnalyzerOutputs/histograms
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/ttbar_${selection}.root ttbar012j/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/ttHbb_${selection}.root ttHbb/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/wwjj_j_${selection}.root wwjj_j/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/wzjj_j_${selection}.root wzjj_j/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/wz_wjj_123j_${selection}.root wz_wjj_123j/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/zzjj_j_${selection}.root zzjj_j/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/zz_zjj_123j_${selection}.root zz_zjj_123j/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/drellyan_${selection}.root zll_123j/${selection}_del*.root

hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/diboson_${selection}.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/wwjj_j_${selection}.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/wzjj_j_${selection}.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/wz_wjj_123j_${selection}.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/zzjj_j_${selection}.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/zz_zjj_123j_${selection}.root

hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/all_bkg_${selection}.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/ttbar_${selection}.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/ttHbb_${selection}.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/diboson_${selection}.root /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/drellyan_${selection}.root


hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/zzhqq_${selection}.root zzhqq/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/z0z0hqq_${selection}.root z0z0hqq/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/z0zThqq_${selection}.root z0zThqq/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/zTzThqq_${selection}.root zTzThqq/${selection}_del*.root

hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/zzhqq_C3_1_${selection}.root zzhqq_C3_1/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/z0z0hqq_C3_1_${selection}.root z0z0hqq_C3_1/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/z0zThqq_C3_1_${selection}.root z0zThqq_C3_1/${selection}_del*.root
hadd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs/zTzThqq_C3_1_${selection}.root zTzThqq_C3_1/${selection}_del*.root

cd /isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/rootOutputs
