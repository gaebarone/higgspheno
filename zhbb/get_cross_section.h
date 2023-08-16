double get_cross_section(const char *process_name) {
  if (process_name == "ttbar_012j") return 88.29;
  else if (process_name == "zll_123j") return 830.4;
  else if (process_name == "zh_zll_hbb_012j") return 0.04718;
  else if (process_name == "ttHbb") return 0.01790922;
  else if (process_name == "wwjj_j") return 1.254;
  else if (process_name == "wzjj_j") return 0.2672;
  else if (process_name == "wz_wjj_123j") return 1.615;
  else if (process_name == "zzjj_j") return 0.0124;
  else if (process_name == "zz_zjj_123j") return 0.4964;
  else return 1;
}