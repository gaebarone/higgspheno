#ifndef GET_CROSS_SECTION_H
#define GET_CROSS_SECTION_H
#include <string>

double get_cross_section(const char *process_name) {
  std::string ttbar012j = "ttbar012j";
  std::string zll_123j = "zll_123j";
  std::string zh_zll_hbb_012j = "zh_zll_hbb_012j";
  std::string ttHbb = "ttHbb";
  std::string wwjj_j = "wwjj_j";
  std::string wzjj_j = "wzjj_j";
  std::string wz_wjj_123j = "wz_wjj_123j";
  std::string zzjj_j = "zzjj_j";
  std::string zz_zjj_123j = "zz_zjj_123j";
  if (process_name == ttbar012j) return 88.29;
  else if (process_name == zll_123j) return 830.4;
  else if (process_name == zh_zll_hbb_012j) return 0.04718;
  else if (process_name == ttHbb) return 0.01805;
  else if (process_name == wwjj_j) return 1.254;
  else if (process_name == wzjj_j) return 0.2672;
  else if (process_name == wz_wjj_123j) return 1.615;
  else if (process_name == zzjj_j) return 0.0124;
  else if (process_name == zz_zjj_123j) return 0.4964;
  else return 1.0;
}
#endif
