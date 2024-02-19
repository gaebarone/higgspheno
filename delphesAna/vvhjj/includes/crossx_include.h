#ifndef CROSSX_INCLUDE_H
#define CROSSX_INCLUDE_H

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
  std::string zzhqq = "zzhqq";
  std::string z0z0hqq = "z0z0hqq";
  std::string z0zThqq = "z0zThqq";
  std::string zTzThqq = "zTzThqq";
  std::string zzhqq_C3_1 = "zzhqq_C3_1";
  std::string z0z0hqq_C3_1 = "z0z0hqq_C3_1";
  std::string z0zThqq_C3_1 = "z0zThqq_C3_1";
  std::string zTzThqq_C3_1 = "zTzThqq_C3_1";
  if (process_name == ttbar012j) return 88.29;
  else if (process_name == zll_123j) return 830.4;
  else if (process_name == zh_zll_hbb_012j) return 0.04718;
  else if (process_name == ttHbb) return 0.01805;
  else if (process_name == wwjj_j) return 1.254;
  else if (process_name == wzjj_j) return 0.2672;
  else if (process_name == wz_wjj_123j) return 1.615;
  else if (process_name == zzjj_j) return 0.0124;
  else if (process_name == zz_zjj_123j) return 0.4964;
  else if (process_name == zzhqq) return 0.0000002118;        // 2.118e-07
  else if (process_name == z0z0hqq) return 0.00000002306;     // 2.306e-08
  else if (process_name == z0zThqq) return 0.0000001004;      // 1.004e-07
  else if (process_name == zTzThqq) return 0.00000008934;     // 8.934e-08
  else if (process_name == zzhqq_C3_1) return 0.000002619;    // 2.619e-06
  else if (process_name == z0z0hqq_C3_1) return 0.000002337;  // 2.337e-06 
  else if (process_name == z0zThqq_C3_1) return 0.0000002032; // 2.032e-07
  else if (process_name == zTzThqq_C3_1) return 0.0000000769; // 7.69e-08
  else return 1.0;
}
#endif

