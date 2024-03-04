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
  std::string z0z0hqq = "z0z0hqq"; 
  std::string z0z0qq = "z0z0qq"; 
  std::string z0zThqq = "z0zThqq"; 
  std::string z0zTqq = "z0zTqq"; 
  std::string zTzThqq = "zTzThqq"; 
  std::string zTzTqq = "zTzTqq"; 
  std::string zzhqq = "zzhqq"; 
  std::string zzqq = "zzqq"; 
  std::string z0z0hqq_C3_1 = "z0z0hqq_C3_1"; 
  std::string z0z0qq_C3_1 = "z0z0qq_C3_1"; 
  std::string z0zThqq_C3_1 = "z0zThqq_C3_1"; 
  std::string z0zTqq_C3_1 = "z0zTqq_C3_1"; 
  std::string zTzThqq_C3_1 = "zTzThqq_C3_1"; 
  std::string zTzTqq_C3_1 = "zTzTqq_C3_1"; 
  std::string zzhqq_C3_1 = "zzhqq_C3_1"; 
  std::string zzqq_C3_1 = "zzqq_C3_1"; 
  if (process_name == ttbar012j) return 88.29; 
  if (process_name == zll_123j) return 830.4; 
  if (process_name == zh_zll_hbb_012j) return 0.04718; 
  if (process_name == ttHbb) return 0.01805; 
  if (process_name == wwjj_j) return 1.254; 
  if (process_name == wzjj_j) return 0.2672; 
  if (process_name == wz_wjj_123j) return 1.615; 
  if (process_name == zzjj_j) return 0.0124; 
  if (process_name == zz_zjj_123j) return 0.4964; 
  else if (process_name == z0z0hqq) return 0.0000000182; 
  else if (process_name == z0z0qq) return 0.0000129348; 
  else if (process_name == z0zThqq) return 0.0000000805; 
  else if (process_name == z0zTqq) return 0.0000780300; 
  else if (process_name == zTzThqq) return 0.0000000899; 
  else if (process_name == zTzTqq) return 0.0001307160; 
  else if (process_name == zzhqq) return 0.0000001895; 
  else if (process_name == zzqq) return 0.0002221200; 
  else if (process_name == z0z0hqq_C3_1) return 0.0000000782; 
  else if (process_name == z0z0qq_C3_1) return 0.0000137934; 
  else if (process_name == z0zThqq_C3_1) return 0.0000001002; 
  else if (process_name == z0zTqq_C3_1) return 0.0000779940; 
  else if (process_name == zTzThqq_C3_1) return 0.0000001089; 
  else if (process_name == zTzTqq_C3_1) return 0.0001285380; 
  else if (process_name == zzhqq_C3_1) return 0.0000002862; 
  else if (process_name == zzqq_C3_1) return 0.0002199600; 
  else return 1.0; 
} 
#endif 
