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
  std::string DY2j3j = "DY2j3j"; 
  std::string hwpwmjj = "hwpwmjj"; 
  std::string hzzjj = "hzzjj"; 
  std::string wpwmhjj = "wpwmhjj"; 
  std::string wpwmjj = "wpwmjj"; 
  std::string zzhjj = "zzhjj"; 
  std::string zzjj = "zzjj"; 
  std::string wpwmhjj_C3_1 = "wpwmhjj_C3_1"; 
  std::string zzhjj_C3_1 = "zzhjj_C3_1"; 
  if (process_name == ttbar012j) return 88.29; 
  if (process_name == zll_123j) return 830.4; 
  if (process_name == zh_zll_hbb_012j) return 0.04718; 
  if (process_name == ttHbb) return 0.01805; 
  if (process_name == wwjj_j) return 1.254; 
  if (process_name == wzjj_j) return 0.2672; 
  if (process_name == wz_wjj_123j) return 1.615; 
  if (process_name == zzjj_j) return 0.0124; 
  if (process_name == zz_zjj_123j) return 0.4964; 
  if (process_name == DY2j3j) return 151.2; 
  else if (process_name == hwpwmjj) return 0.0005606000; 
  else if (process_name == hzzjj) return 0.0000021520; 
  else if (process_name == wpwmhjj) return 0.0000157006; 
  else if (process_name == wpwmjj) return 0.0360800000; 
  else if (process_name == zzhjj) return 0.0000002892; 
  else if (process_name == zzjj) return 0.0005808000; 
  else if (process_name == wpwmhjj_C3_1) return 0.0000198012; 
  else if (process_name == zzhjj_C3_1) return 0.0000004167; 
  else return 1.0; 
} 
#endif 
