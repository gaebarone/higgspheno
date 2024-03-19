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
  std::string wp0wm0hqq = "wp0wm0hqq"; 
  std::string wp0wm0qq = "wp0wm0qq"; 
  std::string wp0wmThqq = "wp0wmThqq"; 
  std::string wp0wmTqq = "wp0wmTqq"; 
  std::string wpTwm0hqq = "wpTwm0hqq"; 
  std::string wpTwm0qq = "wpTwm0qq"; 
  std::string wpTwmThqq = "wpTwmThqq"; 
  std::string wpTwmTqq = "wpTwmTqq"; 
  std::string wpwmhqq = "wpwmhqq";
  std::string wpwmqq = "wpwmqq"; 
  std::string z0z0hqq = "z0z0hqq"; 
  std::string z0z0qq = "z0z0qq"; 
  std::string z0zThqq = "z0zThqq"; 
  std::string z0zTqq = "z0zTqq"; 
  std::string zTzThqq = "zTzThqq"; 
  std::string zTzTqq = "zTzTqq"; 
  std::string zzhqq = "zzhqq"; 
  std::string zzqq = "zzqq"; 
  if (process_name == ttbar012j) return 88.29; 
  if (process_name == zll_123j) return 830.4; 
  if (process_name == zh_zll_hbb_012j) return 0.04718; 
  if (process_name == ttHbb) return 0.01805; 
  if (process_name == wwjj_j) return 1.254; 
  if (process_name == wzjj_j) return 0.2672; 
  if (process_name == wz_wjj_123j) return 1.615; 
  if (process_name == zzjj_j) return 0.0124; 
  if (process_name == zz_zjj_123j) return 0.4964; 
  else if (process_name == wp0wm0hqq) return 0.0000011600; 
  else if (process_name == wp0wm0qq) return 0.0009249840; 
  else if (process_name == wp0wmThqq) return 0.0000025706; 
  else if (process_name == wp0wmTqq) return 0.0028489600; 
  else if (process_name == wpTwm0hqq) return 0.0000027306; 
  else if (process_name == wpTwm0qq) return 0.0030484800; 
  else if (process_name == wpTwmThqq) return 0.0000062106; 
  else if (process_name == wpTwmTqq) return 0.0101523200; 
  else if (process_name == wpwmhqq) return 0.0000126695; 
  else if (process_name == wpwmqq) return 0.0169777600; 
  else if (process_name == z0z0hqq) return 0.0000000211; 
  else if (process_name == z0z0qq) return 0.0000150044; 
  else if (process_name == z0zThqq) return 0.0000000934; 
  else if (process_name == z0zTqq) return 0.0000905148; 
  else if (process_name == zTzThqq) return 0.0000001043; 
  else if (process_name == zTzTqq) return 0.0001516306; 
  else if (process_name == zzhqq) return 0.0000002199; 
  else if (process_name == zzqq) return 0.0002576592; 
  else return 1.0; 
} 
#endif 
