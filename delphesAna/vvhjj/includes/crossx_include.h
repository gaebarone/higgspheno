#ifndef GET_CROSS_SECTION_H 
#define GET_CROSS_SECTION_H 
#include <string> 

double get_cross_section(const char *process_name) {
<<<<<<< HEAD
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
  std::string wwhqq = "wwhqq";
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
  else if (process_name == wwhqq) return 0.0000133; //2.66e-05/2
  else return 1.0;
}
#endif

=======
  std::string ttbar012j = "ttbar012j"; 
  std::string zll_123j = "zll_123j"; 
  std::string zh_zll_hbb_012j = "zh_zll_hbb_012j"; 
  std::string ttHbb = "ttHbb"; 
  std::string wwjj_j = "wwjj_j"; 
  std::string wzjj_j = "wzjj_j"; 
  std::string wz_wjj_123j = "wz_wjj_123j"; 
  std::string zzjj_j = "zzjj_j"; 
  std::string zz_zjj_123j = "zz_zjj_123j"; 
  std::string w+0w-0hqq = "w+0w-0hqq"; 
  std::string w+0w-0qq = "w+0w-0qq"; 
  std::string w+0w-Thqq = "w+0w-Thqq"; 
  std::string w+0w-Tqq = "w+0w-Tqq"; 
  std::string w+Tw-0hqq = "w+Tw-0hqq"; 
  std::string w+Tw-0qq = "w+Tw-0qq"; 
  std::string w+Tw-Thqq = "w+Tw-Thqq"; 
  std::string w+Tw-Tqq = "w+Tw-Tqq"; 
  std::string w+w-hqq = "w+w-hqq"; 
  std::string w+w-qq = "w+w-qq"; 
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
  else if (process_name == w+0w-0hqq) return 0.0000010000; 
  else if (process_name == w+0w-0qq) return 0.0007974000; 
  else if (process_name == w+0w-Thqq) return 0.0000022160; 
  else if (process_name == w+0w-Tqq) return 0.0024560000; 
  else if (process_name == w+Tw-0hqq) return 0.0000023540; 
  else if (process_name == w+Tw-0qq) return 0.0026280000; 
  else if (process_name == w+Tw-Thqq) return 0.0000053540; 
  else if (process_name == w+Tw-Tqq) return 0.0087520000; 
  else if (process_name == w+w-hqq) return 0.0000109220; 
  else if (process_name == w+w-qq) return 0.0146360000; 
  else if (process_name == z0z0hqq) return 0.0000002020; 
  else if (process_name == z0z0qq) return 0.0001437200; 
  else if (process_name == z0zThqq) return 0.0000008946; 
  else if (process_name == z0zTqq) return 0.0008670000; 
  else if (process_name == zTzThqq) return 0.0000009994; 
  else if (process_name == zTzTqq) return 0.0014524000; 
  else if (process_name == zzhqq) return 0.0000021060; 
  else if (process_name == zzqq) return 0.0024680000; 
  else return 1.0; 
} 
#endif 
>>>>>>> 739c677419ee00783f2dd7295120bf00935cf2b4
