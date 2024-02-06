import os
import re
from bs4 import BeautifulSoup

loopsm_directory = '/isilon/data/users/sellis9/mg_conda/MG5_aMC_v3_5_2/EVENTS/vvhjj/loopsm/qq'
hhvbf_directory = '/isilon/data/users/sellis9/mg_conda/MG5_aMC_v3_5_2/EVENTS/vvhjj/hhvbf/qq'
output_file = '/isilon/data/users/sellis9/higgsandmore/delphesAna/common_includes/get_cross_section.h'

#BRANCHING RATIOS
z_ee_BR = 2 * 0.03
z_mumu_BR = 2 * 0.03
h_bb_BR = 0.50

def extract_cross_section(html_content):
    soup = BeautifulSoup(html_content, 'html.parser')
    links = soup.find_all('a', href=True)
    for link in links:
        if link['href'] == "./Events/run_01/run_01_tag_1_banner.txt":
            next_line = link.find_parent().find_next_sibling()
            if next_line:
                return next_line.text.strip()
    return None

#def extract_cross_section_decayed(html_content):
#    soup = BeautifulSoup(html_content, 'html.parser')
#    links = soup.find_all('a', href=True)
#    for link in links:
#        if link['href'] == "./Events/run_01_decayed_1/run_01_decayed_1_tag_1_banner.txt":
#            next_line = link.find_parent().find_next_sibling()
#            if next_line:
#                return next_line.text.strip()
#    return None

def extract_value(input_str):
    if input_str is None:
        return None
    pattern = r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'
    matches = re.findall(pattern, input_str)
    return float(matches[0]) if matches else None

def scientific_to_decimal(scientific_notation):
    if scientific_notation is None:
        return None
    decimal_value = float(scientific_notation)
    decimal_notation = '{:.10f}'.format(decimal_value)
    return float(decimal_notation)

with open(output_file, 'w') as output:
    output.write(f"#ifndef GET_CROSS_SECTION_H \n")
    output.write(f"#define GET_CROSS_SECTION_H \n")
    output.write(f"#include <string> \n")

    output.write(f"\n")

    output.write(f"double get_cross_section(const char *process_name) {{\n")
    output.write(f"  std::string ttbar012j = \"ttbar012j\"; \n")
    output.write(f"  std::string zll_123j = \"zll_123j\"; \n")
    output.write(f"  std::string zh_zll_hbb_012j = \"zh_zll_hbb_012j\"; \n")
    output.write(f"  std::string ttHbb = \"ttHbb\"; \n")
    output.write(f"  std::string wwjj_j = \"wwjj_j\"; \n")
    output.write(f"  std::string wzjj_j = \"wzjj_j\"; \n")
    output.write(f"  std::string wz_wjj_123j = \"wz_wjj_123j\"; \n")
    output.write(f"  std::string zzjj_j = \"zzjj_j\"; \n")
    output.write(f"  std::string zz_zjj_123j = \"zz_zjj_123j\"; \n")

    for subdir in os.listdir(loopsm_directory):
        output.write(f"  std::string {subdir} = \"{subdir}\"; \n")

    for subdir in os.listdir(hhvbf_directory):
        output.write(f"  std::string {subdir} = \"{subdir}\"; \n")

    output.write(f"  if (process_name == ttbar012j) return 88.29; \n")
    output.write(f"  if (process_name == zll_123j) return 830.4; \n")
    output.write(f"  if (process_name == zh_zll_hbb_012j) return 0.04718; \n")
    output.write(f"  if (process_name == ttHbb) return 0.01805; \n")
    output.write(f"  if (process_name == wwjj_j) return 1.254; \n")
    output.write(f"  if (process_name == wzjj_j) return 0.2672; \n")
    output.write(f"  if (process_name == wz_wjj_123j) return 1.615; \n")
    output.write(f"  if (process_name == zzjj_j) return 0.0124; \n")
    output.write(f"  if (process_name == zz_zjj_123j) return 0.4964; \n")

    for subdir in os.listdir(loopsm_directory):
        dir_path = os.path.join(loopsm_directory, subdir)
        if os.path.isdir(dir_path):
            html_file = os.path.join(dir_path, 'crossx.html')
            if os.path.exists(html_file):
                with open(html_file, 'r') as file:
                    html_content = file.read()
                    cross_section = extract_cross_section(html_content)
                    cross_section = extract_value(cross_section)
                    if cross_section is not None:  # Check if cross_section is not None
                        cross_section_decimal = scientific_to_decimal(cross_section)
                        final_cross_section = cross_section_decimal * z_ee_BR * z_mumu_BR * h_bb_BR
                        output.write(f"  else if (process_name == {subdir}) return {final_cross_section:.10f}; \n")
                    else:
                        output.write(f"  else if (process_name == {subdir}) return 1.00; \n")
            else:
                output.write(f"{subdir}, 'crossx.html' file not found\n")

    for subdir in os.listdir(hhvbf_directory):
        dir_path = os.path.join(hhvbf_directory, subdir)
        if os.path.isdir(dir_path):
            html_file = os.path.join(dir_path, 'crossx.html')
            if os.path.exists(html_file):
                with open(html_file, 'r') as file:
                    html_content = file.read()
                    cross_section = extract_cross_section(html_content)
                    cross_section = extract_value(cross_section)
                    if cross_section is not None:  # Check if cross_section is not None
                        cross_section_decimal = scientific_to_decimal(cross_section)
                        final_cross_section = cross_section_decimal * z_ee_BR * z_mumu_BR * h_bb_BR
                        output.write(f"  else if (process_name == {subdir}) return {final_cross_section:.10f}; \n")
                    else:
                        output.write(f"  else if (process_name == {subdir}) return 1.00; \n")
            else:
                output.write(f"{subdir}, 'crossx.html' file not found\n")

    output.write(f"  else return 1.0; \n")
    output.write(f"}} \n")
    output.write(f"#endif \n")


print(f"Output saved to {output_file}")




