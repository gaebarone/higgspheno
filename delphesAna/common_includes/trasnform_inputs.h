#ifndef TRASFORM_INPUTS_H
#define TRASFORM_INPUTS_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>
#include <sstream>
#include "TString.h"

// Simple JSON-like string representation
const std::string jsonString = R"(
{
  "output_names": ["label_ll", "label_bb", "label_cc"],
  "input_names": ["pf_points", "pf_features", "pf_vectors", "pf_mask"],
  "pf_points": {
    "var_names": ["part_deta1", "part_dphi1", "part_deta2", "part_dphi2"],
    "var_infos": {
      "part_deta1": {"median": 0, "norm_factor": 1},
      "part_dphi1": {"median": 0, "norm_factor": 1},
      "part_deta2": {"median": 0, "norm_factor": 1},
      "part_dphi2": {"median": 0, "norm_factor": 1}
    },
    "var_length": 128
  },
  "pf_features": {
    "var_names": ["part_pt_log", "part_e_log", "part_logptrel", "part_logerel", "part_deltaR1", "part_deltaR2", "part_charge", "part_d0", "part_d0err", "part_dz", "part_dzerr", "part_deta1", "part_dphi1", "part_deta2", "part_dphi2", "part_source"],
    "var_infos": {
      "part_pt_log": {"median": 1.7, "norm_factor": 0.7},
      "part_e_log": {"median": 2.0, "norm_factor": 0.7},
      "part_logptrel": {"median": -4.7, "norm_factor": 0.7},
      "part_logerel": {"median": -4.7, "norm_factor": 0.7},
      "part_deltaR1": {"median": 0.2, "norm_factor": 4.0},
      "part_deltaR2": {"median": 0.2, "norm_factor": 4.0},
      "part_charge": {"median": 0, "norm_factor": 1},
      "part_d0": {"median": 0, "norm_factor": 1},
      "part_d0err": {"median": 0, "norm_factor": 1},
      "part_dz": {"median": 0, "norm_factor": 1},
      "part_dzerr": {"median": 0, "norm_factor": 1},
      "part_deta1": {"median": 0, "norm_factor": 1},
      "part_dphi1": {"median": 0, "norm_factor": 1},
      "part_deta2": {"median": 0, "norm_factor": 1},
      "part_dphi2": {"median": 0, "norm_factor": 1},
      "part_source": {"median": 0, "norm_factor": 1}
    },
    "var_length": 128
  },
  "pf_vectors": {
    "var_names": ["part_px", "part_py", "part_pz", "part_energy"],
    "var_infos": {
      "part_px": {"median": 0, "norm_factor": 1},
      "part_py": {"median": 0, "norm_factor": 1},
      "part_pz": {"median": 0, "norm_factor": 1},
      "part_energy": {"median": 0, "norm_factor": 1}
    },
    "var_length": 128
  },
  "pf_mask": {
    "var_names": ["part_mask"],
    "var_infos": {
      "part_mask": {"median": 0, "norm_factor": 1}
    },
    "var_length": 128
  }
}
}
)";

// Basic JSON parser
std::unordered_map<std::string, std::vector<std::string>> parseVectorizedMap(const std::string& json) {
    std::unordered_map<std::string, std::vector<std::string>> vectorizedMap;
    std::istringstream iss(json);

    std::string line;
    while (std::getline(iss, line)) {
        size_t pos = line.find(':');
        if (pos != std::string::npos) {
            std::string key = line.substr(1, pos - 2); // Extract key between double quotes
            std::string values = line.substr(pos + 1);
            
            // Tokenize values and extract variable names
            std::istringstream valuesStream(values);
            std::string value;
            while (std::getline(valuesStream, value, ',')) {
	      vectorizedMap[ (TString(key.c_str()).ReplaceAll("\"","")).Data() ].emplace_back(value.substr(1, value.size() - 2)); // Remove double quotes
            }
        }
    }

    return vectorizedMap;
}

std::unordered_map<std::string, std::map<std::string, double>> parseTensorMap(const std::string& json) {
    std::unordered_map<std::string, std::map<std::string, double>> tensorMap;
    std::istringstream iss(json);

    std::string line;
    while (std::getline(iss, line)) {
        size_t pos = line.find(':');
        if (pos != std::string::npos) {
            std::string key = line.substr(1, pos - 2); // Extract key between double quotes
            std::string values = line.substr(pos + 1);
            
            // Tokenize values and extract variable names and medians
            std::istringstream valuesStream(values);
            std::string value;
            while (std::getline(valuesStream, value, ',')) {
                size_t colonPos = value.find(':');
                std::string varName = value.substr(1, colonPos - 2); // Extract variable name between double quotes
                double median = std::stod(value.substr(colonPos + 1));
                tensorMap[key][varName] = median;
            }
        }
    }

    return tensorMap;
}

/*
int main() {
    // Convert JSON-like string to vectorized map
    std::unordered_map<std::string, std::vector<std::string>> vectorizedMap = parseVectorizedMap(jsonString);

    // Convert JSON-like string to tensor-based map
    std::unordered_map<std::string, std::map<std::string, double>> tensorMap = parseTensorMap(jsonString);

    // Print the results for verification
    std::cout << "Vectorized Map:\n";
    for (const auto& entry : vectorizedMap) {
        std::cout << entry.first << ": ";
        for (const auto& value : entry.second) {
            std::cout << value << " ";
        }
        std::cout << "\n";
    }

    std::cout << "\nTensor Map:\n";
    for (const auto& entry : tensorMap) {
        std::cout << entry.first << ":\n";
        for (const auto& subEntry : entry.second) {
            std::cout << "  " << subEntry.first << ": " << subEntry.second << "\n";
        }
    }

    return 0;
}
*/
#endif
