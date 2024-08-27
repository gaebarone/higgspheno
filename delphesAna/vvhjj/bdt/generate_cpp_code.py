import sys
import argparse
import re

def get_single_booster_cpp_code(booster_tree, branch_id, class_index, indentation_level=0):
    level = booster_tree[branch_id].split(' ')

    booster_code = ""

    if 'leaf' in level[0]:
        # Extract the leaf value and add to the sum
        leaf_value = float(level[0].split('=')[1])
        booster_code += "{0}sum[{1}] += {2};\n".format("  " * indentation_level, class_index, leaf_value)
        return booster_code

    # Parse the branching information
    branch_ids = level[1].split(',')
    yes_branch_id = int(branch_ids[0].split("=")[1])
    no_branch_id = int(branch_ids[1].split("=")[1])
    missing_branch_id = int(branch_ids[2].split("=")[1])

    # Extract feature name and comparison
    condition = level[0]

    # Extract the feature name between brackets
    feature_name_part = condition.split(':[', 1)[1].split(']', 1)[0].strip()
    
    # Remove any extraneous characters around the feature name and threshold
    feature_name = feature_name_part.split('<')[0].split('>')[0].strip()

    # Determine the comparison operator and threshold value
    if '<' in condition:
        comparison = '<'
        threshold_value = condition.split('<', 1)[1].strip()
    else:
        comparison = '>'
        threshold_value = condition.split('>', 1)[1].strip()

    # Ensure no trailing characters in the threshold value
    threshold_value = threshold_value.rstrip(']')

    # Generate the C++ code
    booster_code += "{0}if (sample[\"{1}\"] {2} {3}) {{\n".format("  " * indentation_level, feature_name, comparison, threshold_value)
    booster_code += get_single_booster_cpp_code(booster_tree, yes_branch_id, class_index, indentation_level + 1)
    booster_code += "{0}}} else {{\n".format("  " * indentation_level)
    booster_code += get_single_booster_cpp_code(booster_tree, no_branch_id, class_index, indentation_level + 1)
    booster_code += "{0}}}\n".format("  " * indentation_level)

    return booster_code

def generate_single_booster_cpp_code(booster, class_index):
    booster_tree = {}
    for line in booster:
        branch_id = int(line.split(':')[0].strip())
        booster_tree[branch_id] = line.strip()

    return get_single_booster_cpp_code(booster_tree, 0, class_index, 1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--xgb_dump', type=str, default='dump.raw.txt', help='Raw boosters dump. Created without passing feature map file to XGBoost dump() function.')
    parser.add_argument('--num_classes', type=int, required=True, help='Number of classes this model is classifying')
    parser.add_argument('--output', type=str, required=True, help='Output file and location')

    args = parser.parse_args()

    result = ""

    result += "#include <iostream>\n"
    result += "#include <fstream>\n"
    result += "#include <vector>\n"
    result += "using namespace std;\n\n"
    result += "std::vector<double> xgb_classify( std::map<std::string, double>  &sample) {\n\n"
    result += "  std::vector<double> sum ({0}, 0.0);\n\n".format(args.num_classes)

    booster_counter = 0
    boosters = []
    with open(args.xgb_dump, 'r') as f:
        for line in f:
            if 'booster' in line:
                boosters.append([])
                booster_counter += 1
            else:
                boosters[booster_counter - 1].append(line.strip())

    for index, booster in enumerate(boosters):
        class_index = index % args.num_classes
        result += generate_single_booster_cpp_code(booster, class_index)
        result += "\n\n"

    result += "  return sum;\n"
    result += "}\n\n"

    with open(args.output, 'w') as f:
        f.write(result)
