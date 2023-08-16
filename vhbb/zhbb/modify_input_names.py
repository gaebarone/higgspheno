# Input text file
input_file_path = "input_files.txt"

# Output text file name
output_file_path = "modified_files.txt"

# Read content from the input file and process each line
with open(input_file_path, "r") as input_file:
    input_lines = input_file.read().splitlines()

output_lines = []

for line in input_lines:
    parts = line.split("/")
    if len(parts) >= 3:
        modified_line = "{} histograms/{}/{} {}".format(line, parts[-2], parts[-1], parts[-2])
        output_lines.append(modified_line)

# Write the modified lines to the output text file
with open(output_file_path, "w") as output_file:
    for modified_line in output_lines:
        output_file.write(modified_line + "\n")

print("Modified content saved to", output_file_path)