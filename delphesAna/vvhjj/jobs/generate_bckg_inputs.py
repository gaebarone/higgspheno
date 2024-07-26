import os

process_list = ['DY2j3j', 'ttbar012j', 'ttHbb', 'wwjj_j', 'wzjj_j', 'zzjj_j', 'wz_wjj_123j', 'zz_zjj_123j']
directory_list = ['/isilon/data/common/smondal5/delphesouts_VVHjj/'] * len(process_list)

def generate_file_list(directory_path, output_file):
        with open(output_file, 'w') as file:
                # get all the files in the directory given
                filenames = os.listdir(directory_path)
                for index, filename in enumerate(filenames):
                        # loop through the file name and add them 
                        if os.path.isfile(os.path.join(directory_path, filename)):
                                file.write(directory_path+'/'+filename)
                                if index < len(filenames) - 1:
                                    file.write('\n')

def generate_file_lists(process_list, directory_list):
        for idx, process_name in enumerate(process_list):
                # get the input file lists for each process
                generate_file_list(directory_list[idx]+process_name, 'inputs/'+process_name+'_inputs.txt')
        with open('inputs/all_inputs.txt', 'w') as output:
                # write each file from above to the all_inputs file
                for index, process_name in enumerate(process_list):
                        with open('inputs/'+process_name+'_inputs.txt', 'r') as file:
                                output.write(file.read())
                                if index < len(process_list) - 1:
                                        output.write('\n')

def modify_input_list():
        # use the list of input files to generate a list of arguments to be given to the analyzer
        input_file_path = "inputs/all_inputs.txt"
        output_file_path = "inputs/bckg_inputs.txt"
        with open(input_file_path, "r") as input_file:
            input_lines = input_file.read().splitlines()
        output_lines = []
        for line in input_lines:
            parts = line.split("/")
            if len(parts) >= 3:
                modified_line = "{} /isilon/data/common/sellis9/vvhjj_outputs/{}/{} {}".format(line, parts[-2], parts[-1], parts[-2])
                output_lines.append(modified_line)
        with open(output_file_path, "w") as output_file:
            for modified_line in output_lines:
                output_file.write(modified_line + "\n")

        print("Modified content saved to", output_file_path)

def extract_first_line(input_file, output_file):
    with open(input_file, 'r') as input_text:
        first_line = input_text.readline().strip()
    with open(output_file, 'w') as output_text:
        output_text.write(first_line)

generate_file_lists(process_list, directory_list)
modify_input_list()
#extract_first_line('all_inputs.txt', 'test.txt')
