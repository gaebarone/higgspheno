import os

process_list = ['zh_zll_hbb_012j', 'ttHbb', 'ttbar012j', 'wwjj_j', 'wz_wjj_123j', 'wzjj_j', 'zll_123j', 'zz_zjj_123j', 'zzjj_j']
directory_list = ['/home/trussel1/MG2/Madgraph/MG5_aMC_v2_6_7/Delphes/data_delphes/']+['/isilon/data/common/smondal5/delphesouts_VVHjj/'] * 8

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
                generate_file_list(directory_list[idx]+process_name, process_name+'_inputs.txt')
        with open('all_inputs.txt', 'w') as output:
                # write each file from above to the all_inputs file
                for index, process_name in enumerate(process_list):
                        with open(process_name+'_inputs.txt', 'r') as file:
                                output.write(file.read())
                                if index < len(process_list) - 1:
                                        output.write('\n')

def modify_input_list():
        # use the list of input files to generate a list of arguments to be given to the analyzer
        input_file_path = "all_inputs.txt"
        output_file_path = "modified_inputs.txt"
        with open(input_file_path, "r") as input_file:
            input_lines = input_file.read().splitlines()
        output_lines = []
        for line in input_lines:
            parts = line.split("/")
            if len(parts) >= 3:
                modified_line = "{} histograms/{}/{} {}".format(line, parts[-2], parts[-1], parts[-2])
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
extract_first_line('all_inputs.txt', 'test.txt')