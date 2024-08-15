import uproot
import matplotlib.pyplot as plt
import numpy as np
import os

# sig                                                                                                                                                                                                                                                                                                                                                                                                                                                          
sig = 'wpwmhqq'
sig_scale = 100

sig_file = sig + '.root'
sig_name = sig + ' x ' + str(sig_scale)

# input dir ( .root )                                                                                                                                                                                                                                                                                                                                                                                                                                          
input_dir = '/isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/'

input_files = [
    {'path': input_dir + sig_file, 'name': sig_name, 'color': 'mediumseagreen'},
    {'path': input_dir + 'diboson.root', 'name': 'diboson', 'color': 'orange'},
    {'path': input_dir + 'drellyan.root', 'name': 'drellyan', 'color': 'mediumblue'},
    {'path': input_dir + 'ttHbb.root', 'name': 'ttHbb', 'color': 'mediumpurple'},
    {'path': input_dir + 'ttbar.root', 'name': 'ttbar', 'color': 'tomato'}
]

# output ( .png )                                                                                                                                                                                                                                                                                                                                                                                                                                              
output_dir = os.path.join('/isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/histograms/stacks/', sig)

def read_root_files(input_files):
    histograms = {}

    for file_info in input_files:
        with uproot.open(file_info['path']) as file:
            for key in file.keys():
                if "TH1F" in str(file[key]):
                    hist = file[key].to_numpy()
                    if key in histograms:
                        histograms[key].append((hist, file_info['name'], file_info['color']))
                    else:
                        histograms[key] = [(hist, file_info['name'], file_info['color'])]
    return histograms

def plot_histograms(histograms, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for key, hist_list in histograms.items():
        plt.figure(figsize=(12, 8))
        bin_edges = hist_list[0][0][1]
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        stacked_heights = np.zeros_like(bin_centers)

        for i, (hist, name, color) in enumerate(hist_list):
            hist_values = hist[0].copy()
            hist_values += 1e-10  # Add small constant to avoid zeros
            if name == sig_name:
                hist_values *= sig_scale  # scale signal histogram                                                                                                                                                                                                                                                                                                                                                                                             
            plt.bar(bin_centers, hist_values, width=np.diff(bin_edges), align='center',
                    bottom=stacked_heights, label=name, color=color)
            stacked_heights += hist_values

        plt.minorticks_on()
        plt.grid(True)
        plt.xlabel(' ')
        plt.ylabel('Events')
        plt.yscale('log')
        modified_key = key.rstrip('1')
        plt.title(f'{modified_key.replace(";", "").replace("/", "_")}')
        plt.legend()

        plot_filename = os.path.join(output_dir, f"{key.replace(';', '').replace('/', '_')}.png")
        plt.savefig(plot_filename)
        plt.close()

def main(input_files, output_dir):
    histograms = read_root_files(input_files)
    if not histograms:
        print("No histograms found.")
        return
    plot_histograms(histograms, output_dir)

if __name__ == "__main__":
    main(input_files, output_dir)
