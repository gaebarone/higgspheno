import numpy as np
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import os

# sig
sig_file = ''
sig_name = 'test'

# input dir ( .root )
input_dir = '/isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/'

input_files = [
    {'path': input_dir + 'diboson_HZZJJ.root', 'name': 'diboson', 'color': 'mediumseagreen'},
    {'path': input_dir + 'drellyan_HZZJJ.root', 'name': 'drellyan', 'color': 'lightskyblue'},
    {'path': input_dir + 'ttHbb_HZZJJ.root', 'name': 'ttHbb', 'color': 'mediumpurple'},
    {'path': input_dir + 'ttbar_HZZJJ.root', 'name': 'ttbar', 'color': 'hotpink'}
]

# output ( .png )
output_dir = os.path.join('/isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/histograms/', sig_name)

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
        hep.style.use("CMS")

        fig, ax = plt.subplots()

        bin_edges = hist_list[0][0][1]
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        stacked_heights = np.zeros_like(bin_centers)

        for hist, name, color in hist_list:
            ax.bar(bin_centers, hist[0], width=np.diff(bin_edges), align='center',
                   bottom=stacked_heights, label=name, color=color, alpha=0.8)
            stacked_heights += hist[0]

        ax.set_ylabel("Events")
        ax.set_title(f'{key.replace(";", "").replace("/", "_")}')
        
        plt.legend()
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1])
        
        hep.cms.label()

        plot_filename = os.path.join(output_dir, f"{key.replace(';', '').replace('/', '_')}.png")
        plt.savefig(plot_filename)
        print(f"Saved plot to: {plot_filename}")
        plt.close()

def main(input_files, output_dir):
    histograms = read_root_files(input_files)
    if not histograms:
        print("No histograms found.")
        return
    plot_histograms(histograms, output_dir)

if __name__ == "__main__":
    main(input_files, output_dir)
