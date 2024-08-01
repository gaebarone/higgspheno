import ROOT
import os

def soverb(signal_file_name):
    # Define the directory for the input files
    base_dir = "/isilon/data/users/sellis9/higgsandmore/delphesAna/vvhjj/outputs/"

    # Define the background file name and the full paths
    background_file_name = os.path.join(base_dir, "all_bkg.root")
    signal_file_path = os.path.join(base_dir, f"{signal_file_name}.root")

    # Open the ROOT files
    signal_file = ROOT.TFile.Open(signal_file_path)
    background_file = ROOT.TFile.Open(background_file_name)

    # Check if files are successfully opened
    if not signal_file or not background_file:
        print("Error: Could not open one of the files.")
        return

    # Retrieve histograms from the files
    hist_names = ["reco", "particle", "parton"]
    hist_dict = {}
    for hist_name in hist_names:
        hist_dict[f"hSel_{hist_name}_signal"] = signal_file.Get(f"hSel_{hist_name}")
        hist_dict[f"hSel_{hist_name}_bkg"] = background_file.Get(f"hSel_{hist_name}")
        hist_dict[f"hEff_{hist_name}_signal"] = signal_file.Get(f"hEff_{hist_name}")
        hist_dict[f"hEff_{hist_name}_bkg"] = background_file.Get(f"hEff_{hist_name}")

    # Check if histograms are successfully retrieved
    for key, hist in hist_dict.items():
        if hist is None:
            print(f"Error: Could not retrieve histogram {key}.")
            return

    # Create the output directory if it doesn't exist
    output_dir = os.path.join("soverb/", signal_file_name)
    os.makedirs(output_dir, exist_ok=True)

    # Create a canvas to draw histograms
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)

    line = ROOT.TLine()
    line.SetLineColor(ROOT.kRed)
    line.SetLineStyle(2)  # Dashed line

    # Process and save each histogram
    for hist_name in hist_names:
        # Ratio plots for hSel
        hSel_ratio = hist_dict[f"hSel_{hist_name}_signal"].Clone(f"hSel_{hist_name}_ratio")
        hSel_ratio.Divide(hist_dict[f"hSel_{hist_name}_bkg"])
        max_value = hSel_ratio.GetMaximum()
        hSel_ratio.GetYaxis().SetRangeUser(0, 1.2 * max_value)
        hSel_ratio.SetStats(0)
        hSel_ratio.SetTitle(f"{signal_file_name} {hist_name} s/b (Sel)")
        hSel_ratio.Draw()
        sel_filename = os.path.join(output_dir, f"hSel_{hist_name}_ratio.png")
        canvas.SaveAs(sel_filename)
        print(f"Saved: {sel_filename}")

        # Ratio plots for hEff
        hEff_ratio = hist_dict[f"hEff_{hist_name}_signal"].Clone(f"hEff_{hist_name}_ratio")
        hEff_ratio.Divide(hist_dict[f"hEff_{hist_name}_bkg"])

        # Set the range for hEff histograms
        #hEff_ratio.GetYaxis().SetRangeUser(0, 2)
        max_value = hEff_ratio.GetMaximum()
        hEff_ratio.GetYaxis().SetRangeUser(0, 1.2 * max_value)

        for bin in range(1, hEff_ratio.GetNbinsX() + 1):
            hEff_ratio.SetBinError(bin, 0)

        hEff_ratio.SetStats(0)
        hEff_ratio.SetTitle(f"{signal_file_name} {hist_name} s/b (Eff)")

        hEff_ratio.Draw("HIST")
        line.DrawLine(hEff_ratio.GetXaxis().GetXmin(), 1, hEff_ratio.GetXaxis().GetXmax(), 1)
        eff_filename = os.path.join(output_dir, f"hEff_{hist_name}_ratio.png")
        canvas.SaveAs(eff_filename)
        print(f"Saved: {eff_filename}")

    # Clean up
    del canvas
    signal_file.Close()
    background_file.Close()
    del signal_file
    del background_file

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python soverb.py <signal_file_name>")
        sys.exit(1)
    soverb(sys.argv[1])
