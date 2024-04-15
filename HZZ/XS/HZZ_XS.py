import os
import array
import ROOT
import matplotlib.pyplot as plt

####################//####################//####################//####################

def add_backgrounds(all_backgrounds, fileDir, backgrounds):

    for background in backgrounds:
        background_file = fileDir + "/" + background
        all_backgrounds.Add(background_file)
    
    return all_backgrounds

####################//####################

def save_hists(histograms, output):

    output_file = ROOT.TFile(output, "RECREATE")

    for hist in histograms:
        hist.Write()

    print("saved to output file")
    output_file.Close()

####################//####################
    
def do_weight(event, file, treeInfo, recoBool):

    if(recoBool): #recoBool True if reco

        weight_value = getattr(event, treeInfo.recoWeight)

    else: #recoBool False if gen
        
        weight_value = getattr(event, treeInfo.genWeight)

    #lumi_value = 38 * 1000 # 2016
    #lumi_value = 45 * 1000 # 2017
    lumi_value = 64 * 1000 # 2018
    xsec_value = getattr(event, "xsec")
    totalEvents_value = file.Get("Counters").GetBinContent(40)

    #L1prefiringWeight_value = getattr(event, "L1prefiringWeight")
    #SFcorr_value = getattr(event, "SFcorr")

    totalWeight = weight_value * xsec_value * lumi_value / totalEvents_value

    return totalWeight

####################//####################

def draw(hist, drawOption, xaxis_title, yaxis_title, output_path, name):

    hist.Draw(drawOption)

    hist.GetXaxis().SetTitle(xaxis_title)
    hist.GetYaxis().SetTitle(yaxis_title)

    int = round(hist.Integral(), 2)

    latex = ROOT.TLatex()
    latex.SetTextSize(0.045)
    latex.SetTextFont(42)
    latex.SetTextColor(ROOT.kRed) 
    latex.SetTextAlign(22)
    latex.DrawLatex(5, 1030, f"Total Events: {int}")

    canvas = ROOT.gROOT.FindObject("c1")
    ROOT.gStyle.SetOptStat(00)
    ROOT.gStyle.SetStatX(0.875)
    ROOT.gStyle.SetStatY(0.875)
    canvas.Update()

    if not os.path.exists(output_path):
        os.makedirs(output_path)
    canvas.SaveAs(os.path.join(output_path, f"{name}.png"))

####################//####################

def write1DResults(txtfile, var1, var2, soverb_binContents):


    txtfile.write(f"s/b_{var1.varName}_{var2.varName}: \n")
    #txtfile.write(f"{soverb_binContents} \n")
    txtfile.write(f"VBF s/b: {soverb_binContents[3]} \n")
    txtfile.write(f"ggF s/b: {soverb_binContents[2]} \n")
    txtfile.write(f"ttH s/b: {soverb_binContents[1]} \n")
    txtfile.write(f"VH s/b: {soverb_binContents[0]} \n")

    soverb_binContents = []

    return soverb_binContents

####################//####################

def write2DResults(txtfile, var1, var2, sig_average, allbackground_average):

    txtfile.write(f"2Dratio_{var1.varName}_{var2.varName}: \n")
    txtfile.write(f"sig VBF ratio: {sig_average} \n")
    txtfile.write(f"all backgrounds ratio: {allbackground_average} \n")

    #return genReco_results

####################//####################

class var:

    def __init__(self, varName, varMin, varMax, var_binEdges):

        self.varName = varName
        self.varMin = varMin
        self.varMax = varMax
        self.var_binEdges = var_binEdges

class treeInfo:

    def __init__(self, Skimmed, UnSkimmed, treePass, treeFail, genWeight, recoWeight, XS, Counters):

        self.Skimmed = Skimmed
        self.UnSkimmed = UnSkimmed

        self.treePass = treePass
        self.treeFail = treeFail

        self.genWeight = genWeight
        self.recoWeight = recoWeight

class genSel:

    def __init__(self, sample, GENabsdetajj_min, GENabsdetajj_max, GENmjj_min, GENmjj_max, GENpT4l_min, GENpT4l_max, GENnjets_pt30_eta2p5_min, GENnjets_pt30_eta2p5_max, GENeta4l_min, GENeta4l_max):

        self.sample = sample

        self.GENabsdetajj_min = GENabsdetajj_min
        self.GENabsdetajj_max = GENabsdetajj_max

        self.GENmjj_min = GENmjj_min
        self.GENmjj_max = GENmjj_max

        self.GENpT4l_min = GENpT4l_min
        self.GENpT4l_max = GENpT4l_max

        self.GENnjets_pt30_eta2p5_min = GENnjets_pt30_eta2p5_min
        self.GENnjets_pt30_eta2p5_max = GENnjets_pt30_eta2p5_max

        self.GENeta4l_min = GENeta4l_min
        self.GENeta4l_max = GENeta4l_max
