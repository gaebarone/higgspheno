import os
import array
import ROOT
import matplotlib.pyplot as plt

from HZZ_XS import do_weight

####################//####################

def calcGenRecoratio(hist, i, j):

    bin_content = hist.GetBinContent(i, j)
    
    nbins_x = hist.GetNbinsX()
    nbins_y = hist.GetNbinsY()
    
    sum_row = 0
    sum_col = 0
    
    for x in range(1, nbins_x + 1):
        if x != i:
            sum_row += hist.GetBinContent(x, j)
    
    for y in range(1, nbins_y + 1):
        if y != j:
            sum_col += hist.GetBinContent(i, y)
    
    if sum_row + sum_col == 0:
        ratio = 99
    else: 
        ratio = bin_content / (sum_row + sum_col)
    
    return ratio

####################//####################

def soverb1D(outdir, sig, varReco, sig_2D_hist, background_2D_hist):

    sig = sig.split('/')[0]

    sig_1D_hist = sig_2D_hist.ProjectionY("sig_1D_hist", 1, sig_2D_hist.GetNbinsY() + 1)
    background_1D_hist = background_2D_hist.ProjectionY("background_1D_hist", 1, background_2D_hist.GetNbinsY() + 1)

    sig_1D_hist.Divide(background_1D_hist)
    
    sig_1D_hist.Draw("e")
    sig_1D_hist.GetXaxis().SetTitle(varReco)
    sig_1D_hist.GetYaxis().SetTitle("s/b")
    sig_1D_hist.GetYaxis().SetRangeUser(0, 5)

    canvas = ROOT.gROOT.FindObject("c1")
    ROOT.gStyle.SetOptStat(10)
    ROOT.gStyle.SetStatX(0.875)
    ROOT.gStyle.SetStatY(0.875)

    output_path = os.path.join(outdir, "soverb1D", f"{varReco}")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    canvas.SaveAs(os.path.join(output_path,f"soverb_{sig}.png"))

####################//#################### 
    
def soverb2D(outdir, sig_bin_contents, background_bin_contents, varReco, binNumbers):

    soverbs = [x / y for x, y in zip(sig_bin_contents, background_bin_contents)]

    plt.plot(binNumbers, soverbs, marker='o', linestyle='-')
    plt.xlabel('bin number')
    plt.ylabel('s/b')
    plt.title('s/b')
    plt.ylim(0, 4)
    plt.gca().invert_xaxis()

    output_path = os.path.join(outdir, "soverb2D", f"{varReco}")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    plt.savefig(os.path.join(output_path, f'soverbBins_{varReco}.png'))

####################//#################### 

def makeGenRecoPlots(treeInfo, sig_or_background, hist2D, fileDir, outdir, treePassName, treeFailName, hists, varGen, varReco):

    file = ROOT.TFile.Open(fileDir + "/" + sig_or_background)

    tree_pass = file.Get(treePassName)
    tree_fail = file.Get(treeFailName)

    sig_or_background = sig_or_background.split('/')[0]

    passBool = True

    for event in tree_pass:

        ZZMass_value = getattr(event, "ZZMass")

        passedFiducialSelection_bbf_value = getattr(event, "passedFiducialSelection_bbf")

        if ZZMass_value > 110 or ZZMass_value < 140:
            
            varGen_value = abs(getattr(event, varGen.varName))
            varReco_value = abs(getattr(event, varReco.varName))

            totalWeight = do_weight(event, file, treeInfo, passBool)

            if passedFiducialSelection_bbf_value == 1:

                hist2D.Fill(varGen_value, varReco_value, totalWeight)

            if passedFiducialSelection_bbf_value == 0:

                hist2D.Fill(-1, varReco_value, totalWeight)

    passBool = False

    for event in tree_fail:

        GENmass4l_value = getattr(event, "GENmass4l")

        passedFiducialSelection_bbf_value = getattr(event, "passedFiducialSelection_bbf")
       
        if GENmass4l_value > 110 or GENmass4l_value < 140:
            
            varGen_value = abs(getattr(event, varGen.varName))

            totalWeight = do_weight(event, file, treeInfo, passBool)
            
            if passedFiducialSelection_bbf_value == 1:
               
                hist2D.Fill(varGen_value, -1, totalWeight)

            #if passedFiducialSelection_bbf_value == 0:

                #hist2D.Fill(-1, -1, totalWeight)

    hist2D.Draw("colz")
    hist2D.GetXaxis().SetTitle(varGen.varName)
    hist2D.GetYaxis().SetTitle(varReco.varName)
    hist2D.GetXaxis().SetRangeUser(-1, varGen.varMax)
    hist2D.GetYaxis().SetRangeUser(-1, varReco.varMax)

    lineX = ROOT.TLine(0, hist2D.GetYaxis().GetXmin(), 0, hist2D.GetYaxis().GetXmax())
    lineX.SetLineColor(ROOT.kRed)
    lineX.SetLineWidth(2)
    lineX.Draw()

    lineY = ROOT.TLine(hist2D.GetXaxis().GetXmin(), 0, hist2D.GetXaxis().GetXmax(), 0)
    lineY.SetLineColor(ROOT.kRed)
    lineY.SetLineWidth(2)
    lineY.Draw()

    hists.append(hist2D)

    #'''
    canvas = ROOT.gROOT.FindObject("c1")
    ROOT.gStyle.SetOptStat(10)
    ROOT.gStyle.SetStatX(0.875)
    ROOT.gStyle.SetStatY(0.875)
    canvas.Update()

    
    output_path = os.path.join(outdir, "GenRecoPlots", f"{varReco.varName}",f"{sig_or_background}")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    canvas.SaveAs(os.path.join(output_path,f"{sig_or_background}.png"))
    #'''

    file.Close()
    return hist2D

####################//#################### 

def doGenReco(treeInfo, sig, backgrounds, outdir, fileDir, treePassName, treeFailName, varGen, varReco, results):

    sig_hists = []
    background_hists = []
    sig_bin_contents = []
    background_bin_contents = []

    #### SIG ####

    sig_2D_hist = ROOT.TH2F(f"genReco_{sig}_{varGen.varName}_{varReco.varName}", f"genReco_{sig}_{varGen.varName}_{varReco.varName}",  len(varGen.var_binEdges)-1, varGen.varMin, varGen.varMax, len(varReco.var_binEdges)-1, varReco.varMin, varReco.varMax)
    makeGenRecoPlots(treeInfo, sig, sig_2D_hist, fileDir, outdir, treePassName, treeFailName, sig_hists, varGen, varReco)

    sig_sum = 0
    for i in range (1, sig_2D_hist.GetNbinsX()+1):
        sig_sum = sig_sum + calcGenRecoratio(sig_2D_hist, i, i)
    sig_average = sig_sum / (sig_2D_hist.GetNbinsX())

    sig_bin_content = sig_2D_hist.GetBinContent(sig_2D_hist.FindBin(varGen.varMax/2, varReco.varMax/2))
    sig_bin_contents.append(sig_bin_content)

    #### BACKGROUNDS ####

    allbackground_2D_hist = ROOT.TH2F(f"genReco_allBackground_{varGen.varName}_{varReco.varName}", f"genReco_allBackground_{varGen.varName}_{varReco.varName}", len(varGen.var_binEdges)-1, varGen.varMin, varGen.varMax, len(varReco.var_binEdges)-1, varReco.varMin, varReco.varMax)

    for background in backgrounds:

        background_2D_hist = ROOT.TH2F(f"genReco_{background}_{varGen.varName}_{varReco.varName}", f"genReco_{background}_{varGen.varName}_{varReco.varName}", len(varGen.var_binEdges)-1, varGen.varMin, varGen.varMax, len(varReco.var_binEdges)-1, varReco.varMin, varReco.varMax)
        makeGenRecoPlots(treeInfo, background, background_2D_hist, fileDir, outdir, treePassName, treeFailName, background_hists, varGen, varReco)

        background_bin_content = background_2D_hist.GetBinContent(background_2D_hist.FindBin(varGen.varMax/2, varReco.varMax/2))
        background_bin_contents.append(background_bin_content)

        allbackground_2D_hist.Add(background_2D_hist)

    allbackground_sum = 0
    for i in range (1, allbackground_2D_hist.GetNbinsX()+1):
        allbackground_sum = allbackground_sum + calcGenRecoratio(allbackground_2D_hist, i, i)
    allbackground_average = allbackground_sum / (allbackground_2D_hist.GetNbinsX())
    
    allbackground_2D_hist.Draw("colz")
    allbackground_2D_hist.GetXaxis().SetTitle(varGen.varName)
    allbackground_2D_hist.GetYaxis().SetTitle(varReco.varName)
    allbackground_2D_hist.GetXaxis().SetRangeUser(-1, varGen.varMax)
    allbackground_2D_hist.GetYaxis().SetRangeUser(-1, varReco.varMax)
    lineX = ROOT.TLine(0, allbackground_2D_hist.GetYaxis().GetXmin(), 0, allbackground_2D_hist.GetYaxis().GetXmax())
    lineX.SetLineColor(ROOT.kRed)
    lineX.SetLineWidth(2)
    lineX.Draw()
    lineY = ROOT.TLine(allbackground_2D_hist.GetXaxis().GetXmin(), 0, allbackground_2D_hist.GetXaxis().GetXmax(), 0)
    lineY.SetLineColor(ROOT.kRed)
    lineY.SetLineWidth(2)
    lineY.Draw()
    canvas = ROOT.gROOT.FindObject("c1")
    ROOT.gStyle.SetOptStat(10)
    ROOT.gStyle.SetStatX(0.875)
    ROOT.gStyle.SetStatY(0.875)
    output_path = os.path.join(outdir, "GenRecoPlots", f"{varReco.varName}", "allbackground_2D")
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    canvas.SaveAs(os.path.join(output_path, f"allbackground_2D_hist.png"))

    #### 1D s/b ####

    #soverb1D(outdir, sig, varReco.varName, sig_2D_hist, allbackground_2D_hist)

    #sig_2D_hist.Delete()
    #background_2D_hist.Delete()

    #### 2D s/b ####
        
    #soverb2D(outdir, sig_bin_contents, background_bin_contents, varReco.varName, binNumbers)

    #save_hists(sig_hists, "sig_output.root")
    #save_hists(background_hists, "background_output.root")

    results.write(f"2Dratio_{varGen.varName}_{varReco.varName}: \n")
    results.write(f"sig VBF ratio: {sig_average} \n")
    results.write(f"all backgrounds ratio: {allbackground_average} \n")

####################//####################//####################//####################
