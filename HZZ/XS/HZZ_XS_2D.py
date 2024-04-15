import os
import array
import ROOT
import matplotlib.pyplot as plt

from HZZ_XS import do_weight
from HZZ_XS import draw

####################//####################//####################//####################

def unroll2D(hist2D, unrolledHist): 

    t=1; 

    for i in range (1, hist2D.GetNbinsX()+1) : 
        for j in range (1, hist2D.GetNbinsY()+1) : 
             
             unrolledHist.SetBinContent(t,hist2D.GetBinContent(i,j))
             unrolledHist.SetBinError(t,hist2D.GetBinError(i,j))
                                      
             t=t+1 

    return unrolledHist          

####################//####################//####################//####################

def do2D(treeInfo, fileDir, outdir, sig_or_background, var1, var2, treePass, treeFail, recoBool): # reco = true, gen = false

    outdir = outdir + "/2D"

    file = ROOT.TFile.Open(fileDir + "/" + sig_or_background)
    tPass = file.Get(treePass)
    tFail = file.Get(treeFail)

    sig_or_background = sig_or_background.split('/')[0]

    histname = f"{sig_or_background}_{var1.varName}_vs_{var2.varName}"

    hist = ROOT.TH2F(histname, histname, len(var1.var_binEdges)-1,  array.array('d', var1.var_binEdges), len(var2.var_binEdges)-1, array.array('d', var2.var_binEdges))
    unrolledHist = ROOT.TH1F(f"{histname}_unrolled", f"{histname}_unrolled", hist.GetNbinsX()* hist.GetNbinsY(), 0, hist.GetNbinsX()* hist.GetNbinsY())


    ''' if recoBool == False:

        for event in tPass:
            
            if getattr(event, "passedFiducialSelection_bbf") == 1:

                if 115 < getattr(event, "GENmass4l") < 135:

                    if 115 < getattr(event, "ZZMass") < 135:

                        branch1_value = getattr(event, var1.varName)
                        branch2_value = getattr(event, var2.varName)
                        totalWeight = do_weight(event, file, treeInfo, recoBool) #RECO
                        hist.Fill(branch1_value, branch2_value, totalWeight)

        for event in tFail:

            if  getattr(event, "passedFiducialSelection_bbf") == 1:

                if 115 < getattr(event, "GENmass4l") < 135:

                    branch1_value = getattr(event, var1.varName)
                    branch2_value = getattr(event, var2.varName)
                    totalWeight = do_weight(event, file, treeInfo, recoBool) #RECO
                    hist.Fill(branch1_value, branch2_value, totalWeight)

    if recoBool == True:

        for event in tPass:

            if getattr(event, "passedFiducialSelection_bbf") == 0:

                if 115 < getattr(event, "ZZMass") < 135:

                    branch1_value = getattr(event, var1.varName)
                    branch2_value = getattr(event, var2.varName)
                    totalWeight = do_weight(event, file, treeInfo, recoBool) #RECO
                    hist.Fill(branch1_value, branch2_value, totalWeight) '''


    for event in tPass:

            if 115 < getattr(event, "ZZMass") < 135:
            
                branch1_value = getattr(event, var1.varName)
                branch2_value = getattr(event, var2.varName)
                totalWeight = do_weight(event, file, treeInfo, recoBool) #RECO
                hist.Fill(branch1_value, branch2_value, totalWeight)


#   for i in range(1, hist.GetNbinsX() + 1):
#       for j in range(1, hist.GetNbinsY() + 1):
#           bin_content = hist.GetBinContent(i, j)
#           bin_width_x = hist.GetXaxis().GetBinWidth(i)
#           bin_width_y = hist.GetYaxis().GetBinWidth(j)
        
#           normalized_content = bin_content / (bin_width_x * bin_width_y)
#           hist.SetBinContent(i, j, normalized_content)


    draw(hist, " ", var1.varName, var2.varName, outdir, histname)

    unroll2D(hist, unrolledHist)

    draw(unrolledHist, " ", var1.varName, var2.varName, f"{outdir}_unrolled", f"{histname}_unrolled")

    file.Close()

####################//####################//####################//####################