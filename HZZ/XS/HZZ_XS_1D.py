import os
import array
import ROOT
import matplotlib.pyplot as plt

from HZZ_XS import do_weight
from HZZ_XS import draw

####################//####################//####################//####################

def sel1Gen(file, treeInfo, event, hist, detajj, mjj): # abs deta jj vs. m jj

    detajj_value = getattr(event, detajj.varName)
    mjj_value = getattr(event, mjj.varName)
    totalWeight = do_weight(event, file, treeInfo, True)

    if detajj.var_binEdges[3] < detajj_value < detajj.var_binEdges[4] or mjj.var_binEdges[3] < mjj_value < mjj.var_binEdges[4]: # VBF
    
        hist.Fill(3, totalWeight) 

    if detajj.var_binEdges[2] < detajj_value < detajj.var_binEdges[3] or mjj.var_binEdges[2] < mjj_value < mjj.var_binEdges[3]: # ggF
    
        hist.Fill(2, totalWeight)

    if detajj.var_binEdges[0] < detajj_value < detajj.var_binEdges[1] or mjj.var_binEdges[1] < mjj_value < mjj.var_binEdges[2]: # ttH
    
        hist.Fill(1, totalWeight)

    if detajj.var_binEdges[1] < detajj_value < detajj.var_binEdges[2] or mjj.var_binEdges[0] < mjj_value < mjj.var_binEdges[1]: # VH
    
        hist.Fill(0, totalWeight)

    return hist

####################//####################//####################//####################

def sel2Gen(file, treeInfo, event, hist, ptH, Nj): # pT H vs. N j

    ptH_value = getattr(event, ptH.varName)
    Nj_value = getattr(event, Nj.varName)
    totalWeight = do_weight(event, file, treeInfo, True)

    if ptH.var_binEdges[2] < ptH_value < ptH.var_binEdges[3] or Nj_value == Nj.var_binEdges[2]: # VBF
    
        hist.Fill(3, totalWeight) 

    if ptH.var_binEdges[0] < ptH_value < ptH.var_binEdges[1] or Nj_value == Nj.var_binEdges[0]: # ggH
    
        hist.Fill(2, totalWeight)

    if ptH.var_binEdges[3] < ptH_value < ptH.var_binEdges[4] or Nj_value > Nj.var_binEdges[2]: # ttH
    
        hist.Fill(1, totalWeight)

    if ptH.var_binEdges[1] < ptH_value < ptH.var_binEdges[2] or Nj_value == Nj.var_binEdges[1]: # VH
    
        hist.Fill(0, totalWeight)

    return hist

####################//####################//####################//####################

def sel3Gen(file, treeInfo, event, hist, ptH, etaH): # pT H vs. eta H

    ptH_value = getattr(event, ptH.varName)
    etaH_value = getattr(event, etaH.varName)
    totalWeight = do_weight(event, file, treeInfo, True)

    if ptH.var_binEdges[2] < ptH_value < ptH.var_binEdges[3] or etaH.var_binEdges[1] < etaH_value < etaH.var_binEdges[2]: # VBF
    
        hist.Fill(3, totalWeight) 

    if ptH.var_binEdges[0] < ptH_value < ptH.var_binEdges[1] or etaH.var_binEdges[3] < etaH_value < etaH.var_binEdges[4]: # ggH
    
        hist.Fill(2, totalWeight)

    if ptH.var_binEdges[3] < ptH_value < ptH.var_binEdges[4] or etaH.var_binEdges[0] < etaH_value < etaH.var_binEdges[1]: # ttH
    
        hist.Fill(1, totalWeight)

    if ptH.var_binEdges[1] < ptH_value < ptH.var_binEdges[2] or etaH.var_binEdges[2] < etaH_value < etaH.var_binEdges[3]: # VH
    
        hist.Fill(0, totalWeight)

    return hist

####################//####################//####################//####################


def sel1Reco(file, treeInfo, event, hist, detajj, mjj): # abs deta jj vs. m jj

    detajj_value = getattr(event, detajj)
    mjj_value = getattr(event, mjj)
    totalWeight = do_weight(event, file, treeInfo, True)

    if 3 < detajj_value < 10 or 300 < mjj_value < 1000: # VBF
    
        hist.Fill(3, totalWeight) 

    if 0 < detajj_value < 1.75 or 125 < mjj_value < 200: # ggF
    
        hist.Fill(2, totalWeight)

    if 0 < detajj_value < 1 or 100 < mjj_value < 150: # ttH
    
        hist.Fill(1, totalWeight)

    if 0 < detajj_value < 1 or 60 < mjj_value < 100: # VH
    
        hist.Fill(0, totalWeight)

    return hist

####################//####################//####################//####################

def sel2Reco(file, treeInfo, event, hist, ptH, Nj): # pT H vs. N j

    ptH_value = getattr(event, ptH)
    Nj_value = getattr(event, Nj)
    totalWeight = do_weight(event, file, treeInfo, True)

    if 50 < ptH_value < 100 or 0 < Nj_value < 1: # VBF
    
        hist.Fill(3, totalWeight) 

    if 0 < ptH_value < 50 or 0 < Nj_value < 1: # ggH
    
        hist.Fill(2, totalWeight)

    if 50 < ptH_value < 150 or 2 < Nj_value < 4: # ttH
    
        hist.Fill(1, totalWeight)

    if 25 < ptH_value < 125 or 0 < Nj_value < 2: # VH
    
        hist.Fill(0, totalWeight)

    return hist

####################//####################//####################//####################

def sel3Reco(file, treeInfo, event, hist, ptH, etaH): # pT H vs. eta H

    ptH_value = getattr(event, ptH)
    etaH_value = getattr(event, etaH)
    totalWeight = do_weight(event, file, treeInfo, True)

    if 50 < ptH_value < 100 or 1 < etaH_value < 2: # VBF
    
        hist.Fill(3, totalWeight) 

    if 0 < ptH_value < 50 or 2 < etaH_value < 4: # ggH
    
        hist.Fill(2, totalWeight)

    if 50 < ptH_value < 150 or 0 < etaH_value < 2: # ttH
    
        hist.Fill(1, totalWeight)

    if 25 < ptH_value < 125 or 0 < etaH_value < 3: # VH
    
        hist.Fill(0, totalWeight)

    return hist

####################//####################//####################//####################

def do1D(treeInfo, fileDir, outdir, sig, backgrounds, var1, var2, soverb_binContents):

    outdir = outdir + "/1D"

    sig_file = ROOT.TFile.Open(fileDir + "/" + sig)

    treePass_sig = sig_file.Get("candTree")
    treeFail_sig = sig_file.Get("candTree_failed")

    sig_name = sig.split('/')[0]
    sig_histName = f"1D_{sig_name}_{var1.varName}_{var2.varName}"
   
    sigHist = ROOT.TH1F(sig_histName, sig_histName, 4, 0, 4)
    allbackgroundHist = ROOT.TH1F(f"1D_allbackgrounds_{var1.varName}_{var2.varName}", f"1D_allbackgrounds_{var1.varName}_{var2.varName}", 4, 0, 4)

    for event in treePass_sig:

        if getattr(event, "passedFiducialSelection_bbf") == 1:

            if var1.varName == "GENabsdetajj" and var2.varName == "GENmjj":
                sel1Gen(sig_file, treeInfo, event, sigHist, var1, var2)
            if var1.varName == "GENpT4l" and var2.varName == "GENnjets_pt30_eta2p5":
                sel2Gen(sig_file, treeInfo, event, sigHist, var1, var2)
            if var1.varName == "GENpT4l" and var2.varName == "GENeta4l":
                sel3Gen(sig_file, treeInfo, event, sigHist, var1, var2)

        if getattr(event, "passedFiducialSelection_bbf") == 0:

            if var1.varName == "absdetajj" and var2.varName == "mjj":
                sel1Reco(sig_file, treeInfo, event, sigHist, var1, var2)
            if var1.varName == "ZZPt" and var2.varName == "njets_pt30_eta2p5":
                sel2Reco(sig_file, treeInfo, event, sigHist, var1, var2)
            if var1.varName == "ZZPt" and var2.varName == "ZZEta":
                sel3Reco(sig_file, treeInfo, event, sigHist, var1, var2)

    for event in treeFail_sig:

        if getattr(event, "passedFiducialSelection_bbf") == 1:

            if var1.varName == "GENabsdetajj" and var2.varName == "GENmjj":
                sel1Gen(sig_file, treeInfo, event, sigHist, var1, var2)
            if var1.varName == "GENpT4l" and var2.varName == "GENnjets_pt30_eta2p5":
                sel2Gen(sig_file, treeInfo, event, sigHist, var1, var2)
            if var1.varName == "GENpT4l" and var2.varName == "GENeta4l":
                sel3Gen(sig_file, treeInfo, event, sigHist, var1, var2)

    draw(sigHist,  " ", f"F({var1.varName},{var2.varName})", "Events",  f"{outdir}/{sig_name}", f"events_{sig_histName}")

    for background in backgrounds:

        background_file = ROOT.TFile.Open(fileDir + "/" + background)

        background_name = background.split('/')[0]
        background_histName = f"1D_{background_name}_{var1.varName}_{var2.varName}"

        treePass_background = background_file.Get("candTree")
        treeFail_background = background_file.Get("candTree_failed")

        backgroundHist = ROOT.TH1F(f"{background}_1D", f"{background}_1D", 4, 0, 4)

        for event in treePass_background:

            if getattr(event, "passedFiducialSelection_bbf") == 1:

                if var1.varName == "GENabsdetajj" and var2.varName == "GENmjj":
                    sel1Gen(background_file, treeInfo, event, backgroundHist, var1, var2)
                if var1.varName == "GENpT4l" and var2.varName == "GENnjets_pt30_eta2p5":
                    sel2Gen(background_file, treeInfo, event, backgroundHist, var1, var2)
                if var1.varName == "GENpT4l" and var2.varName == "GENeta4l":
                    sel3Gen(background_file, treeInfo, event, backgroundHist, var1, var2)

            if getattr(event, "passedFiducialSelection_bbf") == 0:

                if var1.varName == "absdetajj" and var2.varName == "mjj":
                    sel1Reco(background_file, treeInfo, event, backgroundHist, var1, var2)
                if var1.varName == "ZZPt" and var2.varName == "njets_pt30_eta2p5":
                    sel2Reco(background_file, treeInfo, event, backgroundHist, var1, var2)
                if var1.varName == "ZZPt" and var2.varName == "ZZEta":
                    sel3Reco(background_file, treeInfo, event, backgroundHist, var1, var2)

        for event in treeFail_background:

            if getattr(event, "passedFiducialSelection_bbf") == 1:

                if var1.varName == "GENabsdetajj" and var2.varName == "GENmjj":
                    sel1Gen(background_file, treeInfo, event, backgroundHist, var1, var2)
                if var1.varName == "GENpT4l" and var2.varName == "GENnjets_pt30_eta2p5":
                    sel2Gen(background_file, treeInfo, event, backgroundHist, var1, var2)
                if var1.varName == "GENpT4l" and var2.varName == "GENeta4l":
                    sel3Gen(background_file, treeInfo, event, backgroundHist, var1, var2)

        allbackgroundHist.Add(backgroundHist)

    draw(backgroundHist,  " ", f"F({var1.varName},{var2.varName})",  "Events",  f"{outdir}/{background_name}", f"events_{background_histName}")
        
    draw(allbackgroundHist,  " ", f"F({var1.varName},{var2.varName})",  "Events",  outdir, f"events_allbackgrounds_{var1.varName}_{var2.varName}")

    sigHist.Divide(allbackgroundHist)

    sigHist.GetXaxis().SetRangeUser(0, 4)
    #sigHist.GetYaxis().SetRangeUser(0, 2)

    draw(sigHist,  " ", f"F({var1.varName},{var2.varName})", "s/b",  outdir, f"soverb_{sig_histName}")

    for i in range(1, sigHist.GetNbinsX() + 1):
        soverb_binContents.append(sigHist.GetBinContent(i))

    sig_file.Close()

    for background in backgrounds:

        background_file.Close()

    return soverb_binContents

####################//####################//####################//####################