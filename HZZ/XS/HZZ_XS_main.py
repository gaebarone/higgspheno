import os
import array
import ROOT
import matplotlib.pyplot as plt

from HZZ_XS import var
from HZZ_XS import treeInfo
from HZZ_XS import genSel

from HZZ_XS import write1DResults
from HZZ_XS import write2DResults

from HZZ_XS_1D import do1D
from HZZ_XS_2D import do2D
from HZZ_XS_GenReco import doGenReco

####################//####################//####################//####################
####################//                 VERSION 00               //####################
####################//####################//####################//#################### 
#'''
outdir = "test"

# var 1 = jj eta
var1Gen_binEdges = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]
var1Reco_binEdges = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5, 9, 9.5, 10]

# var 2 = mjj
var2Gen_binEdges = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000]
var2Reco_binEdges = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000]

# var 3 = higgs pT / (var3Gen = "GENH_pt")
var3Gen_binEdges = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500]
var3Reco_binEdges = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500]

# var 4 = n jets
var4Gen_binEdges = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
var4Reco_binEdges = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

# var 5 = higgs eta / (var5Gen = "GENH_eta")
var5Gen_binEdges = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5]
var5Reco_binEdges = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5]
#'''
####################//####################//####################//####################
####################//                 VERSION 0                //####################
####################//####################//####################//####################
'''
outdir = "ver0"

# var 1 = jj eta
# var 2 = mjj
# var 3 = higgs pT / (var3Gen = "GENH_pt")
# var 4 = n jets
# var 5 = higgs eta / (var5Gen = "GENH_eta")

var1Gen_binEdges = [0, 2.5, 5, 7.5, 10]
var2Gen_binEdges = [0, 250, 500, 750, 1000]
var3Gen_binEdges = [0, 125, 250, 375, 500]
var4Gen_binEdges = [0, 2.5, 5, 7.5, 10]
var5Gen_binEdges = [0, 1.25, 2.5, 3.75, 5]

var1Reco_binEdges = [0, 2.5, 5, 7.5, 10]
var2Reco_binEdges = [0, 250, 500, 750, 1000]
var3Reco_binEdges = [0, 125, 250, 375, 500]
var4Reco_binEdges = [0, 2.5, 5, 7.5, 10]
var5Reco_binEdges = [0, 1.25, 2.5, 3.75, 5]

#VBFGenSel = GenSel(VBF, GENabsdetajj_min, GENabsdetajj_max, GENmjj_min, GENmjj_max, GENpT4l_min, GENpT4l_max, GENnjets_pt30_eta2p5_min, GENnjets_pt30_eta2p5_max, GENeta4l_min, GENeta4l_max)

'''
####################//####################//####################//####################
####################//                 VERSION 1                //####################
####################//####################//####################//####################
'''
outdir = "ver1"

# var 1 = jj eta
# var 2 = mjj
# var 3 = higgs pT / (var3Gen = "GENH_pt")
# var 4 = n jets
# var 5 = higgs eta / (var5Gen = "GENH_eta")

var1Gen_binEdges = [0, 0.5, 1.5, 2.5, 10]
var2Gen_binEdges = [0, 100, 150, 250, 1000]
var3Gen_binEdges = [0, 50, 100, 150, 500]
var4Gen_binEdges = [0, 1, 2, 3, 10]
var5Gen_binEdges = [0, 1, 1.5, 2.25, 5]

var1Reco_binEdges = [0, 0.5, 1.5, 2.25, 10]
var2Reco_binEdges = [0, 100, 150, 250, 1000]
var3Reco_binEdges = [0, 50, 100, 150, 500]
var4Reco_binEdges = [0, 1, 2, 3, 10]
var5Reco_binEdges = [0, 1, 1.5, 2.25, 5]

#VBFGenSel = GenSel(VBF, GENabsdetajj_min, GENabsdetajj_max, GENmjj_min, GENmjj_max, GENpT4l_min, GENpT4l_max, GENnjets_pt30_eta2p5_min, GENnjets_pt30_eta2p5_max, GENeta4l_min, GENeta4l_max)

'''
####################//####################//####################//####################
####################//                 VERSION 2                //####################
####################//####################//####################//####################
'''
outdir = "ver2"

# var 1 = jj eta
# var 2 = mjj
# var 3 = higgs pT / (var3Gen = "GENH_pt")
# var 4 = n jets
# var 5 = higgs eta / (var5Gen = "GENH_eta")

var1Gen_binEdges = [0, 0.5, 1.5, 2.25, 10]
var2Gen_binEdges = [0, 100, 150, 250, 1000]
var3Gen_binEdges = [0, 40, 100, 140, 500]
var4Gen_binEdges = [0, 1, 2, 3, 10]
var5Gen_binEdges = [0, 1.1, 1.6, 2.25, 5]

var1Reco_binEdges = [0, 0.5, 1.5, 2.25, 10]
var2Reco_binEdges = [0, 100, 150, 250, 1000]
var3Reco_binEdges = [0, 40, 100, 140, 500]
var4Reco_binEdges = [0, 1, 2, 3, 10]
var5Reco_binEdges = [0, 1.1, 1.6, 2.25, 5]

#VBFGenSel = GenSel(VBF, GENabsdetajj_min, GENabsdetajj_max, GENmjj_min, GENmjj_max, GENpT4l_min, GENpT4l_max, GENnjets_pt30_eta2p5_min, GENnjets_pt30_eta2p5_max, GENeta4l_min, GENeta4l_max)

'''
####################//####################//####################//####################
####################//####################//####################//####################

# skimmed

var1Gen = var("GENabsdetajj", 0, 10, var1Gen_binEdges)
var1Reco = var("absdetajj", 0, 10, var1Reco_binEdges)

var2Gen = var("GENmjj", 0, 1000, var2Gen_binEdges)
var2Reco = var("mjj", 0, 1000, var2Reco_binEdges)

var3Gen = var("GENpT4l", 0, 500, var3Gen_binEdges)
var3Reco = var("ZZPt", 0, 500, var3Reco_binEdges)

var4Gen = var("GENnjets_pt30_eta2p5", 0, 10, var4Gen_binEdges)
var4Reco = var("njets_pt30_eta2p5", 0, 10, var4Reco_binEdges)

var5Gen = var("GENeta4l", 0, 5, var5Gen_binEdges)
var5Reco = var("ZZEta", 0, 5, var5Reco_binEdges)

fileDir="/Users/spencerellis/Desktop/research/hzz/hzz_xs/data/2018UL"

treePassName = "candTree"
treeFailName = "candTree_failed"

sig ="VBFH125/VBFH125_reducedTree_MC_2018.root"
backgrounds = ["ggH125/ggH125_reducedTree_MC_2018.root", "ttH125/ttH125_reducedTree_MC_2018.root", "WplusH125/WplusH125_reducedTree_MC_2018.root", "WminusH125/WminusH125_reducedTree_MC_2018.root", "ZH125/ZH125_reducedTree_MC_2018.root"]
zz = "ZZTo4l/ZZTo4l_reducedTree_MC_2018.root"

treeInfo = treeInfo(True, False, treePassName, treeFailName, "genHEPMCweight", "overallEventWeight", "xsec", "Counters")
#treeInfo = treeInfo(True, False, treePassName, treeFailName, "genHEPMCweight", "weight", "xsec", "Counters")

####################//####################//####################//####################
'''
# unskimmed

var1Gen = var(" ", 0, 10, var1Gen_binEdges)
var1Reco = var("DiJetDEta", 0, 10, var1Reco_binEdges)

var2Gen = var(" ", 0, 1000, var2Gen_binEdges)
var2Reco = var("DiJetMass", 0, 1000, var2Reco_binEdges)

var3Gen = var("GenHPt", 0, 500, var3Gen_binEdges)
var3Reco = var("ZZPt", 0, 500, var3Reco_binEdges)

var4Gen = var(" ", 0, 10, var4Gen_binEdges)
var4Reco = var("nCleanedJetsPt30", 0, 10, var4Reco_binEdges)

var5Gen = var("GenHRapidity", 0, 5, var5Gen_binEdges)
var5Reco = var("ZZEta", 0, 5, var5Reco_binEdges)

fileDir="/Users/spencerellis/Desktop/research/hzz/hzz_xs/data"

treePassName = "ZZTree/candTree"
treeFailName = "ZZTree/candTree_failed"

sig_unskimmed = "MC_2018_VBFH125.root"
backgrounds_unskimmed = ["MC_2018_ggH125.root", "MC_2018_ttH125.root", "MC_2018_ZH125.root"]
sig_runIII = "RunIII_VBFH125.root"

treeInfo = treeInfo(False, True, treePassName, treeFailName, "overallEventWeight", "overallEventWeight", "xsec", "ZZTree/Counters")
'''
####################//####################//####################//####################

soverb_binContents = []

results = open(f"{outdir}_results.txt", "w")

####################//####################//####################//####################

'''
do1D(treeInfo, fileDir, outdir, sig, backgrounds, var1Gen, var2Gen, soverb_binContents)
write1DResults(results, var1Gen, var2Gen, soverb_binContents)
do1D(treeInfo, fileDir, outdir, sig, backgrounds, var3Gen, var4Gen, soverb_binContents)
write1DResults(results, var3Gen, var4Gen, soverb_binContents)
do1D(treeInfo, fileDir, outdir, sig, backgrounds, var3Gen, var5Gen, soverb_binContents)
write1DResults(results, var3Gen, var5Gen, soverb_binContents)
'''

do2D(treeInfo, fileDir, outdir, sig, var1Gen, var2Gen, treePassName, treeFailName, False) # false for gen
do2D(treeInfo, fileDir, outdir, sig, var3Gen, var4Gen, treePassName, treeFailName, False) # false for gen
#do2D(treeInfo, fileDir, outdir, sig, var3Gen, var5Gen, treePassName, treeFailName, False) # false for gen
do2D(treeInfo, fileDir, outdir, sig, var1Reco, var2Reco, treePassName, treeFailName, True)
do2D(treeInfo, fileDir, outdir, sig, var3Reco, var4Reco, treePassName, treeFailName, True)
#do2D(treeInfo, fileDir, outdir, sig, var3Reco, var5Reco, treePassName, treeFailName, True)

for background in backgrounds:
    do2D(treeInfo, fileDir, outdir, background, var1Gen, var2Gen, treePassName, treeFailName, False) # false for gen
    do2D(treeInfo, fileDir, outdir, background, var3Gen, var4Gen, treePassName, treeFailName, False) # false for gen
    #do2D(treeInfo, fileDir, outdir, background, var3Gen, var5Gen, treePassName, treeFailName, False) # false for gen
    do2D(treeInfo, fileDir, outdir, background, var1Reco, var2Reco, treePassName, treeFailName, True)
    do2D(treeInfo, fileDir, outdir, background, var3Reco, var4Reco, treePassName, treeFailName, True)
    #do2D(treeInfo, fileDir, outdir, background, var3Reco, var5Reco, treePassName, treeFailName, True)

'''
do2D(treeInfo, fileDir, outdir, zz, var1Gen, var2Gen, treePassName, treeFailName, False) # false for gen
do2D(treeInfo, fileDir, outdir, zz, var3Gen, var4Gen, treePassName, treeFailName, False) # false for gen
do2D(treeInfo, fileDir, outdir, zz, var3Gen, var5Gen, treePassName, treeFailName, False) # false for gen
do2D(treeInfo, fileDir, outdir, zz, var1Reco, var2Reco, treePassName, treeFailName, True)
do2D(treeInfo, fileDir, outdir, zz, var3Reco, var4Reco, treePassName, treeFailName, True)
do2D(treeInfo, fileDir, outdir, zz, var3Reco, var5Reco, treePassName, treeFailName, True)
'''
'''
doGenReco(treeInfo, sig, backgrounds, outdir, fileDir, treePassName, treeFailName, var1Gen, var1Reco, results)
doGenReco(treeInfo, sig, backgrounds, outdir, fileDir, treePassName, treeFailName, var2Gen, var2Reco, results)
doGenReco(treeInfo, sig, backgrounds, outdir, fileDir, treePassName, treeFailName, var3Gen, var3Reco, results)
doGenReco(treeInfo, sig, backgrounds, outdir, fileDir, treePassName, treeFailName, var4Gen, var4Reco, results)
doGenReco(treeInfo, sig, backgrounds, outdir, fileDir, treePassName, treeFailName, var5Gen, var5Reco, results)
'''

####################//####################//####################//####################

#do2D(fileDir, outdir, sig_unskimmed, var1Reco, var2Reco, treePassName, treeFailName, True)

#for background_unskimmed in backgrounds_unskimmed:
#    do2D(fileDir, outdir, background_unskimmed, var1Reco, var2Reco, treePassName, treeFailName, True)