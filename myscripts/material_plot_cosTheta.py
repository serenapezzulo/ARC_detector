# Code taken from Brieuc for plotting material scan using Gaudi
# https://raw.githubusercontent.com/BrieucF/LAr_scripts/main/geometry/material_plot.py

from __future__ import print_function
import argparse

from plotstyle import FCCStyle

import ROOT
ROOT.gROOT.SetBatch(ROOT.kTRUE)

def main():
    parser = argparse.ArgumentParser(description='Material Plotter')
    parser.add_argument('--fname', "-f", dest='fname', default="out_material_scan.root", type=str, help="name of file to read")
    parser.add_argument('--angleMax', "-M", dest='angleMax', default=1., type=float, help="maximum angle")
    parser.add_argument('--angleMin', "-m", dest='angleMin', default=0, type=float, help="minimum angle")
    parser.add_argument('--angleBin', "-b", dest='angleBin', default=0.02, type=float, help="angle bin width")
    parser.add_argument('--angleDef', "-d", dest='angleDef', default='cos(#theta)', type=str, help="angle type")
    args = parser.parse_args()

    f = ROOT.TFile.Open(args.fname, "read")
    tree = f.Get("materials")
    nbins = int(tree.GetEntries())
    histDict = {}


    histDict['Rohacell300WP'] = {
                    "x0": ROOT.TH1F("", "",     nbins, args.angleMin, args.angleMax),
                    "lambda": ROOT.TH1F("", "", nbins, args.angleMin, args.angleMax),
                    "depth": ROOT.TH1F("", "",  nbins, args.angleMin, args.angleMax)
                }
    histDict['C4F10_PFRICH'] = {
                    "x0": ROOT.TH1F("", "",     nbins, args.angleMin, args.angleMax),
                    "lambda": ROOT.TH1F("", "", nbins, args.angleMin, args.angleMax),
                    "depth": ROOT.TH1F("", "",  nbins, args.angleMin, args.angleMax)
                }
    histDict['Aerogel_PFRICH'] = {
                    "x0": ROOT.TH1F("", "",     nbins, args.angleMin, args.angleMax),
                    "lambda": ROOT.TH1F("", "", nbins, args.angleMin, args.angleMax),
                    "depth": ROOT.TH1F("", "",  nbins, args.angleMin, args.angleMax)
                }
    histDict['CarbonFibStr'] = {
                    "x0": ROOT.TH1F("", "",     nbins, args.angleMin, args.angleMax),
                    "lambda": ROOT.TH1F("", "", nbins, args.angleMin, args.angleMax),
                    "depth": ROOT.TH1F("", "",  nbins, args.angleMin, args.angleMax)
                }
    histDict['SiliconOptical'] = {
                    "x0": ROOT.TH1F("", "",     nbins, args.angleMin, args.angleMax),
                    "lambda": ROOT.TH1F("", "", nbins, args.angleMin, args.angleMax),
                    "depth": ROOT.TH1F("", "",  nbins, args.angleMin, args.angleMax)
                }
    # go through the eta bins and fill the histograms in the histDict, skipping air
    # keys in the histDict are the material names
    for angleBin, entry in enumerate(tree):
        nMat = entry.nMaterials
        for i in range(nMat):
            if entry.material.at(i) == "Air": continue
            if entry.material.at(i) not in histDict.keys():
                histDict[entry.material.at(i)] = {
                    "x0": ROOT.TH1F("", "",     nbins, args.angleMin, args.angleMax),
                    "lambda": ROOT.TH1F("", "", nbins, args.angleMin, args.angleMax),
                    "depth": ROOT.TH1F("", "",  nbins, args.angleMin, args.angleMax)
                }
            hs = histDict[entry.material.at(i)]
            hs["x0"].SetBinContent(angleBin+1, hs["x0"].GetBinContent(angleBin+1) + entry.nX0.at(i))
            hs["lambda"].SetBinContent(angleBin+1, hs["lambda"].GetBinContent(angleBin+1) + entry.nLambda.at(i))
            hs["depth"].SetBinContent(angleBin+1, hs["depth"].GetBinContent(angleBin+1) + entry.matDepth.at(i))

    axis_titles = ["Number of X_{0}", "Number of #lambda", "Material depth [cm]"]

    # This loop does the drawing, sets the style and saves the pdf files
    for plot, title in zip(["x0", "lambda", "depth"], axis_titles):
        legend = ROOT.TLegend(.75, .75, .94, .94)
        legend.SetLineColor(0)
        ths = ROOT.THStack()
        for i, material in enumerate(histDict.keys()):
            linecolor = 1
            if i >= len(FCCStyle.fillcolors):
                i = i%len(FCCStyle.fillcolors)

            fillcolor = FCCStyle.fillcolors[i]
            histDict[material][plot].SetLineColor(linecolor)
            histDict[material][plot].SetFillColor(fillcolor)
            histDict[material][plot].SetLineWidth(1)
            histDict[material][plot].SetFillStyle(1001)

            ths.Add(histDict[material][plot])
            legend.AddEntry(histDict[material][plot], material, "f")

        ths.SetMaximum(1.5 * ths.GetMaximum())
        cv = ROOT.TCanvas()
        ths.Draw()
        ths.GetXaxis().SetTitle(args.angleDef)
        ths.GetYaxis().SetTitle(title)

        legend.Draw()
        cv.Print(plot + ".root")
        cv.Print(plot + ".png")

        #ths.GetXaxis().SetRangeUser(0, args.angleMax)
        #cv.Print(plot + "pos.pdf")
        #cv.Print(plot + "pos.png")

if __name__ == "__main__":
    FCCStyle.initialize()
    main()


