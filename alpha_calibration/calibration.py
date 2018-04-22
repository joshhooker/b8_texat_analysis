from __future__ import absolute_import, division, print_function, unicode_literals

import ROOT

tfile = ROOT.TFile.Open("../spectra_siCalibration.root")

fileDict = {}
for h in tfile.GetListOfKeys():
  h = h.ReadObj()
  #print h.ClassName(), h.GetName()
  if(h.ClassName()=='TH1F' or h.ClassName()=='TH1D' or h.ClassName()=='TH2F' or h.ClassName()=='TH2D'):
    fileDict[h.GetName()] = h

# print(fileDict)
peaks = ROOT.TSpectrum(5)
nfound = peaks.Search(fileDict['siEForward_d1_q1_ch5'], 2, "", 0.10)
print(nfound)
print(peaks)
# total = ROOT.TF1("total", "gaus(0) + gaus(3) + gaus(6) + gaus(9)")
# total.SetParameter(0, 10); total.SetParameter(1, 400); total.SetParameter(2, 5)
# total.SetParameter(3, 10); total.SetParameter(4, 600); total.SetParameter(5, 5)
# total.SetParameter(6, 10); total.SetParameter(7, 800); total.SetParameter(8, 5)
# total.SetParameter(9, 10); total.SetParameter(10, 1000); total.SetParameter(11, 5)

# fileDict['siEForward_d1_q1_ch5'].Fit("total", "MLL")
fileDict['siEForward_d1_q1_ch5'].Draw()

input()