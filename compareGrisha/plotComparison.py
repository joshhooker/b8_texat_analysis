from __future__ import print_function
from ROOT import TCanvas, TGraph, TGraphErrors, TMultiGraph
from ROOT import gROOT
from math import sin
from array import array

import numpy as np

mg = TMultiGraph()

c1 = TCanvas( 'c1', 'A Simple Graph Example', 200, 10, 700, 500 )

dataset = np.loadtxt("ourData.dat")
grishaData = np.loadtxt("grishaData.dat")

n = len(dataset[:, 0])
nGrisha = len(grishaData[:, 0])

cmEnergy = array('d')
excitationEnergy = array('d')
cs = array('d')
csX = array('d')
csErr = array('d')

excitationEnergyGrisha = array('d')
csGrisha = array('d')

num = 0
for i in range(n):
    if(dataset[i, 0]) < 0.97:
        continue
    excitationEnergy.append(dataset[i, 0] + 1.3)
    cs.append(dataset[i, 2]*1000)
    csX.append(0)
    csErr.append(dataset[i, 3]*1000)
    num += 1

numGrisha = 0
for i in range(nGrisha):
    excitationEnergyGrisha.append(grishaData[i, 0])
    csGrisha.append(grishaData[i, 1])
    numGrisha += 1

gr = TGraphErrors(num, excitationEnergy, cs, csX, csErr)
gr1 = TGraph(numGrisha, excitationEnergyGrisha, csGrisha)

gr1.SetLineColor(2)

# gr1.Draw()
# gr.Draw("psame")

gr.Draw()

# mg.Add(gr)
# mg.Add(gr1)

# mg.Draw("ae")

c1.Modified()
c1.Update()

raw_input()