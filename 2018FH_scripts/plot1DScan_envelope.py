# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
# This code is intended to be used as a NLL plotting macro when doing a discrete profiling, in the same style as the official plot1DScan.py tool    #                                                                               # 
# It generates a plot of negative log-likelihood (NLL) result using two different functions to parametrize the same distribution                    #  
# and the resulting NLL from the discrete profiling (envelope). It is intended to be used within the same CMSSW as combine harvester                #
#                                                                                                                                                   #
# Note: The user must edit the filenames and the outputfile for the plots accordingly                                                               #
#                                                                                                                                                   #
# Questions? Contact:                                                                                                                               #
# Author: Daina Leyva Pernia                                                                                                                        #
# e-mail: daina.leyva.pernia@desy.de                                                                                                                #
# prints 1-sigma uncertainty values in the top right corner of the plot.                                                                            #
# Date: October 16, 2023                                                                                                                            #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  


import ROOT
import numpy as np
import CombineHarvester.CombineTools.plotting as plot
import argparse
from functools import partial
from six.moves import range
import math

# Initialize ROOT in batch mode
ROOT.gROOT.SetBatch(True)

def RemoveGraphXDuplicates(graph):
    i = 0
    while i < graph.GetN() - 1:
        if graph.GetX()[i + 1] == graph.GetX()[i]:
            graph.RemovePoint(i + 1)
        else:
            i += 1
        
def output_name(mass):
    return "scan{mass}GeV"

def findMinCoord(graph, coord):
    min_x = graph.GetX()[0]
    min_y = graph.GetY()[0]
    for i in range(1, gr1.GetN()):
        x = gr1.GetX()[i]
        y = gr1.GetY()[i]
        if y < min_y:
            min_x = x
            min_y = y
    if coord == "x":
        return min_x
    elif coord == "y":
        return min_y
    else:
        print("findMinCoord ERROR, please select x or y")
        return 0


def Eval(obj, x, params):
    return obj.Eval(x[0])

parser = argparse.ArgumentParser()
#parser.add_argument('main', help='Main input file for the scan')
parser.add_argument('--mass', default='500', type=str, help='masspoint')
parser.add_argument('--SR', default='SR2', type=str, help='fit range in "SRX" format')
parser.add_argument('--ycut', type=float, default=6., help='Remove points with y > y-cut')
parser.add_argument('--y-max', type=float, default=8., help='y-axis maximum')
parser.add_argument('--x-min', type=float, default=-2., help='x-axis minimum')
parser.add_argument('--x-max', type=float, default=2., help='x-axis maximum')
#parser.add_argument('--plotall', type=bool, default=False,  help='plot the envelope and the other functions')
parser.add_argument('--plotall', action='store_true', default=False, help='plot the envelope and the other functions')
parser.add_argument('--output', '-o', help='output name without file extension', default="scan")
parser.add_argument('--logo', default='CMS')
#parser.add_argument('--logo-sub', default='Internal')
parser.add_argument('--logo-sub', default='Work in progress')
args = parser.parse_args()

print('--------------------------------------')
print(args.output)
print('--------------------------------------')

# Canvas and plotting misc.
canv = ROOT.TCanvas(args.output, args.output)

# Open the ROOT files and create TGraphs
file1 = ROOT.TFile("higgsCombine.Envelope."+str(args.SR)+".MultiDimFit.mH"+str(args.mass)+".root", "READ")
tree1 = file1.Get("limit")
tree1.Draw("2*(deltaNLL+nll+nll0):r")
gr1 = ROOT.TGraph(tree1.GetSelectedRows(), tree1.GetV2(), tree1.GetV1())
gr1.SetMarkerStyle(20)
gr1.SetMarkerColor(ROOT.kBlack)
gr1.SetLineColor(ROOT.kBlack)
gr1.SetLineWidth(2)
gr1.Sort()
gr1.SetTitle("")


file2 = ROOT.TFile("higgsCombine.fixed_pdf_0."+str(args.SR)+".MultiDimFit.mH"+str(args.mass)+".root", "READ")
tree2 = file2.Get("limit")
tree2.Draw("2*(deltaNLL+nll+nll0):r")
gr2 = ROOT.TGraph(tree2.GetSelectedRows(), tree2.GetV2(), tree2.GetV1())
gr2.SetMarkerStyle(20)
gr2.SetMarkerColor(ROOT.kBlue)
gr2.SetLineColor(ROOT.kBlue)
gr2.SetLineWidth(2)
gr2.Sort()

file3 = ROOT.TFile("higgsCombine.fixed_pdf_1."+str(args.SR)+".MultiDimFit.mH"+str(args.mass)+".root", "READ")
tree3 = file3.Get("limit")
tree3.Draw("2*(deltaNLL+nll+nll0):r")
gr3 = ROOT.TGraph(tree3.GetSelectedRows(), tree3.GetV2(), tree3.GetV1())
gr3.SetMarkerStyle(20)
gr3.SetMarkerColor(ROOT.kRed)
gr3.SetLineColor(ROOT.kRed)
gr3.SetLineWidth(2)
gr3.Sort()

# Shift Y axis to match the minimum value
minYEnvelope = min(gr1.GetY())
shiftValue = -minYEnvelope

for i in range(gr1.GetN()):
    gr1.GetY()[i] += shiftValue
for i in range(gr2.GetN()):
    gr2.GetY()[i] += shiftValue
for i in range(gr3.GetN()):
    gr3.GetY()[i] += shiftValue

# Set canvas and axis titles
canv.cd()
#canv.SetGrid()

# Remove high y points for better visualization
plot.RemoveGraphYAbove(gr1, args.ycut)
plot.RemoveGraphYAbove(gr2, args.ycut)
plot.RemoveGraphYAbove(gr3, args.ycut)

# Set Axis Title
gr1.GetXaxis().SetTitle("r")
gr1.GetYaxis().SetTitle("-2 \Delta lnL")

# Draw the ROOT graphs on the canvas
gr1.Draw('ALP')
if args.plotall == True:
    gr2.Draw('LP same')
    gr3.Draw('LP same')
    #gr1.Draw('LP same')
    

# Access the x and y-axis objects to change the title size
x_axis = gr1.GetXaxis()
y_axis = gr1.GetYaxis()
x_axis.SetTitleSize(0.04) 
y_axis.SetTitleSize(0.04) 

# Create a legend
legend = ROOT.TLegend(0.14, 0.76, 0.4, 0.6)  # (x1, y1, x2, y2) coordinates for the legend
if args.plotall == True:
    legend.AddEntry(gr1, "Envelope", "l")
    legend.AddEntry(gr2, "Quadratic", "l")
    legend.AddEntry(gr3, "Linear", "l")
legend.SetBorderSize(0)  # Remove the border of the legend
legend.Draw()

# Plotting misc
plot.DrawCMSLogo(canv, args.logo, args.logo_sub, 11, 0.045, 0.035, 1.2, cmsTextSize=0.5)
axishist = plot.GetAxisHist(canv)
axishist.SetMaximum(args.y_max)
new_min = axishist.GetXaxis().GetXmin()
new_max = axishist.GetXaxis().GetXmax()
axishist.GetXaxis().SetLimits(new_min, new_max)

# Lines for 1 and 2 sigma uncertainties
line = ROOT.TLine()
line.SetLineColor(16)
yvals = [1., 4.]
spline = ROOT.TSpline3("spline3", gr1)
func_method = partial(Eval, spline)
func_left  = ROOT.TF1('splinefnl', func_method, gr1.GetX()[0], findMinCoord(gr1,"x"), 1)
func_right = ROOT.TF1('splinefnr', func_method, findMinCoord(gr1,"x"), gr1.GetX()[gr1.GetN() - 1], 1)
for yval in yvals:
    # Horizontal lines
    plot.DrawHorizontalLine(canv, line, yval)
    # Vertical lines
    x_coordinate_left = func_left.GetX(yval)
    x_coordinate_right = func_right.GetX(yval)
    line.DrawLine(x_coordinate_left, 0, x_coordinate_left, yval)
    line.DrawLine(x_coordinate_right, 0, x_coordinate_right, yval)
    #print ("r = "+ str(findMinCoord(gr1,"x")) + " + " + str(x_coordinate_left)  + " - " + str(x_coordinate_right))

# Print r for min NLL and 1sigma uncertainty 
rmin = findMinCoord(gr1,"x")
r1sdown = rmin - func_left.GetX(1.)
r1sup = func_right.GetX(1.) - rmin
r2sdown = rmin - func_left.GetX(4.)
r2sup = func_right.GetX(4.) - rmin
print ("rmin = "+ str(rmin))
print ("1sdown = "+ str(r1sdown))
print ("1sup = "+ str(r1sup))
print ("2sdown = "+ str(r2sdown))
print ("2sup = "+ str(r2sup))
textfit = 'r = %.3f{}^{#plus %.3f}_{#minus %.3f}' % (rmin, r1sup, r1sdown)
text_box = ROOT.TPaveText(0.6, 0.75, 0.89, 0.85, "NDCNB")
text_box.SetFillStyle(0)
text_box.SetBorderSize(0)
text_box.SetTextFont(42)
text_box.AddText(textfit)
text_box.Draw()


# Save the plot as a PDF and PNG file
if args.plotall == True:
    name="plotsScan/1Dscan_"+ str(args.mass) +"_"+ str(args.SR)
else:
    name="plotsScan/1Dscan_envelope_"+ str(args.mass) +"_"+ str(args.SR)

canv.SaveAs(name + ".pdf")
canv.SaveAs(name + ".png")