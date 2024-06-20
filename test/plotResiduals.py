import ROOT as r
import math
from src.BTL import BTL
from src.BTLId import BTLId
from array import array
import numpy as np


def plotEtaPt(t):

    heta = r.TH1F('heta', '', 50, -1.5, 1.5)
    hpt = r.TH1F('hpt', '', 50, 26, 50)
    heta.GetXaxis().SetTitle('eta')
    hpt.GetXaxis().SetTitle('pt [GeV]')
    
    t.Project('heta', 'eta', '')
    t.Project('hpt', 'pt', '')
    c = r.TCanvas('can_eta')
    heta.GetYaxis().SetRangeUser(0, 1500)
    heta.SetStats(0)
    heta.Draw()
    
    c.SaveAs('eta.png')          
    c = r.TCanvas('can_pt')
    hpt.GetYaxis().SetRangeUser(0, 2500)
    hpt.SetStats(0)
    hpt.Draw()
    c.SaveAs('pt.png')

def getHisto(t, side, tray, runumber, rutype, module, ptmin, ptmax, var):

    varplot = 'x+xe'
    label = 'R_x [cm]'
    if var != 'x':
        varplot = 'y-ye'
        label = 'R_y [cm]'

    name = 'residual_' + var + '_' + str(side) + '_' + str(tray) + '_' + str(runumber) + '_' + str(rutype) + '_' + str(module) + '_' + str(ptmin) + '_' + str(ptmax)
    h = r.TH1F(name, '', 50, -1.5, 1.5)
    h.GetXaxis().SetTitle(label)
    cut1 = 'side == ' + str(side) + ' && tray == ' + str(tray) + ' && RUNumber == ' + str(runumber) + ' && RUType == ' + str(rutype) + ' && module == ' + str(module)
    cut2 = 'pt >= ' + str(ptmin) + ' && pt < ' + str(ptmax)
    cut = cut1 + ' && ' + cut2
    t.Project(name, varplot, cut)
    # Fit histogram with root #
    h.Fit('gaus','','', -1.0, 1.0)
    # Get Root Fit and Goodness of Fit Parameters #
    f = h.GetFunction('gaus')
    if h.GetEntries() == 0:
        return h, 0, -1, 0
    const,mu,sigma = f.GetParameter(0), f.GetParameter(1), f.GetParameter(2)
    econst,emu,esigma = f.GetParError(0), f.GetParError(1), f.GetParError(2)
   
    return h, mu, sigma, emu

def plotResidual(t, side, tray, runumber, rutype, module, ptmin, ptmax, var, tag):

    name = 'residual_' + var + '_' + str(side) + '_' + str(tray) + '_' + str(runumber) + '_' + str(rutype) + '_' + str(module) + '_' + str(ptmin) + '_' + str(ptmax)
    h, mu, sigma, emu = getHisto(t, side, tray, runumber, rutype, module, ptmin, ptmax, var)
    c = r.TCanvas('can_' + name)
    h.Draw()
    c.SaveAs(name + '_' + tag + '.png')
    return mu, sigma, emu

def makePtPlot(t, side, tray, runnumber, rutype, module, var):
        
        name = 'residualPt_' + var + '_' + str(side) + '_' + str(tray) + '_' + str(runnumber) + '_' + str(rutype) + '_' + str(module)
        h1, mu1, sigma1, emu1 = getHisto(t, side, tray, runnumber, rutype, module, 2.0, 5.0, var)
        h2, mu2, sigma2, emu2 = getHisto(t, side, tray, runnumber, rutype, module, 5.0, 10.0, var)
        h3, mu3, sigma3, emu3 = getHisto(t, side, tray, runnumber, rutype, module, 10.0, 20.0, var)
        h4, mu4, sigma4, emu4 = getHisto(t, side, tray, runnumber, rutype, module, 20.0, 35.0, var)
        h5, mu5, sigma5, emu5 = getHisto(t, side, tray, runnumber, rutype, module, 35.0, 50.0, var)
        c = r.TCanvas('can_' + name)
        h1.Draw()
        h2.Draw('SAME')
        h3.Draw('SAME')
        h4.Draw('SAME')
        h5.Draw('SAME')
        c.SaveAs(name + '.png')

if __name__=='__main__':

    #Configuring the BTL
    R = 114.8
    TrayLength = 300.0
    TrayWidth = 2.0*R*np.sin(3.0*np.pi/180.0)
    TrayStartZ = 1.0
    TrayStartPhi = 5.0*np.pi/180.0
    RULength = 45.0
    ModuleLength = 5.4
    ModuleWidth = 4.0
    rphi_error = 0.1
    z_error = 0.4
    t_error = 0.1
    btl = BTL(R, TrayLength, TrayWidth, TrayStartZ, TrayStartPhi, RULength, ModuleLength, ModuleWidth, rphi_error, z_error, t_error, 9.4)
    
    files = ['output1M.root', 'output0.75M.root', 'output0.5M.root', 'output0.25M.root']
    tags = ['1M', '0.75M', '0.5M', '0.25M']
    #files = ['output1M.root', 'output0.75M.root']
    #tags = ['1M', '0.75M']
    colors = [r.kRed, r.kBlue, r.kGreen, r.kBlack]
    legx = r.TLegend(0.75, 0.7, 0.9, 0.85)
    legx.SetFillColor(r.kWhite)
    legx.SetTextFont(42)
    legx.SetTextSize(0.04)
    legx.SetLineWidth(0)
    legx.SetBorderSize(0)
    legy = r.TLegend(0.75, 0.7, 0.9, 0.85)
    legy.SetFillColor(r.kWhite)
    legy.SetTextFont(42)
    legy.SetTextSize(0.04)
    legy.SetLineWidth(0)
    legy.SetBorderSize(0)
    auxX = r.TH1F('auxX', '', 1, 0, 350)
    auxX.GetXaxis().SetTitle("Module Z [cm]")
    auxX.GetYaxis().SetTitle("X uncertainty [cm]")
    auxX.SetStats(0)
    auxX.SetTitle('X uncertainty vs. module Z')
    auxX.GetYaxis().SetRangeUser(0, 0.5)
    auxY = r.TH1F('auxY', '', 1, 0, 350)
    auxY.GetXaxis().SetTitle("Module Z [cm]")
    auxY.GetYaxis().SetTitle("Y uncertainty [cm]")
    auxY.SetStats(0)
    auxY.SetTitle('Y uncertainty vs. module Z')
    auxY.GetYaxis().SetRangeUser(0, 0.3)
    grxs = []
    grys = []

    for i, fg in enumerate(files):

        f = r.TFile(fg)
        tag = tags[i]
        ccolor = colors[i]
        t = f.Get('hits')
    
        #plotEtaPt(t)
    
        #Make resolution study
        ruType = [0, 1, 2]
        ruNumber = [0, 1]
        module = [0, 3, 6, 9, 12, 15, 18, 21]

        zx = array('d')
        zy = array('d')
        ex = array('d')
        ey = array('d')

        for ruT in ruType:
            for ruN in ruNumber:
                for mod in module:
                    id = BTLId()
                    id.setSide(1)
                    id.setTray(12)
                    id.setRU(ruT, ruN)
                    id.setModule(mod)
                    modulaco = btl.btlNominal.getBTLModule(id)
                    mux, sigmax, emux = plotResidual(t, 1, 12, ruN, ruT, mod, 0, 100, 'x', tag)
                    muy, sigmay, emuy = plotResidual(t, 1, 12, ruN, ruT, mod, 0, 100, 'y', tag)
                    if sigmax > 0:
                        zx.append(modulaco.r[2])
                        ex.append(emux)
                    if sigmay > 0:
                        zy.append(modulaco.r[2])
                        ey.append(emuy)

        
        grx = r.TGraph(len(zx), zx, ex)
        gry = r.TGraph(len(zy), zy, ey)
        legx.AddEntry(grx, tag, "P")
        legy.AddEntry(gry, tag, "P")
        grx.GetXaxis().SetTitle("Module Z [cm]")
        grx.GetYaxis().SetTitle("X uncertainty [cm]")
        gry.GetXaxis().SetTitle("Module Z [cm]")
        gry.GetYaxis().SetTitle("Y uncertainty [cm]")
        grx.SetMarkerColor(ccolor)
        grx.SetMarkerSize(1)
        grx.SetMarkerStyle(20)
        gry.SetMarkerColor(ccolor)
        gry.SetMarkerSize(1)
        gry.SetMarkerStyle(20)
        grx.SetTitle('Uncertainty vs. module Z')
        gry.SetTitle('Uncertainty vs. module Z')
        grxs.append(grx)
        grys.append(gry)
    
    c = r.TCanvas('canx')
    auxX.Draw()
    for gr in grxs:
        gr.Draw('P')
    legx.Draw()
    c.SaveAs('plotx.png')
    cy = r.TCanvas('cany')
    auxY.Draw()
    for gr in grys:
        gr.Draw('P')
    legy.Draw()
    cy.SaveAs('ploty.png')    


    