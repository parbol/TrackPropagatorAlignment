import ROOT as r
import math


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
    const,mu,sigma = f.GetParameter(0), f.GetParameter(1), f.GetParameter(2)
    econst,emu,esigma = f.GetParError(0), f.GetParError(1), f.GetParError(2)
   
    return h, mu, sigma, emu

def plotResidual(t, side, tray, runumber, rutype, module, ptmin, ptmax, var):

    name = 'residual_' + var + '_' + str(side) + '_' + str(tray) + '_' + str(runumber) + '_' + str(rutype) + '_' + str(module) + '_' + str(ptmin) + '_' + str(ptmax)
    h, mu, sigma, emu = getHisto(t, side, tray, runumber, rutype, module, ptmin, ptmax, var)
    c = r.TCanvas('can_' + name)
    h.Draw()
    c.SaveAs(name + '.png')
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


    f = r.TFile('output.root')
    t = f.Get('hits')
    plotEtaPt(t)
    #mu, sigma, emu = plotResidual(t, 1, 0, 0, 0, 0, 0, 100, 'y')
    #makePtPlot(t, 1, 0, 0, 0, 0, 'x')
