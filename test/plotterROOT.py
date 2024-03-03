import ROOT as r
import math







if __name__=='__main__':


    hdist = r.TH1F("hdist", "", 100, 0, 0.04)
    hdist.GetXaxis().SetTitle('d [cm]')
    ht = r.TH1F("ht", "", 100, 0, 0.005)
    ht.GetXaxis().SetTitle('t [ns]')
    hdistpt = r.TProfile("hdistpt", "", 20, 0, 50, 0, 0.04)
    hdistpt.GetXaxis().SetTitle('pt [GeV]')
    hdistpt.GetYaxis().SetTitle('d [cm]')
    htpt = r.TProfile("htpt", "", 20, 0, 50, 0, 0.01)
    htpt.GetXaxis().SetTitle('pt [GeV]')
    htpt.GetYaxis().SetTitle('t [ns]')
    hdisteta = r.TProfile("hdisteta", "", 20, 1.5, 3.0, 0, 0.04)
    hdisteta.GetXaxis().SetTitle('eta')
    hdisteta.GetYaxis().SetTitle('d [cm]')
    hteta = r.TProfile("hteta", "", 20, 1.5, 3.0, 0, 0.01)
    hteta.GetXaxis().SetTitle('eta')
    hteta.GetYaxis().SetTitle('t [ns]')

    for l in open('datos.txt').readlines():
        
        a = l.split()
        xi = float(a[0])
        yi = float(a[1])
        zi = float(a[2])
        ti = float(a[3])
        xin = float(a[4])
        yin = float(a[5])
        zin = float(a[6])
        tin = float(a[7])
        pt = float(a[8])
        eta = float(a[9])
        d = math.sqrt((xi-xin)**2+(yi-yin)**2+(zi-zin)**2)
        hdist.Fill(d)
        ht.Fill(ti-tin)
        hdistpt.Fill(pt, d)
        hdisteta.Fill(eta, d)
        htpt.Fill(pt, ti-tin)
        hteta.Fill(eta, ti-tin)

    canD = r.TCanvas('canD', 'canD')
    canD.cd()
    hdist.Draw("HIST")
    canD.SaveAs('Distance.png')
    canT = r.TCanvas('canT', 'canT')
    canT.cd()
    ht.Draw("HIST")
    canT.SaveAs('Time.png')
    canptD = r.TCanvas('canptD', 'canptD')
    canptD.cd()
    hdistpt.Draw()
    canptD.SaveAs('ptD.png')
    canetaD = r.TCanvas('canetaD', 'canetaD')
    canetaD.cd()
    hdisteta.Draw()
    canetaD.SaveAs('etaD.png')
    canptt = r.TCanvas('canptt', 'canptt')
    canptt.cd()
    htpt.Draw()
    canptt.SaveAs('ptT.png')
    canetat = r.TCanvas('canetat', 'canetat')
    canetat.cd()
    hteta.Draw()
    canetat.SaveAs('etaT.png')