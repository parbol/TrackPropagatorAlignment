import numpy as np
import sys
import optparse
import ROOT as r
from array import array

from src.Tracker import Tracker
from src.Track import Track
from src.BTL import BTL



if __name__ == "__main__":

   
    parser = optparse.OptionParser(usage='usage: %prog [options] path', version='%prog 1.0')
    parser.add_option('-n', '--ntrack', action='store', type=int, dest='nTracks', default=10, help='Number of tracks')
    parser.add_option('-m', '--misalignedGeom', action='store', type='string', dest='misalignedFile', default='misaligned.csv', help='Name of misaligned geometry')
    parser.add_option('-o', '--output', action='store', type='string', dest='outputFile', default='output.root', help='Name of output file.')


    (opts, args) = parser.parse_args()

    #configuring the tracker
    layers = np.linspace(1, 100, 20)
    layersz = np.linspace(130, 270, 5)
    sigma_rphi = 0.01
    sigma_z = 0.01
    tracker = Tracker(layers, layersz, sigma_rphi, sigma_z, 220.0, [0,0,0])
   
    #Configuring the BTL
    R = 120.0
    TrayLength = 300.0
    TrayWidth = 2.0*R*np.sin(3.0*np.pi/180.0)
    TrayStartZ = 1.0
    TrayStartPhi = 5.0*np.pi/180.0
    RULength = 45.0
    ModuleLength = 5.4
    ModuleWidth = 4.0
    rphi_error = 0.1
    z_error = 0.1
    t_error = 0.1
    btlAligned = BTL(R, TrayLength, TrayWidth, TrayStartZ, TrayStartPhi, RULength, ModuleLength, ModuleWidth, rphi_error, z_error, t_error, 9.4)
    btlMisaligned = BTL(R, TrayLength, TrayWidth, TrayStartZ, TrayStartPhi, RULength, ModuleLength, ModuleWidth, rphi_error, z_error, t_error, 9.4)
    btlMisaligned.readGeometry('btlUpdatedGeometry.txt')
    

    pt = array('f', [0])
    phi = array('f', [0])
    eta = array('f', [0])
    charge = array('i', [0])
    type = array('i', [0])
    nEvent = array('i', [0])
    sidenom = array('i', [0])
    traynom = array('i', [0])
    RUTypenom = array('i', [0])
    RUNumbernom = array('i', [0])
    modulenom = array('i', [0])
    sidemis = array('i', [0])
    traymis = array('i', [0])
    RUTypemis = array('i', [0])
    RUNumbermis = array('i', [0])
    modulemis = array('i', [0])
    xenom = array('f', [0])
    yenom = array('f', [0])
    zenom = array('f', [0])
    xemis = array('f', [0])
    yemis = array('f', [0])
    zemis = array('f', [0])
    xmnom = array('f', [0])
    ymnom = array('f', [0])
    zmnom = array('f', [0])
    xmmis = array('f', [0])
    ymmis = array('f', [0])
    zmmis = array('f', [0])
    xgmnom = array('f', [0])
    ygmnom = array('f', [0])
    zgmnom = array('f', [0])
    tree = r.TTree("hits", "hits")
    tree.Branch('pt', pt, 'pt/F')
    tree.Branch('phi', phi, 'phi/F')
    tree.Branch('eta', eta, 'eta/F')
    tree.Branch('charge', charge, 'charge/I')
    tree.Branch('type', type, 'type/I')
    tree.Branch('nEvent', nEvent, 'nEvent/I')
    tree.Branch('sidenom', sidenom, 'sidenom/I')
    tree.Branch('traynom', traynom, 'traynom/I')
    tree.Branch('RUTypenom', RUTypenom, 'RUTypenom/I')
    tree.Branch('RUNumbernom', RUNumbernom, 'RUNumbernom/I')
    tree.Branch('modulenom', modulenom, 'modulenom/I')
    tree.Branch('sidemis', sidemis, 'sidemis/I')
    tree.Branch('traymis', traymis, 'traymis/I')
    tree.Branch('RUTypemis', RUTypemis, 'RUTypemis/I')
    tree.Branch('RUNumbermis', RUNumbermis, 'RUNumbermis/I')
    tree.Branch('modulemis', modulemis, 'modulemis/I')
    tree.Branch('xenom', xenom, 'xenom/F')
    tree.Branch('yenom', yenom, 'yenom/F')
    tree.Branch('zenom', zenom, 'zenom/F')
    tree.Branch('xemis', xemis, 'xemis/F')
    tree.Branch('yemis', yemis, 'yemis/F')
    tree.Branch('zemis', zemis, 'zemis/F')
    tree.Branch('xmnom', xmnom, 'xmnom/F')
    tree.Branch('ymnom', ymnom, 'ymnom/F')
    tree.Branch('zmnom', zmnom, 'zmnom/F')
    tree.Branch('xmmis', xmmis, 'xmmis/F')
    tree.Branch('ymmis', ymmis, 'ymmis/F')
    tree.Branch('zmmis', zmmis, 'zmmis/F')
    tree.Branch('xgmnom', xgmnom, 'xgmnom/F')
    tree.Branch('ygmnom', ygmnom, 'ygmnom/F')
    tree.Branch('zgmnom', zgmnom, 'zgmnom/F')
    counter = 0
    alist = []
    while counter < int(opts.nTracks):

        phi_ = np.random.uniform(0, 2.0*np.pi)
        eta_ = np.random.uniform(-1.6, 1.6)
        pt_ = np.random.uniform(2.0, 20.0)
        charge_ = int(np.sign(np.random.uniform(-1.0, 1.0)))
        track = Track(0, 0, phi_, eta_, pt_, charge_)
        track2 = Track(0, 0, phi_, eta_, pt_, charge_)
        tracker.fullMeasurement(track)
        validAligned = btlAligned.fullMeasurement(track)
        tracker.fullMeasurement(track2)
        validMisaligned = btlMisaligned.fullMeasurement(track2)
        pt[0] = pt_
        phi[0] = phi_
        eta[0] = eta_
        nEvent[0] = counter
        charge[0] = charge_
        if not validAligned and not validMisaligned:
            continue
        if validAligned:
            type[0] = 0
            xenom[0] = track.lxi[len(track.lxi)-1]
            yenom[0] = track.lyi[len(track.lyi)-1]
            zenom[0] = track.lzi[len(track.lzi)-1]
            xmnom[0] = track.lxm[len(track.lxm)-1]
            ymnom[0] = track.lym[len(track.lym)-1]
            zmnom[0] = track.lzm[len(track.lzm)-1]
            xgmnom[0] = track.xi[len(track.xi)-1]
            ygmnom[0] = track.yi[len(track.yi)-1]
            zgmnom[0] = track.zi[len(track.zi)-1]   
            sidenom[0] = (track.subdet[len(track.subdet)-1])[0]
            traynom[0] = track.subdet[len(track.subdet)-1][1]
            RUTypenom[0] = track.subdet[len(track.subdet)-1][2]
            RUNumbernom[0] = track.subdet[len(track.subdet)-1][3]
            modulenom[0] = track.subdet[len(track.subdet)-1][4]

        if validMisaligned:
            type[0] = 1
            xemis[0] = track2.lxi[len(track.lxi)-1]
            yemis[0] = track2.lyi[len(track.lyi)-1]
            zemis[0] = track2.lzi[len(track.lzi)-1]
            xmmis[0] = track2.lxm[len(track.lxm)-1]
            ymmis[0] = track2.lym[len(track.lym)-1]
            zmmis[0] = track2.lzm[len(track.lzm)-1]
            sidemis[0] = track2.subdet[len(track2.subdet)-1][0]
            traymis[0] = track2.subdet[len(track2.subdet)-1][1]
            RUTypemis[0] = track2.subdet[len(track2.subdet)-1][2]
            RUNumbermis[0] = track2.subdet[len(track2.subdet)-1][3]
            modulemis[0] = track2.subdet[len(track2.subdet)-1][4]

        if validAligned and validMisaligned:
            type[0] = 3 
        
        tree.Fill()
        counter = counter + 1
           
    f = r.TFile(opts.outputFile, "RECREATE")
    f.WriteObject(tree, "hits")
    f.Close()

    



 
