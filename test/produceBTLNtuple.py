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
    layers = np.linspace(3, 114, num=20)
    layersz = np.linspace(130, 270, 5)
    sigma_rphi = 0.312/np.sqrt(12.00)
    sigma_z = 0.3
    tracker = Tracker(layers, layersz, sigma_rphi, sigma_z, 220.0, [0,0,0])
   
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
    z_error = 0.1
    t_error = 0.1
    btl = BTL(R, TrayLength, TrayWidth, TrayStartZ, TrayStartPhi, RULength, ModuleLength, ModuleWidth, rphi_error, z_error, t_error, 9.4)
    btl.btlReal.readGeometry('btlUpdatedGeometry.txt')
    
    pt = array('f', [0])
    phi = array('f', [0])
    eta = array('f', [0])
    charge = array('i', [0])
    type = array('i', [0])
    nEvent = array('i', [0])
    side = array('i', [0])
    tray = array('i', [0])
    RUType = array('i', [0])
    RUNumber = array('i', [0])
    module = array('i', [0])
    xe = array('f', [0])
    ye = array('f', [0])
    ze = array('f', [0])
    xge = array('f', [0])
    yge = array('f', [0])
    zge = array('f', [0])
    x = array('f', [0])
    y = array('f', [0])
    z = array('f', [0])
    xg = array('f', [0])
    yg = array('f', [0])
    zg = array('f', [0])
    
    tree = r.TTree("hits", "hits")
    tree.Branch('pt', pt, 'pt/F')
    tree.Branch('phi', phi, 'phi/F')
    tree.Branch('eta', eta, 'eta/F')
    tree.Branch('charge', charge, 'charge/I')
    tree.Branch('type', type, 'type/I')
    tree.Branch('nEvent', nEvent, 'nEvent/I')
    tree.Branch('side', side, 'side/I')
    tree.Branch('tray', tray, 'tray/I')
    tree.Branch('RUType', RUType, 'RUType/I')
    tree.Branch('RUNumber', RUNumber, 'RUNumber/I')
    tree.Branch('module', module, 'module/I')
    tree.Branch('xe', xe, 'xe/F')
    tree.Branch('ye', ye, 'ye/F')
    tree.Branch('ze', ze, 'ze/F')
    tree.Branch('xge', xge, 'xge/F')
    tree.Branch('yge', yge, 'yge/F')
    tree.Branch('zge', zge, 'zge/F')
    tree.Branch('x', x, 'x/F')
    tree.Branch('y', y, 'y/F')
    tree.Branch('z', z, 'z/F')
    tree.Branch('xg', xg, 'xg/F')
    tree.Branch('yg', yg, 'yg/F')
    tree.Branch('zg', zg, 'zg/F')

    counter = 0
    alist = []
    while counter < int(opts.nTracks):

        #phi_ = np.random.uniform(np.pi/2.0 - 0, 2.0*np.pi)
        phi_ = np.random.uniform(5.0*np.pi/180.0, 5.0*np.pi/180.0 + 4.0*np.pi/180.0)
        eta_ = np.random.uniform(0, 0.2)
        pt_ = np.random.uniform(10.0, 50.0)
        charge_ = int(np.sign(np.random.uniform(-1.0, 1.0)))
        track = Track(0, 0, phi_, eta_, pt_, charge_)
        tracker.fullMeasurement(track)
        valid = btl.fullMeasurement(track)
        if not valid:
            continue
        pt[0] = pt_
        phi[0] = phi_
        eta[0] = eta_
        nEvent[0] = counter
        charge[0] = charge_
       
        type[0] = 0
        xe[0] = track.lxi[len(track.lxi)-1]
        ye[0] = track.lyi[len(track.lyi)-1]
        ze[0] = track.lzi[len(track.lzi)-1]
        x[0] = track.lxm[len(track.lxm)-1]
        y[0] = track.lym[len(track.lym)-1]
        z[0] = track.lzm[len(track.lzm)-1]
        xge[0] = track.xi[len(track.xi)-1]
        yge[0] = track.yi[len(track.yi)-1]
        zge[0] = track.zi[len(track.zi)-1]
        xg[0] = track.xm[len(track.xm)-1]
        yg[0] = track.ym[len(track.ym)-1]
        zg[0] = track.zm[len(track.zm)-1]   
        side[0] = (track.subdet[len(track.subdet)-1])[0]
        tray[0] = track.subdet[len(track.subdet)-1][1]
        RUType[0] = track.subdet[len(track.subdet)-1][2]
        RUNumber[0] = track.subdet[len(track.subdet)-1][3]
        module[0] = track.subdet[len(track.subdet)-1][4]

        tree.Fill()
        counter = counter + 1
           
    f = r.TFile(opts.outputFile, "RECREATE")
    f.WriteObject(tree, "hits")
    f.Close()

    



 
