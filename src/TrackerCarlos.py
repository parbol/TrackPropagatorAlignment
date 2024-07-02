from src.Plane import Plane
from src.Module import Module
from src.BTLTray import BTLTray
from src.EulerRotation import EulerRotation
from src.BTLId import BTLId
import numpy as np
import sys


class LayerCarlos:

    def __init__(self, R, ModuleLength, N, ModuleStartPhi, Ly, rphiError, zError, tError):

        self.R = R
        self.ModuleStartPhi = ModuleStartPhi
        self.ModuleLength = ModuleLength
        self.N = N
        self.rphi_error = rphiError
        self.z_error = zError
        self.t_error = tError
        
        self.phiStep = 2.0*np.pi / N

        self.Lx = np.tan(self.phiStep/2.0) * self.R * 2.0
        self.Ly = Ly
        self.modules = []
        for i in range(0, N):
            phi = self.ModuleStartPhi + self.phiStep / 2.0 + i * self.phiStep
            x = self.R * np.cos(self.phiStep/2.0) * np.cos(phi)
            y = self.R * np.cos(self.phiStep/2.0) * np.sin(phi)
            z = 0.0
            vx = np.asarray([np.sin(phi), -np.cos(phi), 0.0])
            vy = np.asarray([0.0, 0.0, 1.0])
            vz = np.asarray([np.cos(phi), np.sin(phi), 0.0])
            euler = EulerRotation()
            euler.setFromVectors(vx, vy, vz)
            mod = Module(np.asarray([x, y, z]), euler, self.Lx, self.Ly)
            self.modules.append(mod)


    def getPhi(self, x, y):
        
        sinphi = (y/np.sqrt(x**2+y**2)) 
        cosphi = (x/np.sqrt(x**2+y**2))
        if cosphi > 1.0:
            cosphi = 1.0
        if cosphi < -1.0:
            cosphi = -1.0
        if sinphi >= 0:
            return np.arccos(cosphi)
        else:
            return -np.arccos(cosphi) + 2.0*np.pi

    def getClosestModule(self, phi):
        
        index = -1
        minPhi = 9999999.0
        for j, i in enumerate(self.modules):
            modphi = self.getPhi(i.x[0], i.x[1])
            if np.abs(modphi-phi) < minPhi:
                minPhi = np.abs(modphi-phi)
                index = j
        return index

    
    def intersection(self, track):

        phi = track.phi
        eta = track.eta
        
        moduleNumber = self.getClosestModule(phi)
        navigationWindow = 2
        navigationList = []
        for i in range(0, navigationWindow+1):
            navigationList.append(i)
            if i != 0:
                navigationList.append(-i)
        for i in navigationList:
            checkModuleNumber = (moduleNumber + i) % self.N
            valid, x, y, z, t = self.modules[checkModuleNumber].intersection(track)
            if valid:
                return True, x, y, z, t
        return False, 0, 0, 0, 0


    def draw(self, ax1, ax2, ax3, ax4, t):

        for m in self.modules:
            m.drawModule(ax1, ax2, ax3, ax4, t, alpha=0.1)


class TrackerCarlos:
    
    def __init__(self, R, ModuleLength, N, ModuleStartPhi, Ly, rphiError, zError, tError):

        ####################################################################################################
        #                                   Representation of a Carlos Tracker                             #
        # R: Radius at which the trays are located
        # ModuleLength: Length in Z of the tray
        # TrayWidth: Width (rphi) of the tray
        # ModuleStartPhi: Start of the tray with respect to 0 in phi
        # rphi_error: uncertainty in rphi
        # z_error: undertainty in z
        # t_error: time uncertainty
        #####################################################################################################
        self.R = R
        self.ModuleStartPhi = ModuleStartPhi
        self.ModuleLength = ModuleLength
        self.N = N
        self.rphi_error = rphiError
        self.z_error = zError
        self.t_error = tError
        self.Ly = Ly
        
        self.layers = []
        for i, r in enumerate(self.R):
            lay = LayerCarlos(r, ModuleLength, N[i], ModuleStartPhi, Ly, rphiError, zError, tError)
            self.layers.append(lay)

    
    def intersection(self, track):

        xt = []
        yt = []
        zt = []
        tt = []
        for l in self.layers:
            valid, x, y, z, t = l.intersection(track)
            if valid:
                xt.append(x)
                yt.append(y)
                zt.append(z)
                tt.append(t)
        if len(xt) != 0:
            return True, xt, yt, zt, tt
        else:
            return False, [], [], [], []
    
        
    def fullMeasurement(self, track):
 
        #Intersection
        status, xt, yt, zt, tt = self.intersection(track)
        if not status:
            return False
        x = np.asarray(xt)
        y = np.asarray(yt)
        z = np.asarray(zt)
        t = np.asarray(tt)
        track.xi = np.concatenate((track.xi, x), axis=0)
        track.yi = np.concatenate((track.yi, y), axis=0)
        track.zi = np.concatenate((track.zi, z), axis=0)
        track.ti = np.concatenate((track.ti, t), axis=0)
        
        return True
      

    def plot_tracker(self, ax1, ax2, ax3, ax4, t='b'):

        for l in self.layers:
            l.draw(ax1, ax2, ax3, ax4, t)
           

             









