from src.Plane import Plane
from src.Module import Module
from src.BTLRU import BTLRU

import numpy as np
import sys



class BTLTray:
    
    def __init__(self, n, side, x, y, z, euler, TrayWidth, TrayLength, RULength, ModuleLength, ModuleWidth, rphiError, zError, tError):

        self.r = np.asarray([x, y, z])
        sinphi = (self.r[1]/np.sqrt(self.r[0]**2+self.r[1]**2)) 
        cosphi = (self.r[0]/np.sqrt(self.r[0]**2+self.r[1]**2))
        if cosphi > 1.0:
            cosphi = 1.0
        if cosphi < -1.0:
            cosphi = -1.0
        if sinphi >= 0:
            self.phi = np.arccos(cosphi)
        else:
            self.phi = -np.arccos(cosphi) + 2.0*np.pi
        normal = np.asarray([self.r[0], self.r[1], 0.0])
        self.n = normal/np.linalg.norm(normal)
        self.trayNumber = n
        self.side = side
        self.TrayWidth = TrayWidth
        self.TrayLength = TrayLength
        self.RULength = RULength
        self.ModuleLength = ModuleLength
        self.ModuleWidth = ModuleWidth 
        self.rphi_error = rphiError
        self.z_error = zError
        self.t_error = tError
        self.eulerAngles = euler
        self.plane = Plane(self.r[0], self.r[1], self.r[2], self.n[0], self.n[1], self.n[2])
        self.RUSpace = (self.TrayLength - 6.0 * self.RULength)/5.0
        self.zFirstRU = np.abs(z) - TrayLength/2.0 + self.RULength/2.0
        if side == -1:
            self.zFirstRU = -1.0 * self.zFirstRU
        self.RUs = []
        ru11 = BTLRU(1, 1, side, x, y, self.zFirstRU, self.eulerAngles, self.TrayWidth, self.RULength, self.ModuleLength, self.ModuleWidth, self.rphi_error, self.z_error, self.t_error)
        ru12 = BTLRU(1, 2, side, x, y, self.zFirstRU + side * (self.RULength + self.RUSpace), self.eulerAngles, self.TrayWidth, self.RULength, self.ModuleLength, self.ModuleWidth, self.rphi_error, self.z_error, self.t_error)
        ru21 = BTLRU(2, 1, side, x, y, self.zFirstRU + side * 2.0*(self.RULength + self.RUSpace), self.eulerAngles, self.TrayWidth, self.RULength, self.ModuleLength, self.ModuleWidth, self.rphi_error, self.z_error, self.t_error)
        ru22 = BTLRU(2, 2, side, x, y, self.zFirstRU + side * 3.0*(self.RULength + self.RUSpace), self.eulerAngles, self.TrayWidth, self.RULength, self.ModuleLength, self.ModuleWidth, self.rphi_error, self.z_error, self.t_error)
        ru31 = BTLRU(3, 1, side, x, y, self.zFirstRU + side * 4.0*(self.RULength + self.RUSpace), self.eulerAngles, self.TrayWidth, self.RULength, self.ModuleLength, self.ModuleWidth, self.rphi_error, self.z_error, self.t_error)
        ru32 = BTLRU(3, 2, side, x, y, self.zFirstRU + side * 5.0*(self.RULength + self.RUSpace), self.eulerAngles, self.TrayWidth, self.RULength, self.ModuleLength, self.ModuleWidth, self.rphi_error, self.z_error, self.t_error)
        self.RUs.append(ru11)
        self.RUs.append(ru12)
        self.RUs.append(ru21)
        self.RUs.append(ru22)
        self.RUs.append(ru31)
        self.RUs.append(ru32)



    def intersection(self, track):
    
        valid, x, y, z, t = self.plane.intersection(track)
        for m in self.RUs:
            d = np.abs(z - m.r[2])
            if d <= m.RULength/2.0:
                valid, x_, y_, z_, t_ = m.intersection(x, y, z, track)
                if valid:
                    return True, x_, y_, z_, t_
        return False, 0, 0, 0, 0



    def draw(self, ax1, ax2, ax3, t):
        
        for m in self.RUs:
            m.draw(ax1, ax2, ax3, t)

             









