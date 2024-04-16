from src.Plane import Plane
from src.BTLModule import BTLModule
from src.BTLId import BTLId

import numpy as np
import sys

class BTLRU:
    
    def __init__(self, btl, x, y, z, euler, TrayWidth, RULength, ModuleLength, ModuleWidth):

        self.r = np.asarray([x, y, z])
        normal = np.asarray([self.r[0], self.r[1], 0.0])
        self.n = normal/np.linalg.norm(normal)
        self.TrayWidth = TrayWidth
        self.RULength = RULength
        self.ModuleLength = ModuleLength
        self.ModuleWidth = ModuleWidth 
        self.eulerAngles = euler
        self.btlId = btl
        self.interspaceY = (self.TrayWidth - 3.0 * self.ModuleWidth)/2.0
        self.interspaceZ = (self.RULength - 8.0 * self.ModuleLength)/7.0
        vn = np.asarray([-self.n[1], self.n[0], 0.0])
        self.vn = vn/np.linalg.norm(vn)

        self.zFirstModule = np.abs(z) - self.RULength/2.0 + self.ModuleLength/2.0   
        if self.btlId.side == -1:
            self.zFirstModule = -1.0 * self.zFirstModule
        
        self.Modules = []
        counter = 1
        for i in range(-1, 2):
            for j in range(0, 8):
                xpos = self.r[0] + i * (self.ModuleWidth + self.interspaceY) * self.vn[0]                                
                ypos = self.r[1] + i * (self.ModuleWidth + self.interspaceY) * self.vn[1]
                zpos = self.zFirstModule + self.btlId.side * j * (self.ModuleLength + self.interspaceZ)
                btl = BTLId()
                btl.setTray(self.btlId.tray)
                btl.setSide(self.btlId.side)
                btl.setRU(self.btlId.RUType+1, self.btlId.RUNumber+1)
                btl.setModule(counter)
                mbtl = BTLModule(btl, xpos, ypos, zpos, self.eulerAngles, self.ModuleLength, self.ModuleWidth)
                self.Modules.append(mbtl)
                counter = counter + 1
    
    
    def intersection(self, x, y, z, track):
        for m in self.Modules:
            dz = np.abs(z-m.r[2])
            dxy = np.sqrt((x-m.r[0])**2+(y-m.r[1])**2)
            if dz < m.ModuleLength/2.0 and dxy < m.ModuleWidth/2.0:
                valid, moduleId, point = m.intersection(x, y, z, track)
                if valid:
                    return True, moduleId, point
        return False, [], []
    

    def draw(self, ax1, ax2, ax3, t):

        #self.Modules[0].draw(ax1, ax2, ax3, t)
        #self.Modules[1].draw(ax1, ax2, ax3, t)
        #self.Modules[2].draw(ax1, ax2, ax3, t)

        for m in self.Modules:
            m.draw(ax1, ax2, ax3, t)

             









