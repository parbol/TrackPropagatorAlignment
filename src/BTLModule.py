from src.Plane import Plane
from src.Module import Module
from src.BTLId import BTLId

import numpy as np
import sys

class BTLModule:
    
    def __init__(self, btl, x, y, z, euler, ModuleLength, ModuleWidth):

        self.r = np.asarray([x, y, z])
        self.ModuleLength = ModuleLength
        self.ModuleWidth = ModuleWidth 
        self.btlId = btl
        self.eulerAngles = euler
        self.module = Module(self.r, self.eulerAngles, self.ModuleWidth, self.ModuleLength)
    
   
    def write(self, f):

        cad = '{side} {tray} {RUType} {RUNumber} {module}'.format(side = str(self.btlId.side),
                                                                  tray = str(self.btlId.tray),
                                                                  RUType = str(self.btlId.RUType),
                                                                  RUNumber = str(self.btlId.RUNumber),
                                                                  module = str(self.btlId.module))
        f.write(cad + '\n')  

    def intersection(self, x, y, z, track):
        valid, x, y, z, t = self.module.intersection(track)
        return valid, self.btlId, [x, y, z, t] 


    def draw(self, ax1, ax2, ax3, t):

        self.module.drawModule(ax1, ax2, ax3, t)


        









