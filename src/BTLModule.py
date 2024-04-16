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
    
   
    def intersection(self, x, y, z, track):
        valid, x, y, z, t = self.module.intersection(track)
        return valid, self.btlId, [x, y, z, t]


    def draw(self, ax1, ax2, ax3, t):

        self.module.drawModule(ax1, ax2, ax3, t)


        









