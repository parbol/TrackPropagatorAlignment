from src.Plane import Plane
from src.Module import Module
import numpy as np
import sys

class BTLModule:
    
    def __init__(self, m, type, n, side, x, y, z, euler, ModuleLength, ModuleWidth, rphiError, zError, tError):

        self.r = np.asarray([x, y, z])
        self.ModuleLength = ModuleLength
        self.ModuleWidth = ModuleWidth 
        self.rphi_error = rphiError
        self.z_error = zError
        self.t_error = tError
        self.moduleNumber = m
        self.type = type
        self.n = n
        self.side = side
        self.eulerAngles = euler
        
        self.module = Module(self.r, self.eulerAngles, self.ModuleWidth, self.ModuleLength)
    
   
    def intersection(x, y, z, track):

        return self.module.intersection(track)


    def draw(self, ax1, ax2, ax3, t):

        self.module.drawModule(ax1, ax2, ax3, t)


        









