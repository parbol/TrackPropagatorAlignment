from src.Module import Module
import numpy as np

class ETL:
    def __init__(self, z, Lx, Ly, ncolumns, nrows, shiftx, shifty, gapx, gapy, innerR, outterR):

        self.firstcenterRightX = shiftx
        self.firstcenterRightY = outterR - shifty
        self.firstcenterLeftX = -shiftx
        self.modules = []        
        for i in range(ncolumns):
            for j in range(nrows):
                x = self.firstcenterRightX + i * (gapx + Lx)
                xleft = self.firstcenterLeftX - i * (gapx + Lx)
                y = self.firstcenterRightY - j * (gapy + Ly)
                xb = x - Lx/2.0
                xbleft = xleft + Lx/2.0
                yb = y - Ly/2.0
                xbm = x + Lx/2.0
                xbmleft = xleft - Lx/2.0
                if np.sqrt(xb**2+yb**2) > innerR and np.sqrt(xbm**2+yb**2) < outterR:
                    m = Module(x, y, z, 0.0, 0.0, 1.0, Lx, Ly)
                    self.modules.append(m)
                if np.sqrt(xbleft**2+yb**2) > innerR and np.sqrt(xbmleft**2+yb**2) < outterR:
                    m = Module(xleft, y, z, 0.0, 0.0, 1.0, Lx, Ly)
                    self.modules.append(m)

    def draw(self, ax1, ax2, ax3, t):

        for m in self.modules:
            m.drawModule(ax1, ax2, ax3, t)

             









