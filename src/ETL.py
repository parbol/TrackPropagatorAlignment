from src.Module import Module
import numpy as np

class ETL:
    def __init__(self, z, Lx, Ly, ncolumns, nrows, shiftx, shifty, gapx, gapy, innerR, outterR):

        self.firstcenterRightX = -shiftx
        self.firstcenterRightY = outterR - shifty
        self.firstcenterLeftX = +shiftx
        self.modules = []        
        for i in range(ncolumns):
            for j in range(nrows):
                x = self.firstcenterRightX - i * (gapx + Lx)
                xleft = self.firstcenterLeftX + i * (gapx + Lx)
                y = self.firstcenterRightY - j * (gapy + Ly)
                xb = x + Lx/2.0
                xbleft = xleft - Lx/2.0
                yb = y - Ly/2.0
                xbm = x - Lx/2.0
                xbmleft = xleft + Lx/2.0
                if np.sqrt(xb**2+yb**2) > innerR and np.sqrt(xbm**2+yb**2) < outterR:
                    r = np.sqrt(y**2+x**2)
                    zt = z + (r-30.0) / np.tan(np.pi/6.0)
                    vnom = [x, y, z]
                    anglenom = [0.0, 0.0, 0.0]
                    v = [x, y, z]
                    angle = [0.0, 0.0, 0.0]
                    m = Module(vnom, anglenom, v, angle, Lx, Ly)
                    self.modules.append(m)
                if np.sqrt(xbleft**2+yb**2) > innerR and np.sqrt(xbmleft**2+yb**2) < outterR:
                    r = np.sqrt(y**2+xleft**2)
                    zt = z + (r-30.0) / np.tan(np.pi/6.0)
                    vnom = [xleft, y, z]
                    anglenom = [0.0, 0.0, 0.0]
                    v = [xleft, y, z]
                    angle = [0.0, 0.0, 0.0]
                    m = Module(vnom, anglenom, v, angle, Lx, Ly)
                    self.modules.append(m)

    def draw(self, ax1, ax2, ax3, t):

        for m in self.modules:
            m.drawModule(ax1, ax2, ax3, t)

             









