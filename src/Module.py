from src.Plane import Plane
import numpy as np

class Module:

    def __init__(self, x0, y0, z0, nx, ny, nz, Lx, Ly):

        self.plane = Plane(x0,y0,z0,nx,ny,nz)
        self.pLL = np.asarray([x0 - Lx/2.0, y0 - Ly/2.0, 0.0]) + self.plane.p
        self.pLR = np.asarray([x0 + Lx/2.0, y0 - Ly/2.0, 0.0]) + self.plane.p
        self.pUL = np.asarray([x0 - Lx/2.0, y0 + Ly/2.0, 0.0]) + self.plane.p
        self.pUR = np.asarray([x0 + Lx/2.0, y0 + Ly/2.0, 0.0]) + self.plane.p


    def drawModule(self, ax1, ax2, ax3, t):

        x_start = [self.pLL[0], self.pLR[0], self.pUR[0], self.pUL[0], self.pLL[0]]
        y_start = [self.pLL[1], self.pLR[1], self.pUR[1], self.pUL[1], self.pLL[1]]
        z_start = [self.pLL[2], self.pLR[2], self.pUR[2], self.pUL[2], self.pLL[2]]

        ax1.plot3D(x_start , z_start, y_start, t)
        ax2.plot(x_start, y_start, t)
        ax3.plot(z_start, y_start, t)


    