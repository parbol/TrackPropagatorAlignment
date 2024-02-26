from src.Plane import Plane
import numpy as np

class Module:

    def __init__(self, x0, y0, z0, nx, ny, nz, Lx, Ly, ):

        self.plane = Plane(x0, y0, z0, nx, ny, nz) 
        self.pLLlocal = np.asarray([-Lx/2.0, -Ly/2.0, 0.0]) 
        self.pLRlocal = np.asarray([Lx/2.0, -Ly/2.0, 0.0]) 
        self.pULlocal = np.asarray([-Lx/2.0, Ly/2.0, 0.0]) 
        self.pURlocal = np.asarray([Lx/2.0, Ly/2.0, 0.0])
        
        self.pLLnominal = self.plane.p + np.asarray(self.plane.rot.dot(self.pLLlocal))[0]
        self.pLRnominal = self.plane.p + np.asarray(self.plane.rot.dot(self.pLRlocal))[0]
        self.pULnominal = self.plane.p + np.asarray(self.plane.rot.dot(self.pULlocal))[0]
        self.pURnominal = self.plane.p + np.asarray(self.plane.rot.dot(self.pURlocal))[0]
   
        self.pLL = self.pLLnominal
        self.pLR = self.pLRnominal
        self.pUL = self.pULnominal
        self.pUR = self.pURnominal
 

    def drawModule(self, ax1, ax2, ax3, t):

        x_start = [self.pLL[0], self.pLR[0], self.pUR[0], self.pUL[0], self.pLL[0]]
        y_start = [self.pLL[1], self.pLR[1], self.pUR[1], self.pUL[1], self.pLL[1]]
        z_start = [self.pLL[2], self.pLR[2], self.pUR[2], self.pUL[2], self.pLL[2]]

        ax1.plot3D(x_start , z_start, y_start, t)
        ax2.plot(x_start, y_start, t)
        ax3.plot(z_start, y_start, t)


    