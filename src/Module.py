from src.Plane import Plane
import numpy as np
import sys

class Module:

    def __init__(self, xnom, anglenom, x, angle, Lx, Ly):

        #In local coordinates the module Z is always parallel to global Z
        #RotX
        self.Lx = Lx
        self.Ly = Ly

        self.xnom = xnom
        matxnom_ = [[1.0, 0, 0],
                [0.0, np.cos(anglenom[0]), -np.sin(anglenom[0])],
                [0.0, np.sin(anglenom[0]), np.cos(anglenom[0])]]
        matxnom = np.asmatrix(matxnom_)
        #RotY
        matynom_ = [[np.cos(anglenom[1]), 0.0, -np.sin(anglenom[1])],
                [0.0, 1.0, 0.0],
                [np.sin(anglenom[1]), 0.0, np.cos(anglenom[1])]]
        matynom = np.asmatrix(matynom_)
        #RotZ
        matznom_ = [[np.cos(anglenom[2]), -np.sin(anglenom[2]), 0.0],
                [np.sin(anglenom[2]), np.cos(anglenom[2]), 0.0],
                [0.0, 0.0, 1.0]]
        matznom = np.asmatrix(matznom_)

        self.rotnom = matxnom.dot(matynom.dot(matznom))
        self.invrotnom = np.linalg.inv(self.rotnom)
        
        znom = np.asarray([0.0, 0.0, 1.0])
        nnom = np.asarray(self.rotnom.dot(znom))[0]   

        #In local coordinates the module Z is always parallel to global Z
        #RotX
        self.x = x
        matx_ = [[1.0, 0, 0],
                [0.0, np.cos(angle[0]), -np.sin(angle[0])],
                [0.0, np.sin(angle[0]), np.cos(angle[0])]]
        matx = np.asmatrix(matx_)
        #RotY
        maty_ = [[np.cos(angle[1]), 0.0, -np.sin(angle[1])],
                [0.0, 1.0, 0.0],
                [np.sin(angle[1]), 0.0, np.cos(angle[1])]]
        maty = np.asmatrix(maty_)
        #RotZ
        matz_ = [[np.cos(angle[2]), -np.sin(angle[2]), 0.0],
                [np.sin(angle[2]), np.cos(angle[2]), 0.0],
                [0.0, 0.0, 1.0]]
        matz = np.asmatrix(matz_)

        self.rot = matx.dot(maty.dot(matz))
        self.invrot = np.linalg.inv(self.rot)

        z = np.asarray([0.0, 0.0, 1.0])
        n = np.asarray(self.rot.dot(z))[0]  

        self.pLLlocal = np.asarray([-Lx/2.0, -Ly/2.0, 0.0]) 
        self.pLRlocal = np.asarray([Lx/2.0, -Ly/2.0, 0.0]) 
        self.pULlocal = np.asarray([-Lx/2.0, Ly/2.0, 0.0]) 
        self.pURlocal = np.asarray([Lx/2.0, Ly/2.0, 0.0])
        
        self.planeNom = Plane(xnom[0], xnom[1], xnom[2], nnom[0], nnom[1], nnom[2]) 
        self.pLLnominal = self.toGlobalNom(self.pLLlocal)
        self.pLRnominal = self.toGlobalNom(self.pLRlocal)
        self.pULnominal = self.toGlobalNom(self.pULlocal)
        self.pURnominal = self.toGlobalNom(self.pURlocal)

        self.plane = Plane(x[0], x[1], x[2], n[0], n[1], n[2]) 
        self.pLL = self.toGlobal(self.pLLlocal)
        self.pLR = self.toGlobal(self.pLRlocal)
        self.pUL = self.toGlobal(self.pULlocal)
        self.pUR = self.toGlobal(self.pURlocal)
 
    def toGlobal(self, v):

        return self.x + np.asarray(self.rot.dot(v))[0]
    
    def toGlobalNom(self, v):

        return self.xnom + np.asarray(self.rotnom.dot(v))[0]
    
    def toLocal(self, v):

        return np.asarray(self.invrot.dot(v - self.x))[0]
    
    def toLocalNom(self, v):

        return np.asarray(self.invrotnom.dot(v - self.xnom))[0]
    
    def isInside(self, p):

        if p[0] < -self.Lx/2.0 or p[0] > self.Lx/2.0:
            return False
        if p[1] < -self.Ly/2.0 or p[1] > self.Ly/2.0:
            return False
        return True
    

    def intersection(self, track):
        
        x, y, z, t = self.plane.intersection(track)
        p = np.asarray([x, y, z])
        plocal = self.toLocal(p)
        if self.isInside(plocal):
            return True, x, y, z, t
        return False, x, y, z, t
    
    def intersectionNom(self, track):
        
        x, y, z, t = self.planeNom.intersection(track)
        p = np.asarray([x, y, z])
        plocal = self.toLocal(p)
        if self.isInside(plocal):
            return True, x, y, z, t
        return False, x, y, z, t

    def drawModule(self, ax1, ax2, ax3, t):

        x_start = [self.pLL[0], self.pLR[0], self.pUR[0], self.pUL[0], self.pLL[0]]
        y_start = [self.pLL[1], self.pLR[1], self.pUR[1], self.pUL[1], self.pLL[1]]
        z_start = [self.pLL[2], self.pLR[2], self.pUR[2], self.pUL[2], self.pLL[2]]

        ax1.plot3D(x_start , z_start, y_start, t)
        ax2.plot(x_start, y_start, t)
        ax3.plot(z_start, y_start, t)


    