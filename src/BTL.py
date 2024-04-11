from src.Plane import Plane
from src.Module import Module
from src.BTLTray import BTLTray
from src.EulerRotation import EulerRotation
import numpy as np
import sys

class BTL:
    
    def __init__(self, R, TrayLength, TrayWidth, TrayStartZ, TrayStartPhi, RULength, ModuleLength, ModuleWidth, rphiError, zError, tError):

        ####################################################################################################
        #                                   Representation of a BTL disk                                   #
        # R: Radius at which the trays are located
        # TrayLength: Length in Z of the tray
        # TrayWidth: Width (rphi) of the tray
        # TrayStartZ: Start of the tray with respect to 0 in Z
        # TrayStartPhi: Start of the tray with respect to 0 in phi
        # RuLength: Readunit length in Z
        # ModuleLength: Module length in Z
        # ModuleWidth: Module width in rphi
        # rphi_error: uncertainty in rphi
        # z_error: undertainty in z
        # t_error: time uncertainty
        #####################################################################################################
        self.R = R
        self.TrayLength = TrayLength
        self.TrayWidth = TrayWidth
        self.TrayStartZ = TrayStartZ
        self.TrayStartPhi = TrayStartPhi
        self.RULength = RULength
        self.ModuleLength = ModuleLength
        self.ModuleWidth = ModuleWidth
        self.rphi_error = rphiError
        self.z_error = zError
        self.t_error = tError
        self.zOfPositiveTrays = self.TrayStartZ + self.Lz/2.0
        self.zOfNegativeTrays = -self.zOfPositiveTrays
        self.anglePerTray = 2.0 * np.arctan((self.TrayWidth/2.0) / self.R)
        self.traySpace = ((np.pi - 2.0 * self.TrayStartPhi) - 18.0 * self.anglePerTray)/17.0
        if self.traySpace < 0: 
            print('The geometry is not valid')
            sys.exit()
        
        positiveTrays = []
        negativeTrays = []
        eulerAngles = []
        for i in range(0, 18):
            phi = self.TrayStartPhi + self.anglePerTray/2.0 + i * (self.anglePerTray+self.traySpace)
            x = self.R * np.cos(phi)
            y = self.R * np.sin(phi)
            zp = self.zOfPositiveTrays
            zm = self.zOfNegativeTrays
            euler = EulerRotation(vx = vx, vy = vy, vz = vz)
            eulerAngles.append(euler)
            positiveTrays.append([x, y, zp])
            negativeTrays.append([x, y, zm])
        
        for i in range(19, 36):
            phi = 2.0 * self.TrayStartPhi + self.TrayStartPhi + self.anglePerTray/2.0 + i * (self.anglePerTray+self.traySpace)
            x = self.R * np.cos(phi)
            y = self.R * np.sin(phi)
            vx = np.asarray([np.sin(phi), -np.cos(phi), 0.0])
            vy = np.asarray([0.0, 0.0, 1.0])
            vz = np.asarray([np.cos(phi), np.sin(phi), 0.0])
            euler = EulerRotation(vx = vx, vy = vy, vz = vz)
            eulerAngles.append(euler)
            zp = self.zOfPositiveTrays
            zm = self.zOfNegativeTrays
            positiveTrays.append([x, y, zp])
            negativeTrays.append([x, y, zm])

        self.pTrays = []
        self.mTrays = []
        for i, t in enumerate(positiveTrays):
            tray = BTLTray(i+1, 1, t[0], t[1], t[2], eulerAngles[i], self.TrayWidth, self.TrayLength, self.RULength, self.ModuleLength, self.ModuleWidth, self.rphi_error, self.z_error, self.t_error)
            self.pTrays.append(tray)
        for i, t in enumerate(negativeTrays):
            tray = BTLTray(i+1, -1, t[0], t[1], t[2], eulerAngles[i], self.TrayWidth, self.TrayLength, self.RULength, self.ModuleLength, self.ModuleWidth, self.rphi_error, self.z_error, self.t_error)
            self.mTrays.append(tray)
                

      

    def intersection(self, track):

        valid, x_, y_, z_, t_ = self.plane.intersection(track)
        if not valid:
            print('Plane intersetion not valid')
            return False, False, [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [-1]
        p = np.asarray([x_, y_, z_])
        status = statusm = True
        x = y = z = t = xnom = ynom = znom = tnom = 0.0
        v = [x, y, z, t]
        vn = [xnom, ynom, znom, tnom]
        for i, m in enumerate(self.modules):
            d = p - m.xnom
            #Is the module close enough?
            if d[0]**2 + d[1]**2 < 100.0:
                status, x, y, z, t = m.intersection(track)      
                statusm, xnom, ynom, znom, tnom = m.intersectionNom(track)
                if status == True:
                    x = np.asarray([x, y, z])
                    xt = m.toGlobalNom(m.toLocal(x))
                    v = [xt[0], xt[1], xt[2], t]
                    vn = [xnom, ynom, znom, tnom]
                    return True, statusm, v, vn, [2]
        return False, statusm, v, vn, [2]


    def fullMeasurement(self, track):
 
        #Intersection
        status, statusm, v, vn, det = self.intersection(track)
        if not status:
            return False
        
        x = np.asarray([v[0]])
        y = np.asarray([v[1]])
        z = np.asarray([v[2]])
        t = np.asarray([v[3]])
        xn = np.asarray([vn[0]])
        yn = np.asarray([vn[1]])
        zn = np.asarray([vn[2]])
        tn = np.asarray([vn[3]])

        track.xi = np.concatenate((track.xi, x), axis=0)
        track.yi = np.concatenate((track.yi, y), axis=0)
        track.zi = np.concatenate((track.zi, z), axis=0)
        track.ti = np.concatenate((track.ti, t), axis=0)
        track.xin = np.concatenate((track.xin, xn), axis=0)
        track.yin = np.concatenate((track.yin, yn), axis=0)
        track.zin = np.concatenate((track.zin, zn), axis=0)
        track.tin = np.concatenate((track.tin, tn), axis=0)
        track.det = track.det + det
        
        #Measurement     
        phi = np.arctan2(y, x)
        r = np.sqrt(y**2 + x**2)
        phin = np.arctan2(yn, xn)
        rn = np.sqrt(yn**2 + xn**2)
        
        sigma_rphi = self.rphi_error / r
        sigma_z = np.zeros(x.shape) + self.z_error
        var_sigma = np.random.normal(0, sigma_rphi)
        phi = phi + var_sigma
        phin = phin + var_sigma
        var_z = np.random.normal(0, sigma_z)
        z = track.zi + var_z
        zn = track.zin + var_z
        x_meas = r * np.cos(phi)
        y_meas = r * np.sin(phi)
        z_meas = z 
        x_measn = rn * np.cos(phin)
        y_measn = rn * np.sin(phin)
        z_measn = zn 
        sigma_t = np.zeros(x.shape) + self.t_error
        var_t = np.random.normal(0, sigma_t)
        t_meas = t + var_t
        t_measn = tn + var_t
        
        track.xm = np.concatenate((track.xm, x_meas), axis=0)
        track.ym = np.concatenate((track.ym, y_meas), axis=0)
        track.zm = np.concatenate((track.zm, z_meas), axis=0)
        track.tm = np.concatenate((track.tm, t_meas), axis=0)
        track.xmn = np.concatenate((track.xmn, x_measn), axis=0)
        track.ymn = np.concatenate((track.ymn, y_measn), axis=0)
        track.zmn = np.concatenate((track.zmn, z_measn), axis=0)
        track.tmn = np.concatenate((track.tmn, t_measn), axis=0)

        return True
    
    
    def XYfromNormalVector(self, n):
    
        siny = -n[0]    
        #Pathological case
        if siny*siny > 1.0:
            if siny > 0.0:
                return np.asarray([0.0, np.pi/2.0, 0.0])
            else:
                return np.asarray([0.0, -np.pi/2.0, 0.0])
        cosy = np.sqrt(1.0 - siny*siny)    
        cosx = n[2]/cosy
        #We set the range of cosx
        if np.abs(cosx) > 1.0:
            cosx = 1.0

        y = np.arccos(cosy)  
        x = np.arccos(cosx)
    
        y = np.arcsin(np.sin(y))
        x = np.arcsin(np.sin(x))


        sinxp = np.abs(np.sin(x))
        sinxm = -sinxp
        cosxp = np.abs(np.cos(x))
        cosxm = -cosxp
        cosyp = np.abs(np.cos(y))
        cosym = -cosyp
    
        if siny >= 0:          
            if self.normalVectorMatches(n, sinxp, siny, cosxp, cosyp):
                x = x
            elif self.normalVectorMatches(n, sinxp, siny, cosxp, cosym):
                y = np.pi-y
            elif self.normalVectorMatches(n, sinxp, siny, cosxm, cosyp):
                x = np.pi-x
            elif self.normalVectorMatches(n, sinxp, siny, cosxm, cosym):
                x = np.pi-x
                y = np.pi-y
            elif self.normalVectorMatches(n, sinxm, siny, cosxp, cosyp):
                x = -x
            elif self.normalVectorMatches(n, sinxm, siny, cosxp, cosym):
                x = -x
                y = np.pi-y
            elif self.normalVectorMatches(n, sinxm, siny, cosxm, cosyp):
                x = np.pi + x
            else:
                x = np.pi + x
                y = np.pi-y
        else:
            if self.normalVectorMatches(n, sinxp, siny, cosxp, cosyp):
                y = -y
            elif self.normalVectorMatches(n, sinxp, siny, cosxp, cosym):
                y = np.pi+y
            elif self.normalVectorMatches(n, sinxp, siny, cosxm, cosyp):
                x = np.pi-x
                y = -y
            elif self.normalVectorMatches(n, sinxp, siny, cosxm, cosym):
                x = np.pi-x
                y = np.pi+y
            elif self.normalVectorMatches(n, sinxm, siny, cosxp, cosyp):
                x = -x
                y = -y
            elif self.normalVectorMatches(n, sinxm, siny, cosxp, cosym):
                x = -x
                y = np.pi+y
            elif self.normalVectorMatches(n, sinxm, siny, cosxm, cosyp):
                x = np.pi + x
                y = -y
            else:
                x = np.pi + x
                y = np.pi+y
        return np.asarray([x,y,0.0])


    def normalVectorMatches(self, n, sinxm, siny, cosxp, cosyp):
        
        tol = 0.001
        A = self.makeXYspecialMatrix(sinxm, siny, cosxp, cosyp)
        z = np.asarray([0.0, 0.0, 1.0])
        vz = np.asarray(A.dot(z))[0]
        if np.abs(vz[0]-n[0]) < tol and np.abs(vz[1]-n[1]) < tol and np.abs(vz[2]-n[2]) < tol:   
            return True
        return False


    def makeXYspecialMatrix(self, sinx, siny, cosx, cosy):

        A = np.asmatrix([[cosy, 0, -siny],
                        [-sinx*siny, cosx, -sinx*cosy],
                        [cosx*siny, sinx, cosx*cosy]])
        return A


    def draw(self, ax1, ax2, ax3, t):

        for m in self.modules:
            m.drawModule(ax1, ax2, ax3, t)

             









