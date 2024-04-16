from src.Plane import Plane
from src.Module import Module
from src.BTLTray import BTLTray
from src.EulerRotation import EulerRotation
from src.BTLId import BTLId
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
        self.zOfPositiveTrays = self.TrayStartZ + self.TrayLength/2.0
        self.zOfNegativeTrays = -self.zOfPositiveTrays
        self.anglePerTray = 2.0 * np.arcsin((self.TrayWidth/2.0) / self.R)
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
            vx = np.asarray([np.sin(phi), -np.cos(phi), 0.0])
            vy = np.asarray([0.0, 0.0, 1.0])
            vz = np.asarray([np.cos(phi), np.sin(phi), 0.0])
            zp = self.zOfPositiveTrays
            zm = self.zOfNegativeTrays
            euler = EulerRotation()
            euler.setFromVectors(vx, vy, vz)
            eulerAngles.append(euler)
            positiveTrays.append([x, y, zp])
            negativeTrays.append([x, y, zm])
        
        for i in range(18, 36):
            phi = np.pi + self.TrayStartPhi + self.anglePerTray/2.0 + (i-18) * (self.anglePerTray+self.traySpace)
            x = self.R * np.cos(phi)
            y = self.R * np.sin(phi)
            vx = np.asarray([np.sin(phi), -np.cos(phi), 0.0])
            vy = np.asarray([0.0, 0.0, 1.0])
            vz = np.asarray([np.cos(phi), np.sin(phi), 0.0])
            euler = EulerRotation()
            euler.setFromVectors(vx, vy, vz)
            eulerAngles.append(euler)
            zp = self.zOfPositiveTrays
            zm = self.zOfNegativeTrays
            positiveTrays.append([x, y, zp])
            negativeTrays.append([x, y, zm])

        self.pTrays = []
        self.mTrays = []
        for i, t in enumerate(positiveTrays):
            btl = BTLId()
            btl.setTray(i+1)
            btl.setSide(1)
            tray = BTLTray(btl, t[0], t[1], t[2], eulerAngles[i], self.TrayWidth, self.TrayLength, self.RULength, self.ModuleLength, self.ModuleWidth)
            self.pTrays.append(tray)
        for i, t in enumerate(negativeTrays):
            btl = BTLId()
            btl.setTray(i+1)
            btl.setSide(-1)
            tray = BTLTray(btl, t[0], t[1], t[2], eulerAngles[i], self.TrayWidth, self.TrayLength, self.RULength, self.ModuleLength, self.ModuleWidth)
            self.mTrays.append(tray)
                
   
    def getClosestTray(self, phi, trays):
        
        index = -1
        minPhi = 9999999.0
        for j, i in enumerate(trays):
            if np.abs(i.phi-phi) < minPhi:
                minPhi = np.abs(i.phi-phi)
                index = j
        return index


    def intersection(self, track):

        phi = track.phi
        eta = track.eta
        trays = []
        side = 0
        if eta >= 0:
            trays = self.pTrays
            side = 1
        else:
            trays = self.mTrays
            side = -1

        trayNumber = self.getClosestTray(phi, trays)
        navigationWindow = 2
        navigationList = []
        for i in range(0, navigationWindow+1):
            navigationList.append(i)
            if i != 0:
                navigationList.append(-i)
        for i in navigationList:
            checkTrayNumber = (trayNumber + i) % 36
            valid, id, point = trays[checkTrayNumber].intersection(track)
            if valid:
                return True, id, point, [3]
        return False, [], [], [3]


    def getModule(self, detInfo):

        tray = detInfo.tray
        side = detInfo.side
        RU = detInfo.RU
        module = detInfo.module
       
        if side > 0:

            return self.pTrays[tray].RUs[RU].Modules[module].module
        
        else:

            return self.mTrays[tray].RUs[RU].Modules[module].module
        

    def fullMeasurement(self, track):
 
        #Intersection
        status, detInfo, v, det = self.intersection(track)
        if not status:
            return False
        x = np.asarray([v[0]])
        y = np.asarray([v[1]])
        z = np.asarray([v[2]])
        t = np.asarray([v[3]])

        track.xi = np.concatenate((track.xi, x), axis=0)
        track.yi = np.concatenate((track.yi, y), axis=0)
        track.zi = np.concatenate((track.zi, z), axis=0)
        track.ti = np.concatenate((track.ti, t), axis=0)
        track.det = track.det + det
        
        module = self.getModule(detInfo)
        
        globalVector = np.asarray([v[0], v[1], v[2]])
        localVector = module.toLocal(globalVector)
        localX = localVector[0]
        localY = localVector[1]
        xunc = np.random.normal(0, self.rphi_error)
        yunc = np.random.normal(0, self.z_error)
        localX += xunc
        localY += yunc
        localVector[0] = localX
        localVector[1] = localY
        newGlobalVector = module.toGlobal(localVector)
        var_t = np.random.normal(0, self.t_error)
        newt = v[3] + var_t
        x_meas = np.asarray([newGlobalVector[0]])
        y_meas = np.asarray([newGlobalVector[1]])
        z_meas = np.asarray([newGlobalVector[2]])
        t_meas = np.asarray([newt])
        
        track.xm = np.concatenate((track.xm, x_meas), axis=0)
        track.ym = np.concatenate((track.ym, y_meas), axis=0)
        track.zm = np.concatenate((track.zm, z_meas), axis=0)
        track.tm = np.concatenate((track.tm, t_meas), axis=0)
       
        return True
      

    def draw(self, ax1, ax2, ax3, t):

        for m in self.pTrays:
            m.draw(ax1, ax2, ax3, t)
        for m in self.mTrays:
            m.draw(ax1, ax2, ax3, t)


             









