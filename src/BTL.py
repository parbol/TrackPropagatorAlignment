from src.Plane import Plane
from src.Module import Module
from src.BTLGeom import BTLGeom
from src.BTLTray import BTLTray
from src.EulerRotation import EulerRotation
from src.BTLId import BTLId
import numpy as np
import sys

class BTL:
    
    def __init__(self, R, TrayLength, TrayWidth, TrayStartZ, TrayStartPhi, RULength, ModuleLength, ModuleWidth, rphiError, zError, tError, X0):
             
        self.btlNominal = BTLGeom(R, TrayLength, TrayWidth, TrayStartZ, TrayStartPhi, RULength, ModuleLength, ModuleWidth, rphiError, zError, tError, X0)
        self.btlReal = BTLGeom(R, TrayLength, TrayWidth, TrayStartZ, TrayStartPhi, RULength, ModuleLength, ModuleWidth, rphiError, zError, tError, X0)

    
    def fullMeasurement(self, track):
 
        #Intersection
        status, detInfo, v, det = self.btlReal.intersection(track)
        if not status:
            return False

        #Checking with the nominal geometry
        moduleNom = self.btlNominal.getModule(detInfo)
        val, v0, v1, v2, v3 = moduleNom.plane.intersection(track)
        if not val:
            return False
        
        x_ = np.asarray([v0])
        y_ = np.asarray([v1])
        z_ = np.asarray([v2])
        t_ = np.asarray([v3])

        xTracker = track.xi[len(track.xi)-1]
        yTracker = track.yi[len(track.yi)-1]
        zTracker = track.zi[len(track.zi)-1]
        vectorTracker = np.asarray([xTracker, yTracker, zTracker])

        track.xi = np.concatenate((track.xi, x_), axis=0)
        track.yi = np.concatenate((track.yi, y_), axis=0)
        track.zi = np.concatenate((track.zi, z_), axis=0)
        track.ti = np.concatenate((track.ti, t_), axis=0)
        track.det = track.det + det
        globalVector = np.asarray([v0, v1, v2])
        localVector = moduleNom.toLocal(globalVector)
        lx = np.asarray([localVector[0]])
        ly = np.asarray([localVector[1]])
        lz = np.asarray([localVector[2]])
        track.lxi = np.concatenate((track.lxi, lx), axis=0)
        track.lyi = np.concatenate((track.lyi, ly), axis=0)
        track.lzi = np.concatenate((track.lzi, lz), axis=0)


        x = np.asarray([v[0]])
        y = np.asarray([v[1]])
        z = np.asarray([v[2]])
        t = np.asarray([v[3]])
        module = self.btlReal.getModule(detInfo)
        globalVector = np.asarray([v[0], v[1], v[2]])
        localVector = module.toLocal(globalVector)

        detInfoList = [detInfo.side, detInfo.tray, detInfo.RUType, detInfo.RUNumber, detInfo.module]
        track.subdet = track.subdet + [detInfoList]
        localX = localVector[0]
        localY = localVector[1]

        #Spatial uncertainty
        xunc = np.random.normal(0, self.btlReal.rphi_error)
        yunc = np.random.normal(0, self.btlReal.z_error)
        #Multiple scattering
        dl = np.linalg.norm(globalVector-vectorTracker)
        betamomentum = track.betamomentum
        xangle, xdisp = self.btlReal.getScatteringMagnitude(dl, betamomentum)
        yangle, ydisp = self.btlReal.getScatteringMagnitude(dl, betamomentum)
        localX += (xunc + xdisp)
        localY += (yunc + ydisp)
        localVector[0] = localX
        localVector[1] = localY
        newGlobalVector = module.toGlobal(localVector)
        
        #Time uncertainty
        var_t = np.random.normal(0, self.btlReal.t_error)
        newt = v[3] + var_t

        x_meas = np.asarray([newGlobalVector[0]])
        y_meas = np.asarray([newGlobalVector[1]])
        z_meas = np.asarray([newGlobalVector[2]])
        lx_meas = np.asarray([localVector[0]])
        ly_meas = np.asarray([localVector[1]])
        lz_meas = np.asarray([localVector[2]])
        t_meas = np.asarray([newt])
        
        track.xm = np.concatenate((track.xm, x_meas), axis=0)
        track.ym = np.concatenate((track.ym, y_meas), axis=0)
        track.zm = np.concatenate((track.zm, z_meas), axis=0)
        track.tm = np.concatenate((track.tm, t_meas), axis=0)
        track.lxm = np.concatenate((track.lxm, lx_meas), axis=0)
        track.lym = np.concatenate((track.lym, ly_meas), axis=0)
        track.lzm = np.concatenate((track.lzm, lz_meas), axis=0)
        if not module.isInside(localVector):
            return False
        return True
      

   

             









