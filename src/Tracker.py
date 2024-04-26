import math
import numpy as np
import matplotlib.pyplot as plt
from random import Random
import sys

class Tracker:
      
    def __init__(self, ri, rz, phi_error, z_error, zsize, centre = [0, 0, 0]):
        ###################################################################################
        #                            A simple tracker model                               #
        # This is a very simplified model of the tracker with the following paramenters:  #
        # ri: a list with the radius of the various layers of the tracker barrel          #
        # rz: a list with the z position of the various layers of the tracker endcap      #
        # rphi_error: resolution uncertainty in the rphi variable                         #
        # z_error: resolution uncertainty in the z variable                               #
        # z_size: determines the end of the barrel in z coordinates                       #
        # centre: centre of the tracker (0,0,0) by default                                #
        ###################################################################################
        self.centre = centre
        self.ri = ri
        self.zp = rz
        self.zm = -1.0 * rz
        self.z = np.concatenate((self.zm, self.zp), axis=0)
        self.rphi_error = phi_error
        self.z_error = z_error
        self.zsize = zsize      
  

    def plot_tracker(self, ax1, ax2, ax3, ax4, fmtb = 'b--', fmte = 'r--'):
        
        theta = np.linspace(0, 2 * np.pi, 20)
        zpace = np.linspace(-self.zsize / 2.0, self.zsize / 2.0, 10)
        
        for zpa in zpace:
            for i in self.ri:
                x = []
                y = []
                z = []
                for th in theta:
                    x.append(i * np.cos(th))
                    y.append(i * np.sin(th))
                    z.append(zpa)
                ax1.plot3D(x, z, y, fmtb, alpha=0.2)
        
        for zpa in self.z:
            for i in self.ri:
                x = []
                y = []
                z = []
                for th in theta:
                    x.append(i * np.cos(th))
                    y.append(i * np.sin(th))
                    z.append(zpa)
                ax1.plot3D(x, z, y, fmte, alpha=0.2)
                    
        for i in self.ri:
            xt = []
            yt = []
            for th in theta:
                xt.append(i * np.cos(th))
                yt.append(i * np.sin(th))
            ax2.plot(xt, yt, fmtb)
        
        for i in self.ri:
            zt = []
            yt = []
            yt2 = []
            for zpa in zpace:
                zt.append(zpa)
                yt.append(i)
                yt2.append(-i) 
            ax3.plot(zt, yt, fmtb, zt, yt2, fmtb)

        for i in self.ri:
            zt = []
            xt = []
            xt2 = []
            for zpa in zpace:
                zt.append(zpa)
                xt.append(i)
                xt2.append(-i) 
            ax4.plot(zt, xt, fmtb, zt, xt2, fmtb)

        for zpa in self.z:
            x = []
            yt = []
            yt2 = []
            for i in self.ri:  
                x.append(zpa)
                yt.append(i)
                yt2.append(-i)
                ax3.plot(x, yt, fmte, x, yt2, fmte)   

        for zpa in self.z:
            x = []
            yt = []
            yt2 = []
            for i in self.ri:  
                x.append(zpa)
                yt.append(i)
                yt2.append(-i)
                ax4.plot(x, yt, fmte, x, yt2, fmte)   

    def convertAngle(self, theta):
        if theta > np.pi:
            return theta - np.pi * 2.0
        elif theta < -np.pi:
            return theta + np.pi * 2.0
        return theta 
 
  
    def intersection(self, track):
                  
        ri2 = self.ri**2
        x_0 = track.x_c
        y_0 = track.y_c
        z_0 = track.dz
               
        delta = -track.rt**2 + x_0**2 + y_0**2 + ri2
        
        ax = 4. * x_0**2 + 4. * y_0**2
        bx = -4. * delta * x_0
        cx = delta**2 - 4. * y_0**2 * ri2
        ay = 4. * x_0**2 + 4. * y_0**2
        by = -4. * delta * y_0
        cy = delta**2 - 4. * x_0**2 * ri2
      
        rootx = bx**2 - 4.0 * ax * cx
        rooty = by**2 - 4.0 * ay * cy
        
        x1_c = np.zeros(rootx.shape)
        x2_c = np.zeros(rootx.shape)
        y1_c = np.zeros(rootx.shape)
        y2_c = np.zeros(rootx.shape)
        tbarrel = []
        for j, troot in enumerate(rootx):
            if rootx[j] < 0:
                if rooty[j] < 0:
                    #raise ValueError('Problem with one intersection')
                    continue
                else:       
                    y1_c[j] = (-by[j] - np.sqrt(rooty[j])) / (2.0 * ay)
                    y2_c[j] = (-by[j] + np.sqrt(rooty[j])) / (2.0 * ay)
                    if abs((y1_c[j] - track.y_c) / track.rt) > 1.0:
                        #raise ValueError('Problem with one intersection')
                        continue
                    if abs((y2_c[j] - track.y_c) / track.rt) > 1.0:
                        #raise ValueError('Problem with one intersection')
                        continue
                    t1 = 1.0 / track.w * self.convertAngle((track.phi + np.arccos((y1_c[j] - track.y_c) / track.rt)))
                    t2 = 1.0 / track.w * self.convertAngle((track.phi + np.arccos((y2_c[j] - track.y_c) / track.rt)))
                    t3 = 1.0 / track.w * self.convertAngle((track.phi - np.arccos((y1_c[j] - track.y_c) / track.rt)))
                    t4 = 1.0 / track.w * self.convertAngle((track.phi - np.arccos((y2_c[j] - track.y_c) / track.rt)))
                    tarr = np.sort(np.array([t1, t2, t3, t4]), axis=None)
                    tv = 99999
                    for i in tarr:
                        if i > 0:
                            tv = i
                            break
                    tbarrel.append(tv)               
            else:
                x1_c[j] = (-bx[j] - np.sqrt(rootx[j])) / (2.0 * ax)
                x2_c[j] = (-bx[j] + np.sqrt(rootx[j])) / (2.0 * ax)
                if abs((x1_c[j] - track.x_c) / track.rt) > 1.0:
                    #raise ValueError('Problem with one intersection')
                    continue
                if abs((x2_c[j] - track.x_c) / track.rt) > 1.0:
                    #raise ValueError('Problem with one intersection')
                    continue
                t1 = 1.0 / track.w * self.convertAngle((track.phi + np.arcsin((x1_c[j] - track.x_c) / track.rt)))
                t2 = 1.0 / track.w * self.convertAngle((track.phi + np.arcsin((x2_c[j] - track.x_c) / track.rt)))
                t3 = 1.0 / track.w * self.convertAngle((track.phi + np.pi - np.arcsin((x1_c[j] - track.x_c) / track.rt)))
                t4 = 1.0 / track.w * self.convertAngle((track.phi + np.pi - np.arcsin((x2_c[j] - track.x_c) / track.rt)))
                tarr = np.sort(np.array([t1, t2, t3, t4]), axis=None)
                tv = 99999
                for i in tarr:
                    if i > 0:
                        tv = i
                        break
                tbarrel.append(tv)
        
        if track.pz == 0:
            tendcap = []
        else:
            tendcapp = []
            tendcapm = []
            for i in self.zp:
                tval = (track.gamma*track.m)/(30.0*track.pz)*(i-track.dz)
                tendcapp.append(tval)
            tendcap = tendcapp
            if tendcapp[0] < 0:
                for i in self.zm:
                    tval = (track.gamma*track.m)/(30.0*track.pz)*(i-track.dz)
                    tendcapm.append(tval)
                tendcap = tendcapm

        t = []
        det = []
        for tj in tbarrel:
            if len(tendcap) != 0 and tj > tendcap[0]:
                break
            x, y, z = track.eval(tj)
            if z < self.centre[2] + self.zsize/2.0 and z > self.centre[2] - self.zsize/2.0:
                t.append(tj)
                det.append(0)
        for tj in tendcap:
            x, y, z = track.eval(tj)
            if np.sqrt(x**2+y**2) < self.ri[-1]:
                t.append(tj)
                det.append(1)

        x, y, z = track.eval(np.array(t))    
        #track.lastT = t[-1]
        return x, y, z, t, det
    
    
    def makeIntersection(self, track):
       
        # for track in tracks:
        x, y, z, t, det = self.intersection(track)
        copyxi = np.copy(track.xi)
        copyyi = np.copy(track.yi)
        copyzi = np.copy(track.zi)
        copyti = np.copy(track.ti)
        track.xi = np.concatenate((track.xi, x), axis=0)
        track.yi = np.concatenate((track.yi, y), axis=0)
        track.zi = np.concatenate((track.zi, z), axis=0)
        track.ti = np.concatenate((track.ti, t), axis=0)
        track.xin = np.concatenate((copyxi, x), axis=0)
        track.yin = np.concatenate((copyyi, y), axis=0)
        track.zin = np.concatenate((copyzi, z), axis=0)
        track.tin = np.concatenate((copyti, t), axis=0)
        track.det = det


    def measure(self, track):
        
        phi = np.arctan2(track.yi, track.xi)
        r = np.sqrt(track.yi**2 + track.xi**2)
        sigma_rphi = self.rphi_error / r
        sigma_z = np.zeros(track.xi.shape) + self.z_error
        phi = phi + np.random.normal(0, sigma_rphi)
        z = track.zi + np.random.normal(0, sigma_z)
        x_meas = r * np.cos(phi)
        y_meas = r * np.sin(phi)
        z_meas = z 

        return x_meas, y_meas, z_meas    
    
    
    def makeMeasurement(self, track):
       
        # for track in tracks:
        x, y, z = self.measure(track)
        copyxm = np.copy(track.xm)
        copyym = np.copy(track.ym)
        copyzm = np.copy(track.zm)
        copytm = np.copy(track.tm)

        track.xm = np.concatenate((track.xm, x), axis=0)
        track.ym = np.concatenate((track.ym, y), axis=0)
        track.zm = np.concatenate((track.zm, z), axis=0)
        track.tm = np.concatenate((track.tm, track.ti), axis=0)
        track.xmn = np.concatenate((copyxm, x), axis=0)
        track.ymn = np.concatenate((copyym, y), axis=0)
        track.zmn = np.concatenate((copyzm, z), axis=0)
        track.tmn = np.concatenate((copytm, track.ti), axis=0)

    def fullMeasurement(self, track):

        self.makeIntersection(track)
        self.makeMeasurement(track)



        
  