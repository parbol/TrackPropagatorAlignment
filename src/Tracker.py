import math
import numpy as np
import matplotlib.pyplot as plt
from random import Random
import sys

class Tracker:
       
    def __init__(self, ri, rz, phi_error, z_error, zsize, centre = [0, 0, 0]):
       
        self.centre = centre
        self.ri = ri
        self.zp = rz
        self.zm = -1.0 * rz
        self.z = np.concatenate((self.zm, self.zp), axis=0)
        self.phi_error = phi_error
        self.z_error = z_error
        self.zsize = zsize      
  
    def plot_tracker(self, ax1, ax2, ax3, fmtb = 'b--', fmte = 'r--'):
        
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

        for zpa in self.z:
            x = []
            yt = []
            yt2 = []
            for i in self.ri:  
                x.append(zpa)
                yt.append(i)
                yt2.append(-i)
                ax3.plot(x, yt, fmte, x, yt2, fmte)   

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
        for tj in tbarrel:
            if tj > tendcap[0]:
                break
            x, y, z = track.eval(tj)
            if z < self.centre[2] + self.zsize/2.0 and z > self.centre[2] - self.zsize/2.0:
                t.append(tj)
        for tj in tendcap:
            x, y, z = track.eval(tj)
            if np.sqrt(x**2+y**2) < self.ri[-1]:
                t.append(tj)

        x, y, z = track.eval(np.array(t))    
        track.lastT = t[-1]
        return x, y, z
    
    
    def plot_intersection(self, track, ax1, ax2, ax3, fmt = 'r*'):
        '''
        Plots the intersection of the tracker and the given track

        Parameters
        ----------
        tracks : list
            list of Track.

        Returns
        -------
        x : list
            3D list with the three cartesian coordinates of the interesection 
            point.

        '''
        # for track in tracks:
        x, y, z = self.intersection(track)
        self.plot_points(x, y, z, ax1, ax2, ax3, fmt)
        return x, y, z    

    def measure(self, track):
        '''
        Simulates the measurement of the intersection point

        Parameters
        ----------
        tracks : list
            list of Track.
        sigma_h: float, optional
            error used for the gaussian. Default value is 0.05

        Returns
        -------
        meast_points_c : list
            3D list with the three cartesian coordinates of the measured 
            interesection point.
        meast_points : float
            theta value for the intersection point, obtained through a gaussian.

        '''
        # for track in tracks:
        x, y, z = self.intersection(track)
        phi = np.arctan2(y, x) + np.random.normal(0, self.phi_error)
        x_meas = self.ri * np.cos(phi)
        y_meas = self.ri * np.sin(phi)
        z_meas = z + np.random.normal(0, self.z_error)
 
        return x_meas, y_meas, z_meas

    def plot_points(self, x, y, z, ax1, ax2, ax3, fmt):

        ax2.plot(x, y, fmt)
        ax1.plot3D(x, z, y, fmt)
        ax3.plot(z, y, fmt)
    
    def plot_meast_points(self, track, ax1, ax2, ax3, fmt = 'bx'):
        '''
        Plots the measured intersection points

        Parameters
        ----------
        tracks : list
            list of Track.

        Returns
        -------
        x : list
            3D list with the three cartesian coordinates of the measured 
            interesection point.
        meast_points : float
            theta value for the intersection point, obtained through a gaussian.

        '''
        # for track in tracks:
        x, y, z = self.measure(track)
        self.plot_points(x, y, z, ax1, ax2, ax3, fmt)
        return x, y, z    

