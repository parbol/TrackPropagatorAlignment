# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 09:55:11 2022

@author: franm
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from random import Random
import sys

class Tracker:
    '''
    Tracker class, which creates a tracker, and is able to do and study some
    interactions with the tracks, as well as plotting themselves and some
    other things.
    '''
    
    def __init__(self, ri, phi_error, z_error, zsize, centre = [0, 0, 0]):
        '''
        Creates the track with the required parameters

        Parameters
        ----------
        ri : list
            list with the radius of the different layers the tracker has.
        centre : tuple, optional
            coordinates of the centre where the tracker is located. The default
            is [0, 0].

        Returns
        -------
        None.

        '''
        self.centre = centre
        self.ri = ri
        self.phi_error = phi_error
        self.z_error = z_error
        self.zsize = zsize
        
  
    def plot_tracker(self, fmt = 'b--'):
        '''
        Plots the tracker
        
        Parameters
        ----------
        fmt : string
            colour and style desired for the tracker representation. Default is
            'b--', which is a blue dashed line
        Returns
        -------
        None.

        '''
     
        # fig = plt.figure(figsize = plt.figaspect(0.3))
 
        # ax1 = fig.add_subplot(1, 3, 2, projection = '3d')
        theta = np.linspace(0, 2 * np.pi, 20)
        zpace = np.linspace(-self.zsize / 2.0, self.zsize / 2.0, 10)
        count = 0
        for zpa in zpace:
            # print(zpa)
            # if count % 2 == 0:
            for i in self.ri:
                x = []
                y = []
                z = []
                for th in theta:
                    x.append(i * np.cos(th))
                    y.append(i * np.sin(th))
                    z.append(zpa)
                ax1.plot3D(x, z, y, fmt)
            count += 1
            
        # comentar la parte de abajo hasta el ax1.plot3D para quitar parte del tracker
        # for i in self.ri:
        #     xt = []
        #     yt = []
        #     # if i % 2 == 1:
        #     for th in theta:
        #         x = []
        #         y = []
        #         z = []
        #         xt.append(i * np.cos(th))
        #         yt.append(i * np.sin(th))
        #         for zpa in zpace:            
        #             x.append(i * np.cos(th))
        #             y.append(i * np.sin(th))
        #             z.append(zpa)
        #         ax1.plot3D(x, z, y, fmt)

        # ax2 = fig.add_subplot(1,3,1)
        for i in self.ri:
            xt = []
            yt = []
            for th in theta:
                xt.append(i * np.cos(th))
                yt.append(i * np.sin(th))
            ax2.plot(xt, yt, fmt)
        
        # ax3 = fig.add_subplot(1,3,3)
        for i in self.ri:
            zt = []
            yt = []
            yt2 = []
            for zpa in zpace:
                zt.append(zpa)
                yt.append(i)
                yt2.append(-i)
                
            ax3.plot(zt, yt, fmt, zt, yt2, fmt)
          
    def convertAngle(self, theta):
        if theta > np.pi:
            return theta - np.pi * 2.0
        elif theta < -np.pi:
            return theta + np.pi * 2.0
        return theta 
 
    def intersection(self, track):
        '''
        Computes the interesection between the tracker and a set of given tracks

        Parameters
        ----------
        tracks : list
            list of Track.

        Returns
        -------
        x_c : float
            x cartesian coordinate of the intersection point.
        y_c : float
            y cartesian coordinate of the intersection point.
        z_c : float
            z cartesian coordinate of the intersection point.

        '''
        # print(len(tracks))
        # for track in tracks:
            
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
        t = []
        for j, troot in enumerate(rootx):
            if rootx[j] < 0:
                if rooty[j] < 0:
                    raise ValueError('Problem with one intersection')
                else:       
                    y1_c[j] = (-by[j] - np.sqrt(rooty[j])) / (2.0 * ay)
                    y2_c[j] = (-by[j] + np.sqrt(rooty[j])) / (2.0 * ay)
                    if abs((y1_c[j] - track.y_c) / track.rt) > 1.0:
                        raise ValueError('Problem with one intersection')
                    if abs((y2_c[j] - track.y_c) / track.rt) > 1.0:
                        raise ValueError('Problem with one intersection')
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
                    t.append(tv)               
            else:
                x1_c[j] = (-bx[j] - np.sqrt(rootx[j])) / (2.0 * ax)
                x2_c[j] = (-bx[j] + np.sqrt(rootx[j])) / (2.0 * ax)
                if abs((x1_c[j] - track.x_c) / track.rt) > 1.0:
                    raise ValueError('Problem with one intersection')
                if abs((x2_c[j] - track.x_c) / track.rt) > 1.0:
                    raise ValueError('Problem with one intersection')
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
                t.append(tv)
        x, y, z = track.eval(np.array(t))    
        return x, y, z
    
    
    def plot_intersection(self, track, fmt = 'r*'):
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
        self.plot_points(x, y, z, fmt)
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

    def plot_points(self, x, y, z, fmt):

        ax2.plot(x, y, fmt)
        ax1.plot3D(x, z, y, fmt)
        ax3.plot(z, y, fmt)
    
    def plot_meast_points(self, track, fmt = 'gx'):
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
        self.plot_points(x, y, z, fmt)
        return x, y, z    

class Track:
    
    def __init__(self, dxy, dz, phi, eta, pt, q):
        '''

        Parameters
        ----------
        rt : float
            desired curvature radius.
        vertex_x : float, optional
            closest x coordinate to the vertex point, or closest approximation 
            point. The default is 0.
        vertex_y : float, optional
            closest y coordinate to the vertex point, or closest approximation 
            point. The default is 0.
        z_0 : float, optional
            initial z value for the track. The default is 0.

        Returns
        -------
        None.

        '''

        ##################### Parametros de traza conocido el vertice #######################
        ##### x(0) -> X position of the vertex
        ##### y(0) -> Y position of the vertex
        ##### z(0) -> Z position of the vertex
        ##### phi -> Direction of the particle in the transverse plane
        ##### eta -> Direction of the particle in the longitudinal plane 
        ##### pt -> transverse momentum
        ##### charge -> the charge
        #################################################################
        ##### In the z coordinate we have: pz = pt * sinh(eta) and z = pt * sinh(eta) / m * t + z(0)
        #################################################################
        ##### In the transverse plane we have:
        ##### Px(t) = Pt cos(wt - phi) and Py(t) = -Pt sin(wt - phi)
        ##### x(t) = Pt/qB [sin(wt - phi) + sin(phi)] + x(0)
        ##### y(t) = Pt/qB [cos(wt - phi) - cos(phi)] + y(0)
        #################################################################
        ##### Relation to the equation of the cirumference:
        ##### (x - x(0) - Pt/qB sin(phi))^2 + (y - y(0) + Pt/qB cos(phi))^2 = (Pt/qB)^2 
        ##### Therefore:
        ##### Radius = Pt/mw
        ##### Center of the circumference X = x(0) + Pt/qB sin(phi)
        ##### Center of the circumference Y = y(0) - Pt/qB cos(phi)
        #################################################################

        ######### Parametros de traza en funcion del punto de maxima aproximacion #############
        ##### dxy -> Distance in the transverse plane from the POCA to the beam line
        ##### dz  -> Distance in z from the the POCA to the origin 
        ##### phi -> Direction of the particle in the transverse plane at the POCA
        ##### eta -> Direction of the particle in the longitudinal plane 
        ##### pt -> transverse momentum
        ##### charge -> the charge
        ##### With this we can easily produce the track
        #################################################################

        ######### How to use the POCA point as vertex #############
        ##### The point of closest approach is that in which the position and momentum
        ##### in the transverse plane are orthogonal (x, y) * (px, py) = 0
        ##### This implies that xpoca * px + ypoca * py = 0 => xpoca * px = -ypoca * py => xpoca / ypoca = -py / px
        ##### But py/px is by the definition the direction of the track at the POCA
        ##### Therefore: xpoca/ypoca = -tan(phi) 
        ##### If we write: xpoca = dxy * cos(phi_poca) and ypoca = dxy * sin(phi_poca) we have that:
        ##### 1 / tan(phi_poca) = -tan(phi) => cotan(phi_poca) = -tan(phi)
        ##### With this angle we can estimate: xpoca = dxy * cos(phi_poca) and ypoca = dxy*sin(phi_poca)
        ##### We can assume that t = 0 is the POCA point and then use the parametrization above 
        ##### pz = pt * sinh(eta) and z = pt * t + z(0) = pt * t + dz
        ##### x(t) = Pt/qB [sin(wt - phi) + sin(phi)] + x(0) = Pt/qB [sin(wt - phi) + sin(phi)] + dxy * cos(phi_poca)
        ##### y(t) = Pt/qB [cos(wt - phi) - cos(phi)] + y(0) = Pt/qB [cos(wt - phi) - cos(phi)] + dxy * sin(phi_poca)

        ##### rho = Pt/qB if we set Pt in GeV, B in Tesla and rho in cm is: (Pt [GeV]/ B[Tesla]) * 1000.0 / 3.0  


        self.dxy = dxy
        self.dz = dz
        self.phi = phi
        self.eta = eta
        self.pt = pt
        self.q = np.sign(q)
        #self.q = q
        self.pz = pt * np.sinh(eta) 

        #Some constants, maybe wise to put them somewhere else
        self.m = 0.140
        b = 3.8
        self.gamma = np.sqrt(self.pt**2 + self.pz**2 + self.m**2)/self.m
        self.w = self.q * 2.7143 / self.gamma

        #POCA points
        self.x_cp = dxy * np.cos(phi + np.pi/2.0)
        self.y_cp = dxy * np.sin(phi + np.pi/2.0)
        self.z_cp = self.dz       

        #Circle description
        self.rt = self.q * pt / b * (1000.0/3.0)
        self.x_c = self.rt * np.sin(phi) + self.x_cp
        self.y_c = -self.rt * np.cos(phi) + self.y_cp
        
        # zpace = np.linspace(-self.zsize / 2.0, self.zsize / 2.0, len(theta)) 
        
    def eval(self, t):
        
        x = self.rt * np.sin(self.w * t - self.phi) + self.x_c       
        y = self.rt * np.cos(self.w * t - self.phi) + self.y_c
        z = 30.0 * self.pz / (self.gamma*self.m) * t + self.dz        
        return x, y, z

    def plot_track(self, fmt = 'r'):
        '''
        Plots track

        Parameters
        ----------
        fmt : string, optional
            describes the colour and style of the track line. The default is 'r'.

        Returns
        -------
        None.

        '''
        # start = np.arctan2(self.displ_y, self.displ_x)
        # theta = np.linspace(0, 2 * np.pi, 100)
        # theta_p = np.linspace(start + np.pi, start + 1 , 100)
        # x = self.rt * np.cos(theta_p) + self.displ_x
        # y = self.rt * np.sin(theta_p) + self.displ_y
        # fig = plt.figure(figsize = plt.figaspect(0.3))
        
        # plt.plot(self.x, self.y, fmt, label = 'Real track')
        # ax2 = fig.add_subplot(1,3,1)
        
        #We have all the ingredients, we just need to propagate the track for a given time
        t = np.linspace(0, 5, 200)
        x, y, z = self.eval(t)
        ax1.plot3D(x, z, y, fmt)
        ax2.plot(x, y, fmt, label = 'Real track')
        ax2.plot(0, 0, 'rx')
        ax3.plot(z, y, fmt)
        # plt.plot(self.displ_x, self.displ_y, 'go')
    
    
        
#%%
from statsmodels.base.model import GenericLikelihoodModel

class MyLikelihood(GenericLikelihoodModel):

     def __init__(self, endog, exog, tracker, **kwds):
         self.tracker = tracker
         self.chi2 = 0
         super(MyLikelihood, self).__init__(endog, exog, **kwds)

     def loglike(self, params):

         dxy = params[0]
         dz = params[1]
         px = params[2]
         py = params[3]
         pz = params[4]
         #phi = params[2]
         #eta = params[3]
         #pt = params[4]
         #ch = np.sign(params[5])
         #ch = params[5]
         phi = np.arctan2(py, px)
         pt = np.sqrt(px*px+py*py) 
         theta = np.arctan2(pt, pz) 
         eta = -np.log(np.tan(theta/2.0))
         myTrack = Track(dxy, dz, phi, eta, pt, 1.0)
         try:
              x, y, z = self.tracker.intersection(myTrack)
         except ValueError as err:   
              return -1.1*self.chi2 
  
         xmeas, ymeas, zmeas = self.endog
         
         self.chi2 = 0
         #print(dxy, dz, phi, eta, pt)
 
         for i in range(0, len(x)):
             chi2Z = (zmeas[i] - z[i])**2 / self.tracker.z_error[i]**2
             normmeas = np.sqrt(xmeas[i]*xmeas[i]+ymeas[i]*ymeas[i])
             norm = np.sqrt(x[i]*x[i]+y[i]*y[i])
             product = xmeas[i] * x[i] + ymeas[i] * y[i]
             anglediff = np.arccos(product/(norm * normmeas))
             chi2Transverse = (anglediff)**2 / self.tracker.phi_error[i] ** 2
             self.chi2 += (chi2Z + chi2Transverse) 
             #print((zmeas[i]-z[i]), self.tracker.z_error[i])     
         #print(-self.chi2)
         return -self.chi2


#%%





#Some global variables
fig = plt.figure(figsize = plt.figaspect(0.3))
ax1 = fig.add_subplot(1, 3, 2, projection = '3d')
ax2 = fig.add_subplot(1, 3, 1)
ax3 = fig.add_subplot(1, 3, 3)


if __name__ == "__main__":

    #configuring the tracker
    layers = np.linspace(1, 100, 20)
    sigma_phi = 2 / layers
    sigma_z = 4 + np.zeros(len(layers))
    tracker = Tracker(layers, sigma_phi, sigma_z, 600, [0,0,0])
    tracker.plot_tracker()

    #An example track
    track = Track(0, 0, np.pi/2.0, 1.0, 15.0, 1.0)
    track.plot_track()
    try: 
        x, y, z = tracker.plot_intersection(track)
    except ValueError as err:
        sys.exit()

    #sys.exit()
    meas = tracker.plot_meast_points(track)

    #Let's help the minimizer by setting good initial parameters
    dxy_init = 0
    dz_init = 0
    phi_init = np.arctan(meas[1][0]/meas[0][0])
    theta_init = np.arctan( np.sqrt(meas[0][0]**2 + meas[1][0]**2)/meas[2][0])
    eta_init = -np.log(np.tan(theta_init / 2.0))

    externalPoint = np.array([meas[0][len(layers)-1], meas[1][len(layers)-1]])
    middlePoint = np.array([meas[0][math.floor(len(layers)/2.0)], meas[1][math.floor(len(layers)/2.0)]])
    pointZero = np.array([0, 0])
    
    middleFlat = (pointZero + externalPoint)/2.0
    h = np.linalg.norm(middleFlat - middlePoint)
    d = np.linalg.norm(middlePoint - pointZero)
    l = np.linalg.norm(middleFlat - pointZero)
    s = (l**2 - h**2)/(2.0*h)
    r = s + h
    estimated_pt = 3.0*3.8/1000.0 * r
   
    print('dxy: ', track.dxy, dxy_init)
    print('dz: ', track.dz, dz_init)
    print('phi: ', track.phi, phi_init)
    print('eta: ', track.eta, eta_init)
    print('pt: ', track.pt, estimated_pt)
    
    px_init = estimated_pt * np.cos(phi_init)
    py_init = estimated_pt * np.sin(phi_init)
    pz_init = estimated_pt * np.sinh(phi_init)

    l = MyLikelihood(meas, np.array([1, 1, 1]), tracker)

    result = l.fit(start_params = (dxy_init, dz_init, px_init, py_init, pz_init))

    print(result.mle_settings)
    print(result.mle_retvals)

    par = result.params
    dxy = par[0]
    dz = par[1]
    phi = np.arctan2(par[3], par[2])
    pt = np.sqrt(par[2]**2+par[3]**2) 
    theta = np.arctan2(pt, par[4]) 
    eta = -np.log(np.tan(theta/2.0))
    #q = par[5]
    print(dxy, track.dxy)
    print(dz, track.dz)
    print(phi, track.phi)
    print(eta, track.eta)
    print(pt, track.pt)
    #print(q, track.q)
    
    track2 = Track(dxy, dz, phi, eta, pt, 1.0)
    track2.plot_track('b')

    plt.show()



 
