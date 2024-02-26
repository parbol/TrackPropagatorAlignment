import math
import numpy as np
import matplotlib.pyplot as plt
from random import Random
import sys

from src.Tracker import Tracker
from src.Track import Track
from src.Plane import Plane
from src.Module import Module
from src.ETL import ETL





if __name__ == "__main__":

    #Some global variables
    fig = plt.figure(figsize = plt.figaspect(0.3))
    ax1 = fig.add_subplot(1, 3, 2, projection = '3d')
    ax2 = fig.add_subplot(1, 3, 1)
    ax3 = fig.add_subplot(1, 3, 3)
    #configuring the tracker
    layers = np.linspace(1, 100, 20)
    layersz = np.linspace(130, 270, 5)
    sigma_phi = 2 / layers
    sigma_z = 4 + np.zeros(len(layers))
    tracker = Tracker(layers, layersz, sigma_phi, sigma_z, 220.0, [0,0,0])
    tracker.plot_tracker(ax1, ax2, ax3)

    #An example track
    
    track = Track(0, 0, np.pi/2.0, 1.7, 10, 1.0)
    p = Plane(0.0, 0.0, 400.0, 0.0, 0.0, 1.0)
    
    m = ETL(350.0, 4.0, 4.0, 30, 50, 1.0, 1.0, 4.0, 4.0, 30.0, 127.0)
    m.draw(ax1, ax2, ax3, 'g')

    x, y, z = tracker.plot_intersection(track, ax1, ax2, ax3, 'g*')
    x, y, z = p.intersection(track)

    track.plot_track(ax1, ax2, ax3, 'black')
    #track2 = Track(0, 0, np.pi/2.0, 0.5, 3.0, 1.0)
    #x, y, z = tracker.plot_intersection(track2, ax1, ax2, ax3, 'g*')
    #track2.plot_track(ax1, ax2, ax3)
    #sys.exit()
    #meas = tracker.plot_meast_points(track)

    #Let's help the minimizer by setting good initial parameters
    #dxy_init = 0
    #dz_init = 0
    #phi_init = np.arctan(meas[1][0]/meas[0][0])
    #theta_init = np.arctan( np.sqrt(meas[0][0]**2 + meas[1][0]**2)/meas[2][0])
    #eta_init = -np.log(np.tan(theta_init / 2.0))

    #externalPoint = np.array([meas[0][len(layers)-1], meas[1][len(layers)-1]])
    #middlePoint = np.array([meas[0][math.floor(len(layers)/2.0)], meas[1][math.floor(len(layers)/2.0)]])
    #pointZero = np.array([0, 0])
    
    #middleFlat = (pointZero + externalPoint)/2.0
    #h = np.linalg.norm(middleFlat - middlePoint)
    #d = np.linalg.norm(middlePoint - pointZero)
    #l = np.linalg.norm(middleFlat - pointZero)
    #s = (l**2 - h**2)/(2.0*h)
    #r = s + h
    #estimated_pt = 3.0*3.8/1000.0 * r
   
    #print('dxy: ', track.dxy, dxy_init)
    #print('dz: ', track.dz, dz_init)
    #print('phi: ', track.phi, phi_init)
    #print('eta: ', track.eta, eta_init)
    #print('pt: ', track.pt, estimated_pt)
    
    #px_init = estimated_pt * np.cos(phi_init)
    #py_init = estimated_pt * np.sin(phi_init)
    #pz_init = estimated_pt * np.sinh(phi_init)

    #l = MyLikelihood(meas, np.array([1, 1, 1]), tracker)

    #result = l.fit(start_params = (dxy_init, dz_init, px_init, py_init, pz_init))

    #print(result.mle_settings)
    #print(result.mle_retvals)

    #par = result.params
    #dxy = par[0]
    #dz = par[1]
    #phi = np.arctan2(par[3], par[2])
    #pt = np.sqrt(par[2]**2+par[3]**2) 
    #theta = np.arctan2(pt, par[4]) 
    #eta = -np.log(np.tan(theta/2.0))
    #q = par[5]
    #print(dxy, track.dxy)
    #print(dz, track.dz)
    #print(phi, track.phi)
    #print(eta, track.eta)
    #print(pt, track.pt)
    #print(q, track.q)
    
    #track2 = Track(dxy, dz, phi, eta, pt, 1.0)
    #track2.plot_track('b')

    plt.show()



 
