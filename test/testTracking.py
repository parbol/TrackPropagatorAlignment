import math
import numpy as np
import matplotlib.pyplot as plt
from random import Random
import sys
import h5py
import pandas as pd


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
    ax1.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax1.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax1.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    
    #configuring the tracker
    layers = np.linspace(1, 100, 20)
    layersz = np.linspace(130, 270, 5)
    sigma_rphi = 0.01
    sigma_z = 0.01
    tracker = Tracker(layers, layersz, sigma_rphi, sigma_z, 220.0, [0,0,0])
    tracker.plot_tracker(ax1, ax2, ax3)
    
    #Configuring the ETL 
    etl = ETL(350.0, 4.0, 4.0, 30, 50, 4.0, 1.0, 2.0, 2.0, 30.0, 127.0, 0.1, 0.001, 0.3)
    etl.draw(ax1, ax2, ax3, 'g')
    
    
    #An example track 
    #track = Track(0, 0, np.pi/2.0, 2.1, 10, 1.0)
    #tracker.fullMeasurement(track)
    #valid = etl.fullMeasurement(track)
    counter = 0
    alist = []
    while counter < 10000:
        phi = np.random.uniform(0, 2.0*np.pi)
        eta = np.random.uniform(1.6, 3.0)
        pt = np.random.uniform(1.0, 50.0)
        track = Track(0, 0, phi, eta, pt, 1.0)
        tracker.fullMeasurement(track)
        valid = etl.fullMeasurement(track)
        if valid:
            n = len(track.xi)
            xi = track.xi[n-1]
            yi = track.yi[n-1]
            zi = track.zi[n-1]
            ti = track.ti[n-1]
            xin = track.xin[n-1]
            yin = track.yin[n-1]
            zin = track.zin[n-1]
            tin = track.tin[n-1]
            a = [xi, yi, zi, ti, xin, yin, zin, tin]
            alist.append(a)
            counter = counter + 1

    df = pd.DataFrame(np.asmatrix(alist))
    df.to_hdf('datos.h5', '/data/d1')
    
    #plt.show()



 
