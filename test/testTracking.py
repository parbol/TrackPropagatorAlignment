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
from src.BTL import BTL






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
    #etl = ETL(350.0, 4.0, 4.0, 30, 50, 4.0, 1.0, 2.0, 2.0, 30.0, 127.0, 0.1, 0.001, 0.3)
    #etl.draw(ax1, ax2, ax3, 'g')
    
    #Configuring the BTL
    R = 120.0
    TrayLength = 300.0
    TrayWidth = 2.0*R*np.sin(3.0*np.pi/180.0)
    TrayStartZ = 1.0
    TrayStartPhi = 5.0*np.pi/180.0
    RULength = 45.0
    ModuleLength = 5.4
    ModuleWidth = 4.0
    rphi_error = 0.1
    z_error = 0.1
    t_error = 0.1
    btl = BTL(R, TrayLength, TrayWidth, TrayStartZ, TrayStartPhi, RULength, ModuleLength, ModuleWidth, rphi_error, z_error, t_error, 9.4)
    btl.draw(ax1, ax2, ax3, 'g')
    #etl.draw(ax1, ax2, ax3, 'g')

    #An example track 
    #track = Track(0, 0, np.pi/2.0, 0.1, 10, 1.0)
    #tracker.fullMeasurement(track)
    #valid = btl.fullMeasurement(track)
    #track.plot_track(ax1, ax2, ax3, 'r')
    counter = 0
    alist = []
    while counter < 10:
        phi = np.random.uniform(0, 2.0*np.pi)
        eta = np.random.uniform(-1.6, 1.6)
        pt = np.random.uniform(1.0, 2.0)
        charge = np.sign(np.random.uniform(-1.0, 1.0))
        track = Track(0, 0, phi, eta, pt, charge)
        tracker.fullMeasurement(track)
        valid = btl.fullMeasurement(track)
        if valid:
            track.plot_track(ax1, ax2, ax3, 'r')
            counter = counter + 1
           
    
    plt.show()



 
