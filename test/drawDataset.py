import numpy as np
import torch 
from torch_geometric.data import Data
import matplotlib.pyplot as plt
import optparse
import math
from src.TrackerCarlos import TrackerCarlos


      
def drawGraph(dataset, ax1, ax2, ax3, ax4, alpha=0.1):
    
    x = dataset['source'].x
    
    
    #ax.plot(x[:,0].numpy(), x[:,1].numpy(), 'g*')
    ax1.plot(x[:,0].numpy(), x[:,1].numpy(), x[:,2].numpy(), 'y*')
    ax2.plot(x[:,0].numpy(), x[:,1].numpy(), 'y*')
    ax2.plot(x[:,2].numpy(), x[:,1].numpy(), 'y*')
    ax2.plot(x[:,2].numpy(), x[:,0].numpy(), 'y*')


    edge_index = dataset['source', 'weight', 'target'].edge_index
    edge_label = dataset['source', 'weight', 'target'].edge_label

    for i, edge in enumerate(torch.t(edge_index)):
       
        edge1 = edge[0].numpy()
        edge2 = edge[1].numpy()
        if edge_label[i] > 0.5:
            x1, y1, z1 = x[edge1,0], x[edge1,1], x[edge1, 2]
            x2, y2, z2 = x[edge2,0], x[edge2,1], x[edge2, 2]
            xg = [] 
            xg.append(x1)
            xg.append(x2)
            yg = []
            yg.append(z1)
            yg.append(z2)
            zg = []
            zg.append(y1)
            zg.append(y2)
            #ax.plot(xg, yg, 'r-')    
            ax1.plot(xg, yg, zg, 'r-')
            ax2.plot(xg, yg, 'r-')
            ax3.plot(zg, yg, 'r-')
            ax4.plot(zg, xg, 'r-')




if __name__=='__main__':

    parser = optparse.OptionParser(usage='usage: %prog [options] path', version='%prog 1.0')
    parser.add_option('-i', '--input', action='store', type='string', dest='inputFile', default='input.pt', help='Input Reference Dataset')
    (opts, args) = parser.parse_args()
    
    
    
    dataset = torch.load(opts.inputFile)
    
    #Some global variables
    #fig = plt.figure(figsize = plt.figaspect(0.3))
    fig = plt.figure(figsize = (8, 8), layout="constrained")
    gs0 = fig.add_gridspec(2, 1, height_ratios=[2,1])
    ax1 = fig.add_subplot(gs0[0], projection = '3d')
    ax1.grid(False)
    #ax1.set_axis_off()
    gs1 = gs0[1].subgridspec(1,3)
    ax2 = fig.add_subplot(gs1[0])
    ax3 = fig.add_subplot(gs1[1])
    ax4 = fig.add_subplot(gs1[2])
    ax1.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax1.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax1.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax1.set_xlabel('x [cm]')
    ax1.set_ylabel('z [cm]')
    ax1.set_zlabel('y [cm]')
    ax2.set_xlabel('x [cm]')
    ax2.set_ylabel('y [cm]')
    ax3.set_xlabel('z [cm]')
    ax3.set_ylabel('y [cm]')
    ax4.set_xlabel('z [cm]')
    ax4.set_ylabel('x [cm]')

    #configuring the tracker
    layers = np.asarray([10, 20, 30, 40, 50, 60, 70, 80, 90])
    N = [7, 8, 9, 10, 11, 12, 13, 14, 15]
    sigma_rphi = 0.01
    sigma_z = 0.01
    tracker = TrackerCarlos(layers, 200, N, 0, 200, sigma_rphi, sigma_z, 0.01)
    tracker.plot_tracker(ax1, ax2, ax3, ax4)

    drawGraph(dataset, ax1, ax2, ax3, ax4)
    
    plt.show()

