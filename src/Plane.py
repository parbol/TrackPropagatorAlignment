import numpy as np
import sys

from src.Track import Track

class Plane:

    def __init__(self, x0, y0, z0, nx, ny, nz):
        
        self.p = np.asarray([x0, y0, z0])
        self.n = np.asarray([nx, ny, nz])
        s = self.norm(self.n)
        if s < 1e-5:
            print('Bad plane definition')
            sys.exit()
        self.n = self.n / s
        self.rotMatrix()


    def norm(self, v):
        return np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    
    def phi(self):

        return np.arctan2(self.n[1], self.n[0])
    
    def theta(self):

        return np.arccos(self.n[2])
    
    def rotMatrix(self):
        
        cosphi = np.cos(self.phi())
        sinphi = np.sin(self.phi())
        costheta = np.cos(self.theta())
        sintheta = np.sin(self.theta())
        mat = [[cosphi*costheta, sinphi*costheta, -sintheta],
               [-sinphi, cosphi, 0],
               [cosphi*sintheta, sinphi*sintheta, costheta]]
              
        self.rot = np.asmatrix(mat)
        self.invrot = np.linalg.inv(self.rot)

    def intersection(self, track):

        trackpoint = np.asarray([track.x_cp, track.y_cp, track.z_cp])
        trackMomentum = np.asarray([track.pt*np.cos(track.phi),
                                    track.pt*np.sin(track.phi),
                                    track.pz])
        newtrackpoint = np.asarray(self.invrot.dot(trackpoint))[0]
        newtrackMomentum = np.asarray(self.invrot.dot(trackMomentum))[0]

        print(self.rot)

        newp = np.asarray((self.invrot.dot(self.p)))[0]
        ztouch = newp[2]

        t = (track.gamma*track.m)/(30.0*newtrackMomentum[2])*(ztouch-trackpoint[2])
        x,y,z = track.eval(t)

        return np.asarray([x,y,z])
        

    