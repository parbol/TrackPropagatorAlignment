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
        mat = [[cosphi*costheta, -sinphi, cosphi*sintheta],
               [sinphi*costheta, cosphi, sinphi*sintheta],
               [-sintheta, 0.0, costheta]]
       
        self.rot = np.asmatrix(mat)
        self.invrot = np.linalg.inv(self.rot)


    def fmin(self, t, a, b, args, track):
        A = self.n[0]
        B = self.n[1]
        C = self.n[2]
        D = -(A*self.p[0]+B*self.p[1]+C*self.p[2])
        r = args[0]
        w = args[1]
        phi = args[2]
        pzct = args[3]
        x0 = args[4]
        y0 = args[5]
        z0 = args[6]
        Delta = A * x0 + B * y0 + C * z0
        return Delta + D + A * r * sin(w*t-phi) + B * r * cos(w*t-phi) + pzct * t
    


    def intersection(self, track):

        args = (track.rt, track.w)

        trackpoint = np.asarray([track.x_cp, track.y_cp, track.z_cp])
        trackMomentum = np.asarray([track.pt*np.cos(track.phi),
                                    track.pt*np.sin(track.phi),
                                    track.pz])
        newtrackpoint = np.asarray(self.invrot.dot(trackpoint))[0]
        newtrackMomentum = np.asarray(self.invrot.dot(trackMomentum))[0]

        
        newp = np.asarray((self.invrot.dot(self.p)))[0]
        ztouch = newp[2]
        t = (track.gamma*track.m)/(29.98*newtrackMomentum[2])*(ztouch-newtrackpoint[2])
        x,y,z = track.eval(t)
        #np.concatenate((track.xi, np.asarray(x)), axis=0)
        #np.concatenate((track.yi, np.asarray(y)), axis=0)
        #np.concatenate((track.zi, np.asarray(z)), axis=0)
        #np.concatenate((track.ti, np.asarray(t)), axis=0)
        #track.xi = np.append(track.xi, [x])
        #track.yi = np.append(track.yi, [y])        
        #track.zi = np.append(track.zi, [z])        
        #track.ti = np.append(track.ti, [t])        
        
        return x, y, z, t

    def intersectionStraight(self, x0, y0, z0, vx, vy, vz):

        A = self.n[0]
        B = self.n[1]
        C = self.n[2]
        D = -(A*self.p[0]+B*self.p[1]+C*self.p[2])
        Delta = A * x0 + B * y0 + C * z0
        Kapa = A * vx + B * vy + C * vz
        t = (-Delta-D)/Kapa
        x = x0 + vx * t
        y = y0 + vy * t
        z = z0 + vz * t
        return x, y, z
    
    def belongsToPlane(self, x, y, z):

        A = self.n[0]
        B = self.n[1]
        C = self.n[2]
        D = -(A*self.p[0]+B*self.p[1]+C*self.p[2])
        k = A * x + B * y + C *z + D
        if np.abs(k) < 1e-3:
            return True
        return False
