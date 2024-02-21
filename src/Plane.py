import numpy as np
import sys



class Plane:

    def __init__(self, x0, y0, z0, nx, ny, nz):
        self.p = np.asarray([x0, y0, z0])
        self.n = np.asarray([nx, ny, nz])
        s = self.norm(self.n)
        if s < 1e-5:
            print('Bad plane definition')
            sys.exit()
        self.n = self.n / s


    def norm(self, v):
        return np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    