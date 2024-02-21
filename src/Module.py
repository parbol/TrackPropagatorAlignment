from src.Plane import Plane


class Module:

    def __init__(self, x0, y0, z0, nx, ny, nz, Lx, Ly):

        self.plane = Plane(x0,y0,z0,nx,ny,nz)
        self.pA = np.asarray([x0 - Lx/2.0, y0 - Ly/2.0])
