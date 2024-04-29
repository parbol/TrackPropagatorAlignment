import numpy as np

from src.BTL import BTL
from src.BTLMisalignment import BTLMisalignment



if __name__ == "__main__":
 
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
    misa = BTLMisalignment(btl, 'misalignemnts.txt', 'btlUpdatedGeometry.txt')
    misa.randomLocalMisalignment(0.5, 0.5)
