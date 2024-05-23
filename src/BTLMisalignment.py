from src.Plane import Plane
from src.Module import Module
from src.BTLTray import BTLTray
from src.EulerRotation import EulerRotation
from src.BTLId import BTLId
import numpy as np
import sys


class BTLMisalignment:
    
    def __init__(self, btl, misalignmentName, btlName):

        self.btl = btl
        self.misalignmentName = misalignmentName
        self.btlName = btlName


    def randomLocalMisalignment(self, sigmax, sigmay):

        trays = self.btl.btlNominal.pTrays + self.btl.btlNominal.mTrays
        f = open(self.misalignmentName, 'w')
        for tray in trays:
            for ru in tray.RUs:
                for module in ru.Modules:
                    displacementX = np.random.normal(0, sigmax)
                    displacementY = np.random.normal(0, sigmay)
                    disp = np.asarray([displacementX, displacementY, 0.0])
                    newpos = module.module.toGlobal(disp)
                    module.updatePosition(newpos, module.module.eulerAngles)
                    cad = '{side} {tray} {RUType} {RUNumber} {module} {dx} {dy}'.format(side = str(module.btlId.side),
                                                                  tray = str(module.btlId.tray),
                                                                  RUType = str(module.btlId.RUType),
                                                                  RUNumber = str(module.btlId.RUNumber),
                                                                  module = str(module.btlId.module),
                                                                  dx = displacementX,
                                                                  dy = displacementY)
                    f.write(cad + '\n')
        f.close()
        self.btl.btlNominal.writeGeometry(self.btlName)
