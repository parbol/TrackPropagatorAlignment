


class BTLId:
    
    def __init__(self):

        self.tray = -1
        self.side = 0
        self.RUType = -1
        self.RUNumber = -1
        self.RU = -1
        self.module = -1

    def setTray(self, tray):
        self.tray = tray
    
    def setSide(self, side):
        self.side = side

    def setRU(self, type, number):
        self.RUType = type
        self.RUNumber = number
        self.RU = self.RUType * 2 + self.RUNumber

    def setModule(self, module):
        self.module = module

        
        









