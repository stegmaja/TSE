from copy import copy
class StepInterp:
    def __init__(self,x,y):
        self.x = copy(x)
        self.y = copy(y)
    def interp1d(self,x0):
        if(x0<=self.x[0]):
            return int(self.y[0])
        elif(x0>=self.x[-1]):
            return int(self.y[-1])
        else:
            return int(self.y[self.x<=x0][-1])