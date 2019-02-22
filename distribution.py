import numpy as np
import internal_forces as ifo

class Slice(object):

    def __init__(self, x, dx):
        self.x = x
        self.dx = dx
        self.vx = ifo.Force(1, np.array([1,0,0]), np.array([x, 0, 0]))
        self.vy = ifo.Force(1, np.array([0,1,0]), np.array([x, 0, 0]))
        self.vz = ifo.Force(1, np.array([0,1,0]), np.array([x, 0, 0]))
        self.mx = 0
        self.my = 0
        self.mz = 0

