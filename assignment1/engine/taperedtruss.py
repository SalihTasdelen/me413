import numpy as np
from .truss import TrussLinear, TrussQuadratic

class TaperedTrussLinear(TrussLinear):

    def __init__(self, E: float, L: float, A0: float, h1: float, h2: float, nElem: int) -> None:
        super().__init__(E, L, A0, nElem)
        self.h1 = h1
        self.dh = h2 - h1

    def getLocalK(self, e: int):
        x1, x2 = self.x_node[self.conn[e]]
        c = (1 + self.dh / self.L / self.h1 * (x1 + x2) / 2)
        return self.localK * c
    
class TaperedTrussQuadratic(TrussQuadratic):

    def __init__(self, E: float, L: float, A0: float, h1: float, h2: float, nElem: int) -> None:
        super().__init__(E, L, A0, nElem)
        self.localK = self.localK
        self.k = self.k / 10
        self.h1 = h1
        self.dh = h2 - h1
        self.localKx1 = self.k * np.array(
            [[ 37, -44,  7],
             [-44,  48, -4],
             [  7,  -4, -3]])
        self.localKx2 = self.k * np.array(
            [[ 36, -32,  -4],
             [-32,  64, -32],
             [ -4, -32,  36]])
        self.localKx3 = self.k * np.array(
            [[ -3,  -4,   7],
             [ -4,  48, -44],
             [  7, -44,  37]])


    def getLocalK(self, e: int):
        x1, x2, x3 = self.x_node[self.conn[e]]
        c = self.dh / self.h1 / self.L
        return self.localK + c * (x1 * self.localKx1 + x2 * self.localKx2 + x3 * self.localKx3)