import numpy as np
from .fem1d import FEM1D

class Truss(FEM1D):

    def __init__(self, E:float, L: float, A: float, nElem: int, nLocalNode: int) -> None:
        """_summary_

        Args:
            E (float): _description_
            L (float): _description_
            A (float): _description_
            nElem (int): _description_
        """        
        super().__init__(L, nElem, nLocalNode)
        self.E = E
        self.A = A

    def getLocalK(self, e: int):
        """Computes the Stiffness matrix for the given element.

        Args:
            e (int): Element ID.

        Raises:
            NotImplementedError: Base Class does not implement the Stiffness Matrix.
        """        
        raise NotImplementedError
    
    def getXField(self, n:int = 20):
        raise NotImplementedError
    
    def getDispField(self, n:int = 20):
        raise NotImplementedError
    
    def getStrainField(self, n:int = 20):
        raise NotImplementedError

    def formKGlobal(self):
        KGlobal = np.zeros((self.nNode, self.nNode))
        for i in range(self.nElem):
            K = self.getLocalK(i)
            KShape = K.shape
            idx = self.conn[i,0]
            KGlobal[idx:idx+KShape[0], idx:idx+KShape[1]] += K
        self.system.KGlobal = KGlobal
        return KGlobal

class TrussLinear(Truss):

    
    def __init__(self, E: float, L: float, A: float, nElem: int) -> None:
        super().__init__(E, L, A, nElem, nLocalNode = 2)
        self.k = self.E * self.A / self.lElem
        self.localK = np.array([[self.k, -self.k], [-self.k, self.k]])

    def getLocalK(self, e: int):
        return self.localK

    def formFGlobal(self, qStart, qEnd, fStart, fEnd):
        FGlobal = np.zeros((self.nNode,1))
        qDiff = (qEnd - qStart) / self.L * self.lElem
        q = np.linspace(qStart, qEnd, self.nNode)

        for i, nodes in enumerate(self.conn):
            FGlobal[nodes[0]] += self.lElem * (q[i] / 2 + qDiff / 6)
            FGlobal[nodes[1]] += self.lElem * (q[i] / 2 + qDiff / 3)
        FGlobal[-1] += fEnd
        FGlobal[0] += fStart
        
        self.system.FGlobal = FGlobal

    def getXField(self, n: int = 20):
        N = n * self.nElem
        x = np.zeros((N, 1))
        ksi = np.linspace(-1, 1, n).reshape((n,1))
        N1 = (1 - ksi) / 2
        N2 = (1 + ksi) / 2
        for i, (node1, node2) in enumerate(self.conn):
            N1x1 = self.x_node[node1] * N1
            N2x2 = self.x_node[node2] * N2
            x[i*n:(i+1)*n] = N1x1 + N2x2
        return x
    
    def getDispField(self, n: int = 20):
        N = n * self.nElem
        uField = np.zeros((N, 1))
        ksi = np.linspace(-1, 1, n).reshape((n,1))
        N1 = (1 - ksi) / 2
        N2 = (1 + ksi) / 2
        for i, (node1, node2) in enumerate(self.conn):
            N1x1 = self.solution[node1] * N1
            N2x2 = self.solution[node2] * N2
            uField[i*n:(i+1)*n] = N1x1 + N2x2
        return uField
    
    def getStrainField(self, n: int = 20):
        N = n * self.nElem
        sField = np.zeros((N, 1))
        B1 = -1/2
        B2 = 1/2
        for i, (node1, node2) in enumerate(self.conn):
            B1x1 = self.solution[node1] * B1 / self.lElem
            B2x2 = self.solution[node2] * B2 / self.lElem
            sField[i*n:(i+1)*n] = B1x1 + B2x2
        return sField
        

class TrussQuadratic(Truss):
    
    def __init__(self, E: float, L: float, A: float, nElem: int) -> None:
        super().__init__(E, L, A, nElem, nLocalNode = 3)
        self.k = self.E * self.A / self.lElem / 3
        self.localK = self.k * np.array(
            [[ 7, -8,  1],
             [-8, 16, -8],
             [1,  -8,  7]])

    def getLocalK(self, e: int):
        return self.localK

    def formFGlobal(self, qStart, qEnd, fStart, fEnd):
        FGlobal = np.zeros((self.nNode,1))
        qDiff = (qEnd - qStart) / self.L * self.lElem
        q = np.linspace(qStart, qEnd, self.nNode)
        
        for i, nodes in enumerate(self.conn):
            f0 = q[i] * self.lElem / 2
            x1, x2, x3 = self.x_node[nodes]
            FGlobal[nodes[0]] += f0 / 3 + qDiff / 30 * (4 * x1 +  2 * x2      -x3)
            FGlobal[nodes[1]] += 4 * f0 / 3 + qDiff / 30 * (2 * x1 + 16 * x2 + 2 * x3)
            FGlobal[nodes[2]] += f0 / 3 + qDiff / 30 * (   -x1 +  2 * x2 + 4 * x3)
        FGlobal[-1] += fEnd
        FGlobal[0] += fStart
        
        self.system.FGlobal = FGlobal

    def getXField(self, n: int = 20):
        N = n * self.nElem
        x = np.zeros((N, 1))
        ksi = np.linspace(-1, 1, n).reshape((n,1))
        N1 = ksi * (ksi - 1) / 2
        N2 = 1 - ksi*ksi
        N3 = ksi * (ksi + 1) / 2
        for i, (node1, node2, node3) in enumerate(self.conn):
            N1x1 = self.x_node[node1] * N1
            N2x2 = self.x_node[node2] * N2
            N3x3 = self.x_node[node3] * N3
            x[i*n:(i+1)*n] = N1x1 + N2x2 + N3x3
        return x

    def getDispField(self, n: int = 20):
        N = n * self.nElem
        uField = np.zeros((N, 1))
        ksi = np.linspace(-1, 1, n).reshape((n,1))
        N1 = ksi * (ksi - 1) / 2
        N2 = 1 - ksi*ksi
        N3 = ksi * (ksi + 1) / 2
        for i, (node1, node2, node3) in enumerate(self.conn):
            N1u1 = self.solution[node1] * N1
            N2u2 = self.solution[node2] * N2
            N3u3 = self.solution[node3] * N3
            uField[i*n:(i+1)*n] = N1u1 + N2u2 + N3u3
        return uField
    
    def getStrainField(self, n: int = 20):
        N = n * self.nElem
        eField = np.zeros((N, 1))
        ksi = np.linspace(-1, 1, n).reshape((n,1))
        B1 = ksi - 1/2
        B2 = - 2*ksi
        B3 = ksi + 1/2
        for i, (node1, node2, node3) in enumerate(self.conn):
            B1e1 = self.solution[node1] * B1 / self.lElem
            B2e2 = self.solution[node2] * B2 / self.lElem
            B3e3 = self.solution[node3] * B3 / self.lElem
            eField[i*n:(i+1)*n] = B1e1 + B2e2 + B3e3
        return eField