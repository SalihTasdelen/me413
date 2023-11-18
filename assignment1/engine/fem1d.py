import numpy as np
from .system import System

class FEM1D:
    def __init__(self, L:float, nElem: int, nLocalNode:int) -> None:
        """Creates a FEM 1D Element that manipulates the system object at
           the background.

        Args:
            L (float): Length of the Truss Element.
            nElem (int): Number of Elements.
            nLocalNode (int): Number of local nodes inside an element.
        """        
        self.nElem = nElem
        self.nLocalNode = nLocalNode
        self.nNode = (nElem - 1) * (nLocalNode - 1) + nLocalNode
        self.L = L
        self.lElem = self.L / self.nElem
        self.conn = self._connect(nLocalNode, nElem)
        self.x_node = np.linspace(0, self.L, self.nNode)
        self.system = System(self.nNode)
        self.solution = None
    
    def solve(self):
        self.solution = self.system.solve()
        return self.solution

    def formKGlobal(self):
        """Forms the 1D K Global Matrix

        Raises:
            NotImplementedError: Base FEM1D Class leaves KGlobal definitoin undefined.
        """        
        raise NotImplementedError
    
    def formFGlobal(self):
        """Forms the 1D F Global Matrix

        Raises:
            NotImplementedError: Base FEM1D Class leaves FGlobal definition undefined.
        """        
        raise NotImplementedError
    
    def formDirichlet(self, u0: float, uL: float):
        """Applies the following boundary conditions;
                u(0) = u0, u(L) = uL

        Args:
            u0 (float): u(0) = u0
            uL (float): u(L) = uL
        """
        D_k = np.array([u0, uL]).reshape((2, 1))
        self.system.form(D_k, u0Known=True, uLKnown=True)

    def formDirichletNeumann(self):
        """Applies the following boundary conditions;
            u(0) = 0, du/dx(L) = 0
        """
        D_k = np.zeros((1,1))
        self.system.form(D_k, u0Known=True, uLKnown=False)

    @staticmethod
    def _connect(n: int = 2, N: int = 1):
        """Generates connectivity information in 1D.

        Args:
            n (int, optional): Number of nodes in an element. Defaults to 2.
            N (int, optional): Number of elements. Defaults to 1.

        Returns:
            np.ndarray: Array of elements with their corresponding node indices.
        """
        conn = np.zeros((N, n), dtype=np.int32)
        elementIdx = np.arange(N*n)
        for i in range(N):
            idx = i*(n - 1)
            conn[i] = elementIdx[idx: idx + n]
        return conn