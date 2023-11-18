import numpy as np

class System:
    def __init__(self, nNode: int) -> None:
        self.nNode = nNode
        self.nDOF = None
        self.nKnowns = None
        self.KGlobal = None
        self.FGlobal = None
        self.K_kk =   None
        self.K_uk =   None
        self.K_uu =   None
        self.F_u  =   None
        self.D_k  =   None
        self.node_u = None
        self.node_k = None

    def solve(self):
        D_u = np.linalg.solve(self.K_uu, self.F_u - self.K_uk @ self.D_k )
        if self.nKnowns == 1:
            return np.vstack((self.D_k[0], D_u))
        if self.nKnowns == 2:
            return np.vstack((self.D_k[0], D_u, self.D_k[-1]))
        return 

    def form(self, D_k: np.ndarray, u0Known = True, uLKnown = False):
        self.nKnowns = int(u0Known) + int(uLKnown)
        self.nDOF = self.nNode - self.nKnowns
        self.D_k = D_k
        self.K_kk = np.zeros((self.nKnowns, self.nKnowns))
        self.K_uk = np.zeros((self.nDOF, self.nKnowns))
        self.K_uu = np.zeros((self.nDOF, self.nDOF))
        self.F_u  = np.zeros((self.nDOF, 1))
        self._formMappings(u0Known, uLKnown)
        self._formMatrices()

    def _formMappings(self, u0Known = True, uLKnown = False):
        # Mapping from Global to Unknowns
        self.node_u = np.arange(self.nNode, dtype=np.int32)
        self.node_k = np.full((self.nNode,), -1, dtype=np.int32)
        # Mark knowns in the unknowns list
        if u0Known:
            self.node_k[0] = 0
            self.node_u -= 1
        if uLKnown:
            self.node_k[-1] = 1
            self.node_u[-1] = -1

    def _formMatrices(self):
        if self.KGlobal is None:
            raise ValueError('KGlobal is undefined.')
        if self.FGlobal is None:
            raise ValueError('FGlobal is undefined.')
        for i in range(self.nNode):
            for j in range(self.nNode):
                if self.node_u[i] != -1:
                    if self.node_u[j] != -1:
                        self.K_uu[self.node_u[i], self.node_u[j]] = self.KGlobal[i,j]
                    
                if self.node_k[i] != -1:
                    if self.node_k[j] != -1:
                        self.K_kk[self.node_k[i], self.node_k[j]] = self.KGlobal[i,j]
                    if self.node_u[j] != -1:
                        self.K_uk[self.node_u[j], self.node_k[i]] = self.KGlobal[i,j]
            if self.node_u[i] != -1:
                self.F_u[self.node_u[i]] = self.FGlobal[i]