import numpy as np
from network_container import NetworkContainer


class Node(NetworkContainer):
    """
    Node in network.
    """

    def set_crd(self,crd):
        """
        Set x,y coordinate.
        """

        self.crd = np.zeros(2)
        self.crd[:] = crd[:]


    def get_crd(self):
        """
        Get x,y coordinate.
        """

        return self.crd

