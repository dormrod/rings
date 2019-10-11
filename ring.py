import numpy as jnp
from network_container import NetworkContainer


class Ring(NetworkContainer):
    """
    Ring in network.
    """

    def __init__(self,id,infinite=False):
        """
        Initialise with flag for infinite ring.
        """

        self.infinite = infinite
        super(Ring, self).__init__(id)


    def merge(self,other,nid0,nid1,id):
        """
        Generate ring by merging with other, keeping continuous path of nodes.
        """

        # Check whether infinite
        if self.infinite or other.infinite:
            infinite = True
        else:
            infinite = False

        # Initialise ring
        ring = Ring(id,infinite)

        # Trace path around nodes (fiddly)
        for i,n in enumerate(self.nodes):
            if n == nid0: posa0 = i
            elif n == nid1: posa1 = i
        for i,n in enumerate(other.nodes):
            if n == nid0: posb0 = i
            elif n == nid1: posb1 = i
        if abs(posa1-posa0)==1:
            if posa1<posa0: a_forwards = True
            else: a_forwards = False
        else:
            if posa0==0: a_forwards = True
            else: a_forwards = False
        if abs(posb1-posb0) == 1:
            if posb1-posb0<0: b_forwards = False
            else: b_forwards = True
        else:
            if posb0==0: b_forwards = False
            else: b_forwards = True
        if a_forwards:
            for i in range(self.num_nodes-1):
                ring.add_node(self.nodes[(posa0+i)%self.num_nodes])
        else:
            for i in range(self.num_nodes-1):
                ring.add_node(self.nodes[(posa0+self.num_nodes-i)%self.num_nodes])
        if b_forwards:
            for i in range(other.num_nodes-1):
                ring.add_node(other.nodes[(posb1+i)%other.num_nodes])
        else:
            for i in range(other.num_nodes-1):
                ring.add_node(other.nodes[(posb1+other.num_nodes-i)%other.num_nodes])

        return ring

