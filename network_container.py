import numpy as np


class NetworkContainer:
    """
    Base class for node and ring objects.
    """

    def __init__(self,id):
        """
        Initialise with id.
        """

        # Setup containers
        self.nodes = [] # Associated nodes
        self.rings = [] # Associated rings

        # Setup book-keeping variables
        self.id = id
        self.active = True
        self.num_nodes = 0
        self.num_rings = 0


    def add_node(self,nid):
        """
        Add node id.
        """

        self.nodes.append(nid)
        self.num_nodes += 1


    def del_node(self,nid):
        """
        Delete node id.
        """

        self.nodes = [n for n in self.nodes if n != nid]
        self.num_nodes = len(self.nodes)


    def unique_nodes(self):
        """
        Make nodes unique
        """

        self.nodes = list(set(self.nodes))
        self.num_nodes = len(self.nodes)


    def add_ring(self,rid):
        """
        Add ring id.
        """

        self.rings.append(rid)
        self.num_rings += 1


    def swap_ring(self,delId,addId):
        """
        Swap instances of ring id for another.
        """

        self.rings = list(set([i if i != delId else addId for i in self.rings]))
        self.num_rings = len(self.rings)


    def get_nodes(self):
        """
        Get node ids. 
        """

        return self.nodes


    def get_edges(self):
        """
        Get node pair ids forming edges.
        """

        edges = [[self.nodes[i],self.nodes[(i+1)%self.num_nodes]] for i in range(self.num_nodes)]
        return edges


    def check_edge(self,nid0,nid1):
        """
        Check node ids adajacent.
        """

        pos0 = -1
        for i,n in enumerate(self.nodes):
            if n==nid0:
                pos0=i
                break

        if pos0==-1: return False
        if self.nodes[(pos0+1)%self.num_nodes]==nid1: return True
        elif self.nodes[(pos0+self.num_nodes-1)%self.num_nodes]==nid1: return True
        else: return False


    def get_rings(self):
        """
        Get ids of rings
        """

        return self.rings


    def clear(self):
        """
        Remove all associated objects and deactivate.
        """

        self.nodes = []
        self.rings = []
        self.num_nodes = 0
        self.num_rings = 0
        self.active = False


    def clear_rings(self):
        """
        Remove all associated rings.
        """

        self.rings = []
        self.num_rings = 0

