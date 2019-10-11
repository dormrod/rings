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


    def add_ring(self,rid):
        """
        Add ring id.
        """

        self.rings.append(rid)
        self.num_rings += 1


    def get_nodes(self):
        """
        Get node ids. 
        """

        return self.nodes


    # def get_edges(self):
    #     """
    #     Get edge ids.
    #     """
    #
    #     return self.edges


    def clear(self):
        """
        Remove all associated objects and deactivate.
        """

        self.nodes = []
        self.rings = []
        self.num_nodes = 0
        self.num_rings = 0
        self.active = False

