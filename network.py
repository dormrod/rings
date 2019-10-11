import numpy as np
from node import Node
from ring import Ring
from scipy.spatial import Delaunay


class Network:
    """
    Class to store nodes, edges and rings.
    """

    def __init__(self,node_crds,cell_dim):
        """
        Initialise with node coordinates, cell information.
        """

        # Set cell dimensions
        self.cell_dim = cell_dim

        # Set id logs
        self.next_node_id = 0
        self.next_ring_id = 0

        # Initialise nodes with given coordinates
        self.nodes = []
        for c in node_crds:
            n = Node(self.next_node_id)
            n.set_crd(c)
            self.nodes.append(n)
            self.next_node_id += 1

        # Initialise empty rings container
        self.rings = []


    def construct(self,edges):
        """
        Construct network with specified edges.
        """

        # Loop over edges and add node connections
        for edge in edges:
            # --- Little hack to remove periodicity ---#
            crd0 = self.nodes[edge[0]].get_crd()
            crd1 = self.nodes[edge[1]].get_crd()
            dx = np.abs(crd0[0]-crd1[0])
            dy = np.abs(crd0[1]-crd1[1])
            mic = (self.cell_dim[1]-self.cell_dim[0])/2.0
            if dx>mic or dy>mic:
                accept = False
            else:
                accept = True
            # --- Little hack over ---#
            if accept:
                self.nodes[edge[0]].add_node(edge[1])
                self.nodes[edge[1]].add_node(edge[0])

        # Deactivate any nodes which are 1-coordinate
        while True:
            clean = True
            for n in self.nodes:
                if n.active:
                    if n.num_nodes == 1:
                        clean = False
                        nid0 = n.id
                        nids = n.get_nodes()
                        for nid1 in nids:
                            self.nodes[nid1].del_node(nid0)
                        self.nodes[nid0].clear()
            if clean:
                break


    def triangulate(self):
        """
        Triangulate network using Delaunay.
        """

        # Triangulate and rings
        delaunay = Delaunay(self.get_node_crds())
        for simplex in delaunay.simplices:
            r = Ring(self.next_ring_id)
            for nid0 in simplex:
                r.add_node(nid0)
                self.nodes[nid0].add_ring(r.id)
            self.nodes[simplex[0]].add_node(simplex[1])
            self.nodes[simplex[0]].add_node(simplex[2])
            self.nodes[simplex[1]].add_node(simplex[0])
            self.nodes[simplex[1]].add_node(simplex[2])
            self.nodes[simplex[2]].add_node(simplex[0])
            self.nodes[simplex[2]].add_node(simplex[1])
            self.next_ring_id += 1
            self.rings.append(r)


    def get_node_crds(self):
        """
        Get node coordinates.
        """

        crds = []
        for n in self.nodes:
            crds.append(n.get_crd())
        return np.array(crds)


    def get_edges(self):
        """
        Get node ids that make up edges.
        """

        edges = []
        for n in self.nodes:
            if n.active:
                nid0 = n.id
                for nid1 in n.get_nodes():
                    if nid0<nid1:
                        edges.append([nid0,nid1])
        return np.array(edges,dtype=int)


    def get_rings(self):
        """
        Get node ids that make up rings.
        """

        rings = []
        for r in self.rings:
            if r.active:
                rings.append(np.array(r.get_nodes(),dtype=int))
        return rings