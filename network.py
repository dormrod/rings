
from collections import Counter
import numpy as np
from scipy.spatial import Delaunay

from node import Node
from ring import Ring
from plot_network import Plot



class Network:
    """
    Class to store nodes, edges and rings.
    """

    def __init__(self,
                 node_crds,
                 cell_dim,
                 periodic: bool = False):
        """
        Initialise with node coordinates, cell information.
        :param node_crds: an Nx2 numpy array of (x, y) positions.
        :param cell_dim: the dimensions of the containing unit cell.
        :param periodic: whether the box is periodic or not.
        """

        # Set cell dimensions
        self.periodic = periodic
        self.cell_dim = cell_dim
        self.pbc = np.array([self.cell_dim[0, 1] - self.cell_dim[0, 0],
                             self.cell_dim[1, 1] - self.cell_dim[1, 0]])
        # Minimum image convention length of half
        # the box size.
        self.mic = self.pbc / 2.0

        # Set id logs
        
        self.next_ring_id = 0

        # If the node coordinates are wrong, show now.
        assert len(node_crds.shape) == 2, "Node coordinates must be an"+\
                                          "Nx2 array. This array is of shape"+\
                                          f"{node_crds.shape}"
        assert node_crds.shape[1] == 2, "Node coordinates must be an"+\
                                        "Nx2 array. This array is of shape"+\
                                         f"{node_crds.shape}"
        # Initialise nodes with given coordinates
        self.nodes = []
        for next_node_id, crd in enumerate(node_crds):
            node = Node(next_node_id)
            node.set_crd(crd)
            self.nodes.append(node)

        # Initialise empty rings container
        self.rings = []

    def construct(self, edges):
        """
        Construct network with specified edges.
        :param edges: a list of (a, b) pairs with
        a and b being indicies of nodes within the network.
        """

        # Loop over edges and add node connections
        if self.periodic:
            for edge in edges:
                self.nodes[edge[0]].add_node(edge[1])
                self.nodes[edge[1]].add_node(edge[0])
        else:
            # If aperiodic explicitly remove any connections
            # greater than the minimum image convention.
            for edge in edges:
                crd0 = self.nodes[edge[0]].get_crd()
                crd1 = self.nodes[edge[1]].get_crd()
                dx = np.abs(crd0[0] - crd1[0])
                dy = np.abs(crd0[1] - crd1[1])
                if not (dx > self.mic[0] or dy > self.mic[1]):
                    self.nodes[edge[0]].add_node(edge[1])
                    self.nodes[edge[1]].add_node(edge[0])

        # Deactivate any nodes which are 0/1-coordinate
        while True:
            clean = True
            for n in self.nodes:
                if n.active:
                    if n.num_nodes == 0:
                        clean = False
                        self.nodes[n.id].clear()
                    elif n.num_nodes == 1:
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
        Create a Delauney triangulation of the network
        using scipy.spatial.Delauney
        :return None: but modifies the member variables.
        """

        # Get active node coordinates
        active_crds = []
        active_map = []
        for i, n in enumerate(self.nodes):
            if n.active:
                active_map.append(i)
                active_crds.append(n.get_crd())
        active_crds = np.array(active_crds)
        num_active = active_crds.shape[0]

        # Make periodic images if required
        if self.periodic:
            central_crds = np.copy(active_crds)
            active_crds = np.zeros((9*num_active, 2))
            active_crds[:num_active, :] = central_crds
            counter = num_active
            for y in [-1, 0, 1]:
                for x in [-1, 0, 1]:
                    if y != 0 or x != 0:
                        next_counter = counter+num_active
                        active_crds[counter:next_counter, 0] = central_crds[: ,0] + x * self.pbc[0]
                        active_crds[counter:next_counter, 1] = central_crds[:, 1] + y * self.pbc[1]
                        counter += num_active

        # Triangulate active nodes and find rings
        delaunay = Delaunay(active_crds)
        added_simplices = set()
        for simplex in delaunay.simplices:
            if np.any(simplex < num_active):
                simplex = np.sort(simplex % num_active).astype(int)
                nids = [active_map[simplex[i]] for i in range(3)]
                # Check which simplices we've seen before, and skip
                # if we've seen it.
                simplex_id = tuple(sorted(nids))
                if simplex_id in added_simplices:
                    continue

                r = Ring(self.next_ring_id)
                for nid in nids:
                    r.add_node(nid)
                    self.nodes[nid].add_ring(r.id)
                    for other_nid in nids:
                        # Connect the other node ids to this,
                        # but don't connect to ourselves.
                        if nid == other_nid:
                            continue
                        self.nodes[nid].add_node(other_nid)
                self.next_ring_id += 1
                self.rings.append(r)
                added_simplices.add(simplex_id)
        
        for n in self.nodes:
            n.unique_nodes()

        # Add perimeter ring
        if not self.periodic:
            self.find_perimeter()

    def find_perimeter(self):
        """
        Find infinite ring corresponding to face.
        """

        # Count the number of times each edge comes up in
        # the list of rings.
        edge_count = Counter([tuple(sorted(e)) for r in self.rings
                              for e in r.get_edges()
                              if r.active])

        # Collect edges used only once
        perimeter_edges = [edge for edge, count in edge_count.items()
                           if count == 1]

        # Trace path around perimeter
        edge = perimeter_edges.pop(0)
        perimeter = [edge[1]]
        while len(perimeter_edges) > 0:
            for i, edge in enumerate(perimeter_edges):
                if perimeter[-1] in edge:
                    edge = perimeter_edges.pop(i)
                    if perimeter[-1] == edge[0]:
                        perimeter.append(edge[1])
                    else:
                        perimeter.append(edge[0])

        # Add perimeter ring to the class members.
        r = Ring(self.next_ring_id, True)
        for nid in perimeter:
            r.add_node(nid)
            self.nodes[nid].add_ring(r.id)
        self.rings.append(r)
        self.next_ring_id += 1

    def map(self, other):
        """
        Find differences between networks and map other onto self.
        """

        # Find edges in both networks
        self_edges = [(e[0], e[1]) for e in self.get_edges()]
        other_edges = [(e[0], e[1]) for e in other.get_edges()]

        # Find edges unique to both networks
        unique_self_edges = []
        for e in self_edges:
            if e not in other_edges and e not in unique_self_edges:
                unique_self_edges.append(e)

        unique_other_edges = []
        for e in other_edges:
            if e not in self_edges and e not in unique_other_edges:
                unique_other_edges.append(e)
        if len(unique_other_edges) > 0: # Should be empty
            return False
        delete_edges = unique_self_edges

        # Delete required edges and merge rings
        for edge in delete_edges:
            nid0, nid1 = edge
            rids = list(set(self.nodes[nid0].get_rings()).intersection(self.nodes[nid1].get_rings()))
            if len(rids) > 2:
                red_rids = []
                for rid in rids:
                    if self.rings[rid].check_edge(nid0, nid1):
                        red_rids.append(rid)
                rids = red_rids
            if len(rids) != 2:
                print(nid0, nid1)
                print(self.nodes[nid0].rings, self.nodes[nid1].rings)
                print(rids)
                plot = Plot(nodes=True, rings=True)
                plot(self)
                return
            r = self.rings[rids[0]].merge(self.rings[rids[1]], 
                                          nid0,
                                          nid1,
                                          self.next_ring_id)
            self.nodes[nid0].del_node(nid1)
            self.nodes[nid1].del_node(nid0)
            for nid2 in r.get_nodes():
                self.nodes[nid2].swap_ring(rids[0], r.id)
                self.nodes[nid2].swap_ring(rids[1], r.id)
                self.rings[rids[0]].clear()
                self.rings[rids[1]].clear()
            self.rings.append(r)
            self.next_ring_id += 1

    def get_node_crds(self):
        """
        Returns an Nx2 numpy array of node coordinates,
        if they have been changed by the process.
        :return coords: Nx2 array of coordinates, in order of node ids.
        """
        return np.array([n.get_crd() for n in self.nodes])

    def get_edges(self,unique=True):
        """
        Get node ids that make up edges.
        :param unique: Whether to double count [a, b] and [b, a].
        If false, returns an array with only unique edges, which
        have a < b.
        :return edges: an Nx2 array of edge connections.
        """

        edges = []
        for n in self.nodes:
            if n.active:
                nid0 = n.id
                for nid1 in n.get_nodes():
                    edge_code = [nid0, nid1]
                    if unique:
                        edge_code = sorted(edge_code)
                    edges.append(edge_code)
        return np.array(edges, dtype=int)

    def get_rings(self, infinite=False):
        """
        Get node ids that make up rings.
        :param infinite: returns the node ids of the "infinite"
        perimiter ring.
        """

        rings = []
        for r in self.rings:
            if r.active and r.infinite == infinite:
                rings.append(np.array(r.get_nodes(), dtype=int))
        return rings
