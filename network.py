import numpy as np
from node import Node
from ring import Ring
from scipy.spatial import Delaunay
from plot_network import Plot
import matplotlib.pyplot as plt

class Network:
    """
    Class to store nodes, edges and rings.
    """

    def __init__(self,node_crds,cell_dim,periodic=False):
        """
        Initialise with node coordinates, cell information.
        """

        # Set cell dimensions
        self.periodic = periodic
        self.cell_dim = cell_dim
        self.pbc =  np.array([self.cell_dim[0,1]-self.cell_dim[0,0],self.cell_dim[1,1]-self.cell_dim[1,0]])
        self.mic = self.pbc / 2.0

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
        if self.periodic:
            for edge in edges:
                self.nodes[edge[0]].add_node(edge[1])
                self.nodes[edge[1]].add_node(edge[0])
        else:
            # If aperiodic explicitly remove any connections greater than wrapping
            for edge in edges:
                crd0 = self.nodes[edge[0]].get_crd()
                crd1 = self.nodes[edge[1]].get_crd()
                dx = np.abs(crd0[0]-crd1[0])
                dy = np.abs(crd0[1]-crd1[1])
                if dx>self.mic[0] or dy>self.mic[1]:
                    accept = False
                else:
                    accept = True
                if accept:
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
        Triangulate network using Delaunay.
        """

        # Get active node coordinates
        active_crds = []
        active_map = []
        for i,n in enumerate(self.nodes):
            if n.active:
                active_map.append(i)
                active_crds.append(n.get_crd())
        active_crds = np.array(active_crds)
        num_active = active_crds.shape[0]

        # Make periodic images if required
        if self.periodic:
            central_crds = np.copy(active_crds)
            active_crds = np.zeros((9*num_active,2))
            active_crds[:num_active,:] = central_crds
            counter = num_active
            for i,y in enumerate([-1,0,1]):
                for j,x in enumerate([-1,0,1]):
                    if y!=0 or x!=0:
                        active_crds[counter:counter+num_active,0] = central_crds[:,0]+x*self.pbc[0]
                        active_crds[counter:counter+num_active,1] = central_crds[:,1]+y*self.pbc[1]
                        counter += num_active
        # self.fig = plt.figure()
        # self.ax = self.fig.add_subplot(111)
        # self.ax.scatter(active_crds[:,0], active_crds[:,1],marker="o",s=2,c='k',zorder=1)
        # for i,c in enumerate(active_crds):
        #     self.ax.text(c[0],c[1],i,size=6)


        # Triangulate active nodes and find rings
        delaunay = Delaunay(active_crds)
        if self.periodic:
            added_simplices = []
            for simplex in delaunay.simplices:
                if np.any(simplex<num_active):
                    nid3 = simplex[0]
                    nid4 = simplex[1]
                    nid5 = simplex[2]
                    simplex = np.sort(simplex%num_active).astype(int)
                    nid0 = active_map[simplex[0]]
                    nid1 = active_map[simplex[1]]
                    nid2 = active_map[simplex[2]]
                    simplex_id = "#{}#{}#{}".format(nid0,nid1,nid2)
                    if simplex_id not in added_simplices:
                        # self.ax.plot([active_crds[nid3,0],active_crds[nid4,0]],[active_crds[nid3,1],active_crds[nid4,1]])
                        # self.ax.plot([active_crds[nid3,0],active_crds[nid5,0]],[active_crds[nid3,1],active_crds[nid5,1]])
                        # self.ax.plot([active_crds[nid5,0],active_crds[nid4,0]],[active_crds[nid5,1],active_crds[nid4,1]])
                        r = Ring(self.next_ring_id)
                        r.add_node(nid0)
                        r.add_node(nid1)
                        r.add_node(nid2)
                        self.nodes[nid0].add_ring(r.id)
                        self.nodes[nid1].add_ring(r.id)
                        self.nodes[nid2].add_ring(r.id)
                        self.nodes[nid0].add_node(nid1)
                        self.nodes[nid0].add_node(nid2)
                        self.nodes[nid1].add_node(nid0)
                        self.nodes[nid1].add_node(nid2)
                        self.nodes[nid2].add_node(nid0)
                        self.nodes[nid2].add_node(nid1)
                        self.next_ring_id += 1
                        self.rings.append(r)
                        added_simplices.append(simplex_id)
        else:
            for simplex in delaunay.simplices:
                nid0 = active_map[simplex[0]]
                nid1 = active_map[simplex[1]]
                nid2 = active_map[simplex[2]]
                r = Ring(self.next_ring_id)
                r.add_node(nid0)
                r.add_node(nid1)
                r.add_node(nid2)
                self.nodes[nid0].add_ring(r.id)
                self.nodes[nid1].add_ring(r.id)
                self.nodes[nid2].add_ring(r.id)
                self.nodes[nid0].add_node(nid1)
                self.nodes[nid0].add_node(nid2)
                self.nodes[nid1].add_node(nid0)
                self.nodes[nid1].add_node(nid2)
                self.nodes[nid2].add_node(nid0)
                self.nodes[nid2].add_node(nid1)
                self.next_ring_id += 1
                self.rings.append(r)
        for n in self.nodes:
            n.unique_nodes()

        # Add perimeter ring
        if not self.periodic:
            self.find_perimeter()


    def find_perimeter(self):
        """
        Find infinite ring corresponding to face.
        """

        # Brute force get unique edges
        edge_count = {}
        for n in self.nodes:
            if n.active:
                nid0 = n.id
                for nid1 in n.get_nodes():
                    if nid0<nid1:
                        edge_code = (nid0,nid1)
                        edge_count[edge_code] = 0

        # Count number of times edges used in rings
        for r in self.rings:
            if r.active:
                for e in r.get_edges():
                    nid0,nid1 = e
                    if nid0<nid1:
                        edge_code = (nid0,nid1)
                    else:
                        edge_code = (nid1,nid0)
                    edge_count[edge_code] += 1

        # Collect edges used only once
        perimeter_edges = []
        for edge,count in edge_count.items():
            if count == 1:
                perimeter_edges.append(edge)

        # Trace path around perimeter
        edge = perimeter_edges.pop(0)
        perimeter = [edge[1]]
        while len(perimeter_edges)>0:
            for i,pedge in enumerate(perimeter_edges):
                if perimeter[-1] in pedge:
                    edge = perimeter_edges.pop(i)
                    if perimeter[-1] == edge[0]:
                        perimeter.append(edge[1])
                    else:
                        perimeter.append(edge[0])

        # Add perimeter ring
        r = Ring(self.next_ring_id,True)
        for nid in perimeter:
            r.add_node(nid)
            self.nodes[nid].add_ring(r.id)
        self.rings.append(r)
        self.next_ring_id += 1


    def map(self,other):
        """
        Find differences between networks and map other onto self.
        """

        # Find edges in both networks
        self_edges = [(e[0],e[1]) for e in self.get_edges()]
        other_edges = [(e[0],e[1]) for e in other.get_edges()]

        # Find edges unique to both networks
        unique_self_edges = []
        for e in self_edges:
            if e not in other_edges and e not in unique_self_edges:
                unique_self_edges.append(e)

        unique_other_edges = []
        for e in other_edges:
            if e not in self_edges and e not in unique_other_edges:
                unique_other_edges.append(e)
        if len(unique_other_edges)>0: # Should be empty
            return False
        delete_edges = unique_self_edges

        # Delete required edges and merge rings
        for edge in delete_edges:
            nid0,nid1 = edge
            rids=list(set(self.nodes[nid0].get_rings()).intersection(self.nodes[nid1].get_rings()))
            if (len(rids)>2):
                red_rids = []
                for rid in rids:
                    if self.rings[rid].check_edge(nid0,nid1):
                        red_rids.append(rid)
                rids = red_rids
            if len(rids)!=2:
                print(nid0,nid1)
                print(self.nodes[nid0].rings,self.nodes[nid1].rings)
                print(rids)
                plot=Plot(nodes=True,rings=True)
                plot(self)
                return
            r = self.rings[rids[0]].merge(self.rings[rids[1]],nid0,nid1,self.next_ring_id)
            self.nodes[nid0].del_node(nid1)
            self.nodes[nid1].del_node(nid0)
            for nid2 in r.get_nodes():
                self.nodes[nid2].swap_ring(rids[0],r.id)
                self.nodes[nid2].swap_ring(rids[1],r.id)
                self.rings[rids[0]].clear()
                self.rings[rids[1]].clear()
            self.rings.append(r)
            self.next_ring_id += 1
            # plot=Plot(nodes=True,rings=True)
            # plot(self,ms=20,save=False)


    def get_node_crds(self):
        """
        Get node coordinates.
        """

        crds = []
        for n in self.nodes:
            crds.append(n.get_crd())
        return np.array(crds)


    def get_edges(self,unique=True):
        """
        Get node ids that make up edges.
        """

        edges = []
        if unique:
            for n in self.nodes:
                if n.active:
                    nid0 = n.id
                    for nid1 in n.get_nodes():
                        if nid0<nid1:
                            edges.append([nid0,nid1])
        else:
            for n in self.nodes:
                if n.active:
                    nid0 = n.id
                    for nid1 in n.get_nodes():
                        edges.append([nid0,nid1])

        return np.array(edges,dtype=int)


    def get_rings(self,infinite=False):
        """
        Get node ids that make up rings.
        """

        rings = []
        for r in self.rings:
            if r.active and r.infinite==infinite:
                rings.append(np.array(r.get_nodes(),dtype=int))
        return rings


