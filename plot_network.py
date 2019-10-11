import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import matplotlib.pylab as pylab
import numpy as np
import sys

class Plot:
    """
    Plot procrystalline lattice.
    """


    def __init__(self,nodes=False,cnxs=False,rings=False,periodic=False):
        """
        Initialise with plot options.
        """

        # Get options
        self.nodes = nodes
        self.cnxs = cnxs
        self.rings = rings
        self.periodic = periodic

        # Set up empty figure
        params = {"figure.figsize": (6, 6)}
        pylab.rcParams.update(params)
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax.set_axis_off()


    def __call__(self,network,**kwargs):
        """
        Plot procrystalline lattice.
        """

        # Images to generate based on periodicity
        if self.periodic:
            images = [-1,0,1]
        else:
            images = [0]

        # Unpack options
        self.lw = kwargs.get("lw",1.0)
        self.lc = kwargs.get("lc","k")
        self.ms = kwargs.get("ms",10)
        self.mc = kwargs.get("mc","k")
        save = kwargs.get("save",False)

        # Load data
        self.load_network(network)

        # Add images to plot
        for y in images:
            for x in images:
                self.plot_rings(x,y)
                self.plot_nodes(x,y)
                self.plot_cnxs(x,y)

        # Display
        if save:
            plt.savefig("config.png",dpi=400)
        plt.show()


    def load_network(self,network):
        """
        Load coordinates and network information.
        """

        # Node coordination, periodicity and node crds
        self.pbc = network.cell_dim
        self.mic = self.pbc/2
        self.node_crds = network.get_node_crds()
        self.mean_ring_size = 6 # self.node_cnxs = []

        # Make connections accounting for periodicity
        self.node_cnxs = network.get_edges()
        self.node_cnx_crds = np.zeros((self.node_cnxs.shape[0],4))
        crd_i = np.zeros(2)
        crd_j = np.zeros(2)
        for i,p in enumerate(self.node_cnxs):
            crd_i[:] = self.node_crds[p[0]][:]
            crd_j[:] = self.node_crds[p[1]][:]
            x = crd_j[0]-crd_i[0]
            y = crd_j[1]-crd_i[1]
            if x>self.mic[0]: x-=self.pbc[0]
            elif x<-self.mic[0]: x+=self.pbc[0]
            if y>self.mic[1]: y-=self.pbc[1]
            elif y<-self.mic[1]: y+=self.pbc[1]
            self.node_cnx_crds[i,0] = crd_i[0]
            self.node_cnx_crds[i,1] = crd_i[0]+x
            self.node_cnx_crds[i,2] = crd_i[1]
            self.node_cnx_crds[i,3] = crd_i[1]+y

        # Make rings accounting for periodicity
        self.node_rings = network.get_rings()
        self.ring_crds = []
        self.max_ring_size = 0
        for ring in self.node_rings:
            crds = np.zeros((ring.size,2))
            for i,j in enumerate(ring):
                crds[i,:] = self.node_crds[j,:]
            for i in range(1,ring.size):
                x = crds[i,0] - crds[i-1,0]
                y = crds[i,1] - crds[i-1,1]
                if x>self.mic[0]: x -= self.pbc[0]
                elif x<-self.mic[0]: x += self.pbc[0]
                if y>self.mic[1]: y -= self.pbc[1]
                elif y<-self.mic[1]: y += self.pbc[1]
                crds[i,0] = crds[i-1,0] + x
                crds[i,1] = crds[i-1,1] + y
            self.ring_crds.append(crds)
            if ring.size > self.max_ring_size:
                self.max_ring_size = ring.size
        self.init_ring_colours(self.mean_ring_size,self.max_ring_size)


    def plot_nodes(self,x_shift,y_shift):
        """
        Plot image nodes as scatter. 
        """

        if not self.nodes: return # Bounce if option not selected

        self.ax.scatter(self.node_crds[:,0]+x_shift*self.pbc[0],self.node_crds[:,1]+y_shift*self.pbc[1],
                        marker="o",s=self.ms,c=self.mc,zorder=1)

        # for i,c in enumerate(self.node_crds):
        #     self.ax.text(c[0],c[1],i,size=8)


    def plot_cnxs(self,x_shift,y_shift):
        """
        Plot image connections as lines. 
        """

        if not self.cnxs: return # Bounce if option not selected

        self.node_cnx_crds[:,:2] += x_shift*self.pbc[0]
        self.node_cnx_crds[:,2:] += y_shift*self.pbc[1]
        for cnx_crd in self.node_cnx_crds:
            self.ax.plot(cnx_crd[:2],cnx_crd[2:],c=self.lc,lw=self.lw)
        self.node_cnx_crds[:,:2] -= x_shift*self.pbc[0]
        self.node_cnx_crds[:,2:] -= y_shift*self.pbc[1]


    def plot_rings(self,x_shift,y_shift):
        """
        Plot rings as polygons.
        """

        if not self.rings: return # Bounce if option not selected

        patches = []
        colours = []
        for ring in self.ring_crds:
            ring[:,0] += x_shift*self.pbc[0]
            ring[:,1] += y_shift*self.pbc[1]
            patches.append(Polygon(np.array(ring), True))
            colours.append(self.ring_colours[ring[:,0].size])
            ring[:,0]-=x_shift*self.pbc[0]
            ring[:,1]-=y_shift*self.pbc[1]
        self.ax.add_collection(PatchCollection(patches,facecolor=colours,linewidths=self.lw,edgecolor="k",zorder=0))


    def init_ring_colours(self,av_ring_size=6,max_ring_size=10):
        """
        Initialise colouring for rings.
        """

        map_lower = cm.get_cmap('Blues_r', 128)
        map_upper = cm.get_cmap('Reds', 128)
        map_mean=cm.get_cmap("Greys")
        map_lower=ListedColormap(map_lower(np.arange(20,100)))
        map_upper=ListedColormap(map_upper(np.arange(20,100)))

        norm_lower=Normalize(vmin=av_ring_size-3,vmax=av_ring_size)
        norm_upper=Normalize(vmin=av_ring_size,vmax=av_ring_size+6)
        colour_mean=map_mean(50)
        self.ring_colours=[]
        for i in range(max_ring_size+1):
            if i < 3:
                self.ring_colours.append("white")
            elif np.abs(i-av_ring_size)<1e-6:
                self.ring_colours.append(colour_mean)
            elif i<av_ring_size:
                self.ring_colours.append(map_lower(norm_lower(i)))
            else:
                self.ring_colours.append(map_upper(norm_upper(i)))


if __name__ == "__main__":

    prefix = sys.argv[1]
    sample = int(sys.argv[2])

    if len(sys.argv) <= 3:
        nodes = True
        cnxs = False
        rings = True
        periodic = False
        save = False
    else:
        flags = sys.argv[3]
        nodes = 'n' in flags
        cnxs = 'c' in flags
        rings = 'r' in flags
        periodic = 'p' in flags
        save = 's' in flags

    plot=Plot(nodes=nodes,cnxs=cnxs,rings=rings,periodic=periodic)
    plot(prefix,sample,save=save)