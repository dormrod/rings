import numpy as np
from network import Network
from plot_network import Plot

if __name__ == "__main__":

    crds = np.genfromtxt("./coords.dat",delimiter=",",usecols=(1,2)).astype(float)
    cnxs = np.genfromtxt("./edges.dat",delimiter=",").astype(int)

    target_network = Network(crds,np.array([0,100]))
    target_network.construct(cnxs)

    tri_network = Network(crds,np.array([0,100]))
    tri_network.construct(cnxs)
    tri_network.triangulate()

    plot_tri = Plot(nodes=True,cnxs=False,rings=True)
    plot_tri(tri_network)

    plot_target = Plot(nodes=True,cnxs=True,rings=False)
    plot_target(target_network)

    tri_network.map(target_network)

    plot_map = Plot(nodes=True,cnxs=False,rings=True)
    plot_map(tri_network)
