#%%mini program to generate a 3D grid of nodes
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation



### other
#plt.clf()

## gleiche value
# total_nodes_3d        = [216]#[512]#[1331,729,512]#
# value3                = [2.8,2.8,2.8,2.8,2.8]
value3                = [1.8,1.8,1.8,1.8,1.8,1.8]
# gleiche node unterschiedlich value
# total_nodes_3d        = [729,729,729,729]#
# value3                = [1.3,1.5,1.8,2.8]#
# gleiche degree  range
total_nodes_3d        = [64,216,512,729,1331,2744]#64,216,
skin_list3d=[]
av_clust3             = []
nom_value3            = []
density_3d            = []
av_path3d_list        = []
max_degree3           = []
min_degree3           = []
av_degree3            = []
list_degree_sequence3 = []
list_clustering_seq3  = []
list_clust_bins3d     = []
path_len_list         = []
list_number_edges_3d  = np.zeros(len(total_nodes_3d))
list_number_nodes_3d  = np.zeros(len(total_nodes_3d))
diameter3             = {}
freq_percent3         = {}
#%% loop over node number list
for pp in range(len(total_nodes_3d)):
    bedingung3             = np.zeros((total_nodes_3d[pp]))
    distances3=[]
#%% meshgrid setup position
    scale_x, scale_y, scale_z = 1.0, 1.0, 1.0
    xn, yn, zn = round(np.cbrt(total_nodes_3d[pp]),0),round(np.cbrt(total_nodes_3d[pp]),0),round(np.cbrt(total_nodes_3d[pp]),0) ##number of nodes in x,y ,z directions
    länge =   round(np.cbrt(total_nodes_3d[pp]),0)-1 #xn-1 #                           # länge des würfels
    x = scale_x*np.linspace(0, länge, int(xn), endpoint=True)                  ##x-coordinates
    y = scale_y*np.linspace(0, länge, int(yn), endpoint=True)                  ##y-coordinates
    z = scale_z*np.linspace(0, länge, int(zn), endpoint=True)                  ##z-coordinates
    xv, yv ,zv = np.meshgrid(x, y, z) 
    xvr = np.reshape(xv,xv.size)
    yvr = np.reshape(yv,yv.size)
    zvr = np.reshape(zv,zv.size)
    random1=np.random.normal(loc=0,scale=(länge/int(xn))/2 ,size=total_nodes_3d[pp])    #=(länge/int(xn))/2      # random variations of distance
    random2=np.random.normal(loc=0,scale=(länge/int(xn))/2 ,size=total_nodes_3d[pp])          # random variations of distance
    random3=np.random.normal(loc=0,scale=(länge/int(xn))/2 ,size=total_nodes_3d[pp])          # random variations of distance
    pos3={}
    pos = {}                                                                    # Build position dictionary from x, y, z coordinates
    for i, (x_val, y_val, z_val) in enumerate(zip(xvr, yvr, zvr)):
    # not random pertubation
        pos[i] = [x_val,y_val,z_val]
    # # # random pertubation
        # pos[i] = [x_val+random1[i], y_val+random2[i], z_val+random3[i]]
        # pos3[i] = [x_val+random1[i], y_val+random2[i]]
        # xvr[i] += random1[i]
        # yvr[i] += random2[i]
        # zvr[i] += random3[i]
    xv = np.reshape(xvr,xv.shape)
    yv = np.reshape(yvr,yv.shape)
    zv = np.reshape(zvr,zv.shape)
#%%
    G3 = nx.path_graph(xvr.size)
    G3.remove_edges_from(list(G3.edges))
    edge_list=[]
    for node_ref in G3.nodes(): #loop over each node
        x_node = pos[node_ref][0]
        y_node = pos[node_ref][1]
        z_node = pos[node_ref][2]
        dist = np.sqrt((xv-x_node)**2+(yv-y_node)**2+(zv-z_node)**2) #this will be the distance matrix for node_ref
        distances3.append(np.concatenate((dist[0], dist[1], dist[2]), axis=None))
        li, lj, lk = np.where(dist < value3[pp])
        for i in range(len(li)):
             node_nr = int(li[i]* (yn * zn)  + lj[i] * zn + lk[i])
             if node_ref != node_nr:  # Avoid self-loops
                 edge_list.append([node_ref, node_nr])
             
    G3.add_edges_from(edge_list)
#%% circle bedingung
    c=0
    for i in range(len(pos)):
        bedingung3[i] = np.sqrt((pos[i][0]-max(x)/2)**2+(pos[i][1]-max(y)/2)**2+(pos[i][2]-max(z)/2)**2)
        if bedingung3[i] >= max(x)/2 + 0.01*max(x)/2:
            G3.remove_node(i)
        if bedingung3[i] < max(x)/2 + 0.01*max(x)/2:
            if bedingung3[i]>=(max(x)/2 + 0.01*max(x)/2)-value3[pp]:
                c = c+1
    #%% Analyse methoden
    number_of_nodes = G3.number_of_nodes()
    list_number_nodes_3d[pp] = number_of_nodes
    number_of_edges = G3.number_of_edges()
    list_number_edges_3d[pp] = number_of_edges
    nom_value3.append(value3[pp]/(länge/xn))
###degree unique
    degree_sequence3 = sorted((d for n, d in G3.degree()), reverse=True)
    list_degree_sequence3.append(degree_sequence3)
    max_degree3.append(max(degree_sequence3))
    min_degree3.append(min(degree_sequence3))
    av_degree3.append(sum(degree_sequence3)/number_of_nodes)
###clustering
    list_clustering_seq3.append(list(nx.clustering(G3).values()))
    av_clust3.append(nx.average_clustering(G3))
    
###path length
    av_path3d_list.append(nx.average_shortest_path_length(G3))#/total_nodes_3d[pp]
    shortest_path_lengths3 = dict(nx.all_pairs_shortest_path_length(G3))
    path_len_list.append(shortest_path_lengths3)
    diameter3[pp]     = max(nx.eccentricity(G3, sp=shortest_path_lengths3).values())
    path_lengths3 = np.zeros(int(diameter3[pp] + 1))
    for pls in shortest_path_lengths3.values():
        pl, cnts = np.unique(list(pls.values()), return_counts=True)
        path_lengths3[pl] += cnts
    freq_percent3[pp] = 100 * path_lengths3[1:] / path_lengths3[1:].sum()
### other
    density_3d.append(nx.density(G3))
    print(len(G3.nodes))
    skin_list3d.append(c)
#%% plot graph  with different layouts
# plt.clf()
# pos = nx.spring_layout(G3)
# # pos = nx.shell_layout(G3)
# # pos = nx.kamada_kawai_layout(G3)
# nx.draw_networkx(G3,pos=pos,with_labels=False,
#                   node_size=20,node_color='r',edge_color='grey')
# nx.draw_networkx_nodes(G3,pos=pos,node_size=20,node_color='r')

# plt.savefig('kawai 3D with edges.svg')
#%% plot by actual positions
# ### Extract node and edge positions from the layout
# node_xyz = np.array([pos[v] for v in sorted(G3)])
# edge_xyz = np.array([(pos[u], pos[v]) for u, v in G3.edges()])

# # Create the 3D figure
# fig = plt.figure(1)
# ax = fig.add_subplot(111, projection="3d")
# ax.view_init(elev=9, azim=-85, roll=-3)
# # Plot the nodes - alpha is scaled by "depth" automatically
# ax.scatter(*node_xyz.T, s=80, color='red')

# # Plot the edges
# for vizedge in edge_xyz:
#       ax.plot(*vizedge.T, color="tab:gray", linewidth=0.2)

# # ax.plot(*lon_3d.T, color='black')
# def _format_axes(ax):
#     """Visualization options for the 3D axes."""
#     # Turn gridlines off
#     ax.grid(True)
#     # Suppress tick labels
#     for dim in (ax.xaxis, ax.yaxis, ax.zaxis):
#         dim.set_ticks([])
#     # Set axes labels
#     ax.set_xlabel("x")
#     ax.set_ylabel("y")
#     ax.set_zlabel("z")
#     # ax.set_xlim(0,4)
# # Adding labels to nodes
# # for node, (x, y, z) in pos.items():
# #       ax.text(x, y, z, str(node), color="r", fontsize=12)


# _format_axes(ax)
# fig.tight_layout()
# plt.show()

# def update(i, fig, ax):
#     ax.view_init(elev=20., azim=i)
#     return fig, ax
# anim = FuncAnimation(fig, update, frames=np.arange(0, 360, 2), repeat=True, fargs=(fig, ax))
# anim.save('3d cube.gif', dpi=80, writer='Pillow', fps=24)
#safe fig
# plt.savefig('longest edge 3d value= 1.0 95°.svg')
#%%
# lon_3d= np.array([pos[0], pos[4]])
# lon_3d= np.array([pos[0], pos[21]])

# lon_3d= np.array([pos[0], pos[37]])
# lon_3d= np.array([pos[0], pos[27]])

#%%
# plt.clf()
# plt.figure(figsize=(6,5))
# plt.plot(value3,max_degree3,'r-.', label='max degree'  )
# plt.plot(value3,av_degree3 ,'r-' , label='average degree')
# plt.plot(value3,min_degree3,'r--', label='min degree'  )
# plt.xlabel('max. distance value', fontdict={"size": 12})
# plt.ylabel('degree', fontdict={"size": 12})
# plt.legend()

# # plt.savefig('3d degree vs value.svg')
#%%
# plt.clf()
# plt.figure(figsize=(6,5))
# plt.plot(list_number_nodes_3d,max_degree3,'r-.'  )
# plt.plot(list_number_nodes_3d,av_degree3 ,'r-' , label='3D')
# plt.plot(list_number_nodes_3d,min_degree3,'r--')
# plt.legend()
# # # plt.xlim(0,620)

# plt.xlabel('number of nodes', fontdict={"size": 12})
# plt.ylabel('degree', fontdict={"size": 12})
# plt.savefig('same dist all Ds degree vs N.svg')