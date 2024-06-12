#%%mini program to generate a 2D grid of nodes
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

# plt.clf()
# N3d                   = [216]#[64]#,64,64]#
# total_nodes_2d        = [100]#,100,100]
# value2                = [1.2,2.8,0.8]
# vergleich gleiche value
# N3d                   = [1331,1000,512]#
# total_nodes_2d        = [676,484,256]
# value2                = [2.8,2.8,2.8,2.8,2.8]
value2                = [1.8,1.8,1.8,1.8,1.8,1.8]

# vergleich gleiche nodes unterschiedliche values
# N3d                   = [729,729,729,729]#
# total_nodes_2d        = [400,400,400,400]#
# value2                = [0.8,1.1,1.8,2.8]#
#v vergleich gleiche degree range
N3d                   = [64,216,512,729,1331,2744]#64,216,
total_nodes_2d        = [25,100,256,361,729,1369]#25,100,
value2                = [1.5,1.4,1.2,1.16,1.01, 1]
av_clust2             = []
density_2d            = []
abstände2            = [] 
skin_list2d=[]
av_path2d_list        = []
max_degree2           = []
min_degree2           = []
av_degree2            = []
list_degree_sequence2 = []
list_clustering_seq2  = []
list_clust_bins2d     = []
path_len2_list        = []

diameter2     = {}
freq_percent2 = {}
node_list_corner=[]
list_number_edges_2d  = np.zeros(len(total_nodes_2d))
list_number_nodes_2d  = np.zeros(len(total_nodes_2d))
#%% loop over total node list
for pp in range(len(total_nodes_2d)):
    distances2=[]
    c=0
    bedingung2 = np.zeros(total_nodes_2d[pp])
#%% mesh grid und pos ini
    scale_x, scale_y = 1.0, 1.0
    xn, yn = (round(np.sqrt(total_nodes_2d[pp]),0), round(np.sqrt(total_nodes_2d[pp]),0)) ##number of nodes in x,y directions
    länge     =   round(np.cbrt(N3d[pp]),0)-1  # xn-1#                       # last coordinate of graph 
    x = scale_x*np.linspace(0,länge,int(xn), endpoint=True)                    ##x-coordinates
    y = scale_y*np.linspace(0,länge,int(xn), endpoint=True)                    ##y-coordinates 
    xv, yv = np.meshgrid(x, y)                                                  ## meshgrid der node positions
    xvr = np.reshape(xv,xv.size)                                                ## liste der x coordinate des meshgrid
    yvr = np.reshape(yv,yv.size)                                                ## liste der y coordinate des meshgrid 
    random1=np.random.normal(loc=0,scale=(länge/int(xn))/2,size=total_nodes_2d[pp])   #(länge/int(xn))/2       ## random variations of distance
    random2=np.random.normal(loc=0,scale=(länge/int(xn))/2,size=total_nodes_2d[pp])          ## random variations of distance
    pos2={}                                                                     ##built position dictionary from x,y coordinate
    for i in range(len(xvr)):                                                   ##loop over length of mesh coord list
        pos2[i] = [xvr[i],yvr[i]]                                               ## set node position dict not random
## wenn random position        
        # pos2[i] = [xvr[i]+random1[i],yvr[i]+random2[i]]
        # xvr[i] += random1[i]
        # yvr[i] += random2[i]
    xv = np.reshape(xvr,xv.shape)
    yv = np.reshape(yvr,yv.shape)
#%% nodes und edges erzeugen
    G2 = nx.path_graph(xvr.size)
    G2.remove_edges_from(list(G2.edges))
#%% circle bedingung
    for i in range(len(pos2)):
        bedingung2[i] = np.sqrt((pos2[i][0]-max(x)/2)**2+(pos2[i][1]-max(y)/2)**2)
        if bedingung2[i] >= max(x)/2 + 0.01*max(x)/2:
            node_list_corner.append(i)
            G2.remove_node(i)
        if bedingung2[i] < max(x)/2 + 0.01*max(x)/2:
            if bedingung2[i]>=(max(x)/2 + 0.01*max(x)/2)-value2[pp]:
                c = c+1
#%%                
    node_list=[]
    for node_ref in G2.nodes(): #loop over each node
        x_node = pos2[node_ref][0] # x pos of referenz node
        y_node = pos2[node_ref][1] # y pos of referenz node
        dist = np.sqrt((xv-x_node)**2+(yv-y_node)**2) #this will be the distance matrix for node_ref
        distances2.append(np.concatenate((dist[0], dist[1]), axis=None))
        li, lj = np.where(dist < value2[pp])
        for i in range(len(li)):
            node_nr = int(li[i]*xn+lj[i])
            if node_ref != node_nr:
                node_list.append([node_ref, node_nr])
    G2.add_edges_from(node_list)

    #%% degree analysis
    number_of_nodes = G2.number_of_nodes()
    list_number_nodes_2d[pp] = number_of_nodes
    number_of_edges = G2.number_of_edges()
    list_number_edges_2d[pp] = number_of_edges
    
    degree_sequence2 = sorted((d for n, d in G2.degree()), reverse=True)
    list_degree_sequence2.append(degree_sequence2)
    max_degree2.append(max(degree_sequence2))
    min_degree2.append(min(degree_sequence2))
    av_degree2.append(sum(degree_sequence2)/number_of_nodes)
    #clustering_sequence2, clust_bins2d, patch = plt.hist(nx.clustering(G2).values())
    list_clustering_seq2.append(list(nx.clustering(G2).values()))
    av_clust2.append(nx.average_clustering(G2))
    # path len
    av_path2d_list.append(nx.average_shortest_path_length(G2))#/(total_nodes_2d[pp])
    shortest_path_lengths2 = dict(nx.all_pairs_shortest_path_length(G2))
    path_len2_list.append(shortest_path_lengths2)
    diameter2[pp]     = max(nx.eccentricity(G2, sp=shortest_path_lengths2).values())
    path_lengths2 = np.zeros(diameter2[pp] + 1, dtype=int)
    for pls in shortest_path_lengths2.values():
        pl, cnts = np.unique(list(pls.values()), return_counts=True)
        path_lengths2[pl] += cnts
    freq_percent2[pp] = 100 * path_lengths2[1:] / path_lengths2[1:].sum()
### other
    density_2d.append(nx.density(G2))
    print(len(G2.nodes))
    abstände2.append(dist[0][0]-dist[1][0])
    skin_list2d.append(c)
#%% plot graph
# plt.clf()

# pos2 = nx.spring_layout(G2)
# pos2 = nx.shell_layout(G2)
# pos2 = nx.kamada_kawai_layout(G2)
# nx.draw_networkx(G2,pos=pos2,with_labels=False,
#                   node_size=20,node_color='b',edge_color='grey')
# nx.draw_networkx_nodes(G2,pos=pos2,node_size=20,node_color='b')

# plt.savefig('kawai 2D .svg')
#%% plot by actual positions
# plt.clf()
# nx.draw_networkx(G2,pos=pos2, with_labels=True,node_size=10,node_color='blue',edge_color='grey',label='circular 2D Graph') ##here, labels correspond to node_id
# nx.draw_networkx_nodes(G3,pos=pos3,node_size=16,node_color='red',label='spherical 3D Graph') ##here, labels correspond to node_id

# nx.draw_networkx_nodes(G2,pos=pos2,nodelist=(list(G2.nodes)),node_size=16,node_color='blue',label='circular 2D Graph') ##here, labels correspond to node_id
# nx.draw_networkx_nodes(G2,pos=pos2,nodelist=node_list_corner,node_size=(3.5),node_color='k', label='removed nodes')
# plt.gca().set_aspect('equal', adjustable='box')

# # nx.draw_networkx_edges(G2,pos=pos2,edgelist=[(0,10)],edge_color='k')
# plt.plot(max(x)/2,max(y)/2,'xr',label='center point')

# plt.tight_layout()
# plt.legend(loc='upper right')
# nx.draw_networkx_edges(G2,pos=pos2,edgelist=(list(G2.edges.data())),arrowsize=3,edge_color='grey')

# plt.savefig('circle perturbated .svg')

#%%
# plt.clf()
# plt.figure(figsize=(6,5))
# plt.plot(value2,max_degree2,'b-.', label='max degree'  )
# plt.plot(value2,av_degree2 ,'b-' , label='average degree')
# plt.plot(value2,min_degree2,'b--', label='min degree'  )
# plt.xlabel('max. distance value', fontdict={"size": 12})
# plt.ylabel('degree', fontdict={"size": 12})
# plt.legend()
# plt.savefig('degree vs value.svg')

#%%
# plt.clf()
# plt.figure(figsize=(6,5))
# plt.plot(list_number_nodes_2d,max_degree2,'b-.' )
# plt.plot(list_number_nodes_2d,av_degree2 ,'b-' , label='2D')
# plt.plot(list_number_nodes_2d,min_degree2,'b--' )
# plt.legend()
# plt.savefig('degree vs value.svg')