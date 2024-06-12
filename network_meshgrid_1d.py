#%% initialize
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

#%% listen ini
# plt.clf()

# total_nodes_1d        = [60]#,66,66]#
# value1                = [2.8,2.8,3.8]
# N3d                   = [216]#,64,64]#
# # gleiche value
# total_nodes_1d        = [510,350,170]#
value1                = [1.8,1.8,1.8,1.8,1.8,1.8]
# value1                = [2.8,2.8,2.8,2.8,2.8,2.8]
# N3d                   = [1331,1000,512]#
#vergleich gleiche nodes unterschiedliche values
# N3d                   = [729,729,729,729]#
# total_nodes_1d        = [290,290,290,290]
# value1                = [1,0.6, 0.39,0.2]
# # vergleich gleiche degree ranges
N3d                     = [27,216,512,729,1331,2744]
total_nodes_1d          = [19,64,200,290,540,1050]
value1                = [0.32,0.29,0.25,0.24,0.19,0.13]
abstände1=[]
av_clust1             = []
skin_list1d           = []
density_1d            = []
av_path1d_list        = []
list_degree_sequence1 = []
max_degree1           = []
min_degree1           = []
av_degree1            = []
list_clustering_seq1  = []
list_clust_bins1d     = []
path_len1_list        = []
diameter1     = {}
freq_percent1 = {}
list_number_edges_1d  = np.zeros(len(total_nodes_1d))
list_number_nodes_1d  = np.zeros(len(total_nodes_1d))
#%% loop over node number list
for pp in range(len(total_nodes_1d)):
    bedingung1   = np.zeros((total_nodes_1d[pp]))
    #%% meshgrid setup position
    scale_x     = 1.0
    xn          = (total_nodes_1d[pp])                                                   ##number of nodes in x directions
    länge       =   round(np.cbrt(N3d[pp]),0)-1 #xn-1# 
    x           = scale_x*np.linspace(0,länge, int(xn), endpoint=True)                   ##x-coordinates 
    xv          = np.meshgrid(x)                                                         ## meshgrid of node position
    xvr         = list(xv[0][:])                                                          ## mesh als liste umgewandelt
    pos1={}                                                                     ## pos dict ini
    random= np.random.normal(loc=0,scale= 0,size=total_nodes_1d[pp]) # (länge/int(xn))/2 random variations of distance
    for i in range(len(xvr)):                                                   ## loop over pos
        pos1[i]= [0,xvr[i]+random[i]]                                         ## set each position for the graph with
        xv[0][i] += random[i]                                                   ## important to pertubate the meshgrid of position for the dist calculation
    #%% graph add nodes and egdes
    G = nx.path_graph(len(xvr))                                                 ## create graph length mesh list
    G.remove_edges_from(list(G.edges))                                          ## remove all edges    
    node_list=[]                                                                ## list of nodes ini for add.edges.from
    for node_ref in G.nodes():                                                  ## loop over each node
        x_node = pos1[node_ref][1]                                              ## get x coord to calculate dist
        dist = (abs(xv-x_node))                                                 ## all distances to node_ref
        li,lj = np.where(dist < value1[pp])                                     ## find elements smaler than value
        for i in range(len(li)):                                                ## loop over indices li 
            node_nr = int(li[i]*xn+lj[i])                                       ## 
            if node_ref != node_nr:
                node_list.append([node_ref, node_nr])
    
    G.add_edges_from(node_list)
    c=0
#%% circle bedingung
    for i in range(len(pos1)):
        bedingung1[i] = np.sqrt((pos1[i][1]-max(x)/2)**2)
        if bedingung1[i] >= max(x)/2 +0.01*max(x)/2:
            G.remove_node(i)
        if bedingung1[i] < max(x)/2 + 0.01*max(x)/2:
            if bedingung1[i]>=(max(x)/2 + 0.01*max(x)/2)-value1[pp]:
                c = c+1
    #%% degree analyse
    number_of_nodes = G.number_of_nodes()
    list_number_nodes_1d[pp] = number_of_nodes
    number_of_edges = G.number_of_edges()
    list_number_edges_1d[pp] = number_of_edges

#degree unique
    
    degree_sequence1 = sorted((d for n, d in G.degree()), reverse=True)
    #dmax2 = max(degree_sequence1)
    list_degree_sequence1.append(degree_sequence1)
    max_degree1.append(max(degree_sequence1))
    min_degree1.append(min(degree_sequence1))
    av_degree1.append(sum(degree_sequence1)/number_of_nodes)
    #clustering_sequence1, clust_bins1d, patch = plt.hist(nx.clustering(G).values())
    list_clustering_seq1.append(list(nx.clustering(G).values()))
    #list_clust_bins1d.append(clust_bins1d)
    av_clust1.append(nx.average_clustering(G))
    shortest_path_lengths1 = dict(nx.all_pairs_shortest_path_length(G))
    av_path1d_list.append(nx.average_shortest_path_length(G))#/(total_nodes_2d[pp])
    path_len1_list.append(shortest_path_lengths1)
    diameter1[pp]     = max(nx.eccentricity(G, sp=shortest_path_lengths1).values())
    path_lengths1 = np.zeros(diameter1[pp] + 1, dtype=int)
    for pls in shortest_path_lengths1.values():
        pl, cnts = np.unique(list(pls.values()), return_counts=True)
        path_lengths1[pl] += cnts
    freq_percent1[pp] = 100 * path_lengths1[1:] / path_lengths1[1:].sum()
### other
    density_1d.append(nx.density(G))
    print(len(G.nodes))
    abstände1.append(dist[0][0]-dist[0][1])
    skin_list1d.append(c)

#%% plot graph
# plt.clf()

# pos1 = nx.spring_layout(G)
# pos1 = nx.shell_layout(G)
# pos1 = nx.kamada_kawai_layout(G)
# nx.draw_networkx(G,pos=pos1,with_labels=False,
#                   node_size=20,node_color='y',edge_color='grey')
# nx.draw_networkx_nodes(G,pos=pos1,node_size=20,node_color='y')

# plt.savefig('kawai 123D .svg')
#%% plot by actual positions
# plt.clf()
# nx.draw_networkx(G ,pos=pos1,node_size=10,node_color='red',label='linear 1D Graph',with_labels=False) ##here, labels correspond to node_id!
# nx.draw_networkx_nodes(G ,pos=pos1,node_size=10,node_color='orange',label='linear 1D Graph') ##here, labels correspond to node_id!
# # plt.grid(visible=True)
# # plt.legend()
# # # plt.xlim(-1,1)
# # plt.ylim(0,3)
# nx.draw_networkx_edges(G ,pos=pos1,edgelist=[(0,60),(0,39),(0,21)],edge_color='k',arrows=True, connectionstyle='arc3,rad=0.7') ##here, labels correspond to node_id!
# # plt.show()
# plt.savefig('1D alongest beide .svg')

#%%
# plt.clf()
# plt.figure(figsize=(6,5))
# plt.plot(0,0,color='k',linestyle='-.',label='max. degree')
# plt.plot(0,0,color='k',linestyle='-',label='average degree')
# plt.plot(0,0,color='k',linestyle='--',label='min. degree')
# plt.plot(list_number_nodes_1d,max_degree1,color='orange',linestyle='-.'  )
# plt.plot(list_number_nodes_1d,av_degree1 ,color='orange',linestyle='-', label='1D')
# plt.plot(list_number_nodes_1d,min_degree1,color='orange',linestyle='--' )

# plt.legend()
# plt.savefig('1D degree vs value.svg')