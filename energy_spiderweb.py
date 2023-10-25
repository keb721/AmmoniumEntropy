################################################
#                                              #
#   CODE TO DETERMINE AVERAGE OXYGEN HASHES    #
#                                              #
################################################

# Import packages
# ---------------

import glob as gb
import numpy as np
import networkx as nx
import matplotlib.cm as cm
import scipy.constants as const
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll


# Parameters
# ----------

geometry = ['211'] #['112']
runs     = ['run1'] # Only deal with one of the repeats - can change if useful?
bases    = ['b3', 'b4', 'b5', 'b6', 'b7', ''] #['b2', 'b3', ''] #['b2', 'b3', 'b4', 'b5', 'b6', 'b7', '']
N_locs   = [[1], [2], [3], [12], [17], [20], [26], [31]]  # [[], [1], [2], [3], [20]]

conns_suffix  = "_connections.txt"        # Need to be careful with fenceposting
configs_pre   = "strucs"

plt.rcParams.update({"text.usetex": True})

# =============== #
# Define function #
# =============== #

def compare_hashes(conns_file, configs_file):
    ''' Function to compare the different abundance of hydrogen locations. '''

    # Get the connectivity matrix - produced by a FORTRAN code so need to deal with fenceposting
    
    conns = np.genfromtxt(conns_file)
    conns = conns[:, 2:]-1
    G     = nx.Graph()

    for i in range(len(conns)):
        # Create a graph of connections - because it's an undirected graph don't need to worry about duplicates
        for j in range(4):
            G.add_edge(i, conns[i, j])
        
    # First want to have a graph of the relative abundances of different hashes    
    # Divide a complete hash from data into information from each different oxygen
    
    delim = [5]+[4]*(len(conns)-1)
    data  = np.genfromtxt(configs_file, autostrip=True, dtype=str, comments="#",  delimiter=delim)
    
    # Create an aribitrary ordering of hashes for ease of analyis - count number of occurences
    # Can skip 0000 as no electron probability
    hashmap = ["0001", "0010", "0100", "1000", "1100", "1010", "1001", "0110", 
               "0101", "0011", "0111", "1011", "1101", "1110", "1111"]
    Ohash   = np.zeros((len(conns), 15))

    for a in range(len(conns)):
        for i, h in enumerate(hashmap):
            Ohash[a, i] = len(np.where(data[:, a]==h)[0])
            
    Ohash = Ohash/len(data)
    
    # Get probability of finding electron at each location
    total = np.zeros((len(conns), 4))
    for i in range(len(hashmap)):
        for j in range(4):
            if hashmap[i][j]=='1':
                total[:, j] += Ohash[:, i]
    
    return G, total, conns

        
def H_map(conns_file, configs_file, ions, basal=-1, lineres=100, show=False, save=True, pdf=True):
    ''' Plot a network showing the relative H density of a run. Ions is a dictionary of the locations
        of ions ONLY (i.e. cannot put in H2Os)'''

    G, total, conns = compare_hashes(conns_file, configs_file)
    # ======== #
    # Plotting #
    # ======== #
    
    fig, ax = plt.subplots(1, 1,figsize=(12.8, 9.6))
    plt.subplots_adjust(left=0.01, right=0.975, bottom=0.025, top=0.975, wspace=0.1, hspace=0.1)
    
    # Get positions of the graph - might need to try out a few positions
        
    pos = nx.spectral_layout(G)
    
    Oxygens = list(np.linspace(0, len(conns)-1, len(conns), dtype=int))

    F_locs   = []
    OH_locs  = [] 
    H3O_locs = []
    N_locs   = []

    for key in ions.keys():
        loc = key - 1
        Oxygens.remove(loc)
        if ions[key] == 'F':
            F_locs.append(loc)
        elif ions[key] == 'OH':
            OH_locs.append(loc)
        elif ions[key] == 'H3O':
            H3O_locs.append(loc)
        elif ions[key] == 'N':
            N_locs.append(loc)
        else:
            print('ERROR: Incorrectly defined ion. Note: do not include H2O')
            exit

    F_locs   = np.asarray(F_locs)
    OH_locs  = np.asarray(OH_locs) 
    H3O_locs = np.asarray(H3O_locs)
    N_locs   = np.asarray(N_locs)

    nx.draw_networkx_nodes(G, pos, nodelist=Oxygens, node_color='tab:orange', label="O")
    nx.draw_networkx_nodes(G, pos, nodelist=list(np.linspace(len(conns), np.max(conns), int(np.max(conns)-len(conns)+1), dtype=int)), node_color='tab:orange', alpha=0.5, label="ghost")
    nx.draw_networkx_nodes(G, pos, nodelist=F_locs, label="F", node_color='gold')
    nx.draw_networkx_nodes(G, pos, nodelist=OH_locs, label="OH", node_color='m')
    nx.draw_networkx_nodes(G, pos, nodelist=H3O_locs, label="H3O", node_color='darkcyan')
    nx.draw_networkx_nodes(G, pos, nodelist=N_locs, label="N", node_color='g')
    nx.draw_networkx_nodes(G, pos, nodelist=[basal-1], label="base", node_color='tab:grey')

    ax  = plt.gca()

    # Now need to draw connections
    cmap = plt.get_cmap('seismic')
    for i in range(len(conns)):
        for j in range(4):
            x = np.linspace(pos[i][0], pos[int(conns[i, j])][0], lineres)
            y = np.linspace(pos[i][1], pos[int(conns[i, j])][1], lineres)
            z = np.linspace(total[i, j], 1-total[i, j], lineres)
            
            tot = np.array([x, y]).T.reshape(-1, 1, 2)
            tot = np.concatenate([tot[:-1], tot[1:]], axis=1) 
            lc  = mcoll.LineCollection(tot, array=z, cmap=cmap, norm=plt.Normalize(0.0, 1.0), linewidth=3)
            ax.add_collection(lc)



    sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0.0, 1.0))
    sm.set_array(np.linspace(0, 1, lineres))
    cb = plt.colorbar(sm, label=r"$P_H$")
    cb.ax.tick_params(labelsize=20)
    cb.set_label(r"$P_H$", size=25)
    #plt.title(configs_file.split("strucs_")[1]+ " without O only")
    plt.legend(fontsize=20, loc=(0.8, 0.55))

    filename = configs_file.split('.')[0]
    plt.title(filename.split("strucs_")[1])
    if save:
        if pdf:
            plt.savefig("bond"+filename.split("strucs")[1]+".pdf")
        else:
            plt.savefig("bond"+filename.split("strucs")[1]+".png")
    if show:
        plt.show()
    if save and not show:
        plt.close(fig)
    return None



def comp_H_map(conns_file, configs_file, no_N_file, ions, basal=-1, lineres=100, show=False, save=True, pdf=True):
    ''' Plot a network showing the relative H density of a run compared to a run with no ions. 
        Ions is a dictionary of the locations of ions ONLY (i.e. cannot put in H2Os)'''

    G, total, conns = compare_hashes(conns_file, configs_file)
    _, no_N_tot, _  = compare_hashes(conns_file, no_N_file)

    # ======== #
    # Plotting #
    # ======== #
    
    fig, ax = plt.subplots(1, 1,figsize=(12.8, 9.6))
    plt.subplots_adjust(left=0.01, right=0.975, bottom=0.025, top=0.975, wspace=0.1, hspace=0.1)

    # Get positions of the graph - might need to try out a few positions
    
    if (conns_file.split('_')[0]=='112'):
#        pos = nx.spring_layout(G) # I think this is slightly worse than kamada_kawai but not certain
        pos = nx.kamada_kawai_layout(G)  # Best of only 2 ok graphs for 112 - spectral gives overlaps so bad
    else:
        pos = nx.spectral_layout(G) # Gives a nice graph for 211 - need to experiment for others
    
    Oxygens = list(np.linspace(0, len(conns)-1, len(conns), dtype=int))

    F_locs   = []
    OH_locs  = [] 
    H3O_locs = []
    N_locs   = []

    for key in ions.keys():
        loc = key - 1
        Oxygens.remove(loc)
        if ions[key] == 'F':
            F_locs.append(loc)
        elif ions[key] == 'OH':
            OH_locs.append(loc)
        elif ions[key] == 'H3O':
            H3O_locs.append(loc)
        elif ions[key] == 'N':
            N_locs.append(loc)
        else:
            print('ERROR: Incorrectly defined ion. Note: do not include H2O')
            exit

    F_locs   = np.asarray(F_locs)
    OH_locs  = np.asarray(OH_locs) 
    H3O_locs = np.asarray(H3O_locs)
    N_locs   = np.asarray(N_locs)

    nx.draw_networkx_nodes(G, pos, nodelist=Oxygens, node_color='tab:orange', label="O")
    nx.draw_networkx_nodes(G, pos, nodelist=list(np.linspace(len(conns), np.max(conns), int(np.max(conns)-len(conns)+1), dtype=int)), node_color='tab:orange', alpha=0.5, label="ghost")
    nx.draw_networkx_nodes(G, pos, nodelist=F_locs, label="F", node_color='gold')
    nx.draw_networkx_nodes(G, pos, nodelist=OH_locs, label="OH", node_color='m')
    nx.draw_networkx_nodes(G, pos, nodelist=H3O_locs, label="H3O", node_color='darkcyan')
    nx.draw_networkx_nodes(G, pos, nodelist=N_locs, label="N", node_color='g')
    nx.draw_networkx_nodes(G, pos, nodelist=[basal-1], label="base", node_color='tab:grey')

    ax  = plt.gca()

    diff_tot = total - no_N_tot

    # Now need to draw connections
    cmap = plt.get_cmap('seismic')
    for i in range(len(conns)):
        for j in range(4):
            x = np.linspace(pos[i][0], pos[int(conns[i, j])][0], lineres)
            y = np.linspace(pos[i][1], pos[int(conns[i, j])][1], lineres)
            z = np.linspace(diff_tot[i, j], 0-diff_tot[i, j], lineres)
            
            tot = np.array([x, y]).T.reshape(-1, 1, 2)
            tot = np.concatenate([tot[:-1], tot[1:]], axis=1) 
            lc  = mcoll.LineCollection(tot, array=z, cmap=cmap, norm=plt.Normalize(-1.0, 1.0), linewidth=4)
            ax.add_collection(lc)


    sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-1.0, 1.0))
    sm.set_array(np.linspace(0, 1, lineres))    
    cb = plt.colorbar(sm, label=r"$\Delta P_H$ (cf. pure ice)")
    cb.ax.tick_params(labelsize=20)
    cb.set_label(r"$\Delta P_H$ (cf. pure ice)", size=25)
    plt.legend(fontsize=20, loc=(0.8, 0.55))

    plt.savefig("test_bonds.png", dpi=300)
    if save:
        filename = configs_file.split('.')[0]
        if pdf:
            plt.savefig("comp_bond"+filename.split("strucs")[1]+".pdf")
        else:
            plt.savefig("comp_bond"+filename.split("strucs")[1]+".png")
    if show:
        plt.show()
    if save and not show:
        plt.close(fig)
    return None
