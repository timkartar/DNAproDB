import networkx as nx
from matplotlib import pyplot as plt
import numpy as np
import sre_yield
from copy import deepcopy
from math import cos, sin
from get_helix_coords import process_resid
import time
import random
import os, sys

FIG_PATH = ''
MEDIA_PATH = ''

plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
style_dict = {}
arrow_dict = {}

log = []

chem_components = dict(np.load(os.path.dirname(os.path.abspath(__file__)) + "/modified_parents.npz",allow_pickle=True))
#python /home/aricohen/Desktop/rnascape/run.py uploads/7vnv-assembly1.cif 7vnv 1 rnascape uploads/7vnv-assembly1.cif.out

"""
Get index of a nucleotide from chain id and resid
"""
def getIndex(target_chid, target_resid, ids, chids):
    matching_ids_with_indices = [(index, sublist) for index, sublist in enumerate(ids) if str(sublist[1]) == str(target_resid)]
    for tuple in matching_ids_with_indices:
        cur_index = tuple[0]
        if chids[cur_index] == target_chid:
            return cur_index

"""
Use DSSR to get LW annotations rather than rnascape file upload
"""
def getCustomMarker(pos, item):
    marker = item[1][pos]
    m = item[0][0]
    
    if pos == 0 :
        t = item[0][1]
    else:
        t = item[0][2]
    
    d = t - m
    d = d/np.linalg.norm(d)
    v = [0,1]
    angle = np.rad2deg(np.arccos(np.dot(d,v)))
    tanangle = np.rad2deg(np.arctan(d[1]/d[0]))

    if marker == ">":
        #print(angle, tanangle, item[1])

        if tanangle < 0 and angle <= 90:
            return (3,0, angle)
        if tanangle < 0 and angle > 90:
            return (3,0, -1*angle)
        
        if tanangle >= 0 and angle > 90:
            return (3,0, angle)
        if tanangle >= 0 and angle <= 90:
            return (3,0, -1*angle)
    
    if marker == "s":
        if tanangle < 0 and angle <= 90:
            return (4,0, angle + 45)
        if tanangle < 0 and angle > 90:
            return (4,0, -1*angle + 45)
        
        if tanangle >= 0 and angle > 90:
            return (4,0, angle + 45)
        if tanangle >= 0 and angle <= 90:
            return (4,0, -1*angle + 45)
    else:
        return marker

def getBasePairingEdges(dssrout, dssrids, points):
    
    #bp types: DSSR [ct][MWm][+-][MWm]
    dssr_bp_types = list(sre_yield.AllStrings('[WMm][WMm]'))
    #['MM', 'WM', 'mM', 'MW', 'WW', 'mW', 'Mm', 'Wm', 'mm']


    bp_marker_types = ['oo','so', '>o', 'os', 'ss', '>s', 'o>', 's>', '>>']
    #marker_bp_types = list(sre_yield.AllStrings(bp_marker_types))
    
    bp_map = {}
    for item in dssr_bp_types:
        bp_map[item] = bp_marker_types[dssr_bp_types.index(item)]
    
    # print(bp_map)
    
    magnification = max(1, min(len(dssrids)/40, 10))
    edges = []
    bp_markers = []

    for item in dssrout['pairs']:
        i1 = dssrids.index(item['nt1'])
        i2 = dssrids.index(item['nt2'])
        edges.append((i1, i2))
        style_dict[(i1, i2)] = ':'#'dashed'
        arrow_dict[(i1, i2)] = 0.001*magnification

        v = points[i1] - points[i2]
        d = np.linalg.norm(v)

        if d < 2: ## do not show bp type for too small edges
            continue
        p = points[i2]+v/2 
        typ = item['DSSR'][1]+ item['DSSR'][3]
        
        direc = v / d
        SCALAR = 0.8
        p2 = p + (-1 * direc * SCALAR)
        p1 = p + (direc * SCALAR)
        p=[p,p1,p2] # use p if hoog/hoog or sugar/sugar, otherwise p1 and p2


        if "." in typ: # DSSR couldn't determine properly
            continue
        if typ == "WW" and item['LW'][0] == 'c': #do not show cis watson crick pairs
            continue
        if item['DSSR'][0] == 'c':
            orient = 'k'
        else:
            orient = 'w'
        bp_markers.append([p, bp_map[typ], orient, item['DSSR'][0]+typ])

        
    return edges, bp_markers, bp_map


def getBasePairingEdgesDssrLw(dssrout, dssrids, points):
    dssr_lw_bp_types = list(sre_yield.AllStrings('[WHS][WHS]'))
    # dssr_lw_bp_types.remove("WW")
    # ['HW', 'SW', 'WH', 'HH', 'SH', 'WS', 'HS', 'SS']
    #dssr_lw_markers = ['so', '<o', 'os', 'ss', '<s', 'o>', 's>', '<>']
    dssr_lw_markers = ['oo','so', '>o', 'os', 'ss', '>s', 'o>', 's>', '>>']
    bp_map = {}
    for item in dssr_lw_bp_types:
        bp_map[item] = dssr_lw_markers[dssr_lw_bp_types.index(item)]
    magnification = max(1, min(len(dssrids)/40, 10))
    edges = []
    bp_markers = []

    for item in dssrout['pairs']:
        i1 = dssrids.index(item['nt1'])
        i2 = dssrids.index(item['nt2'])
        edges.append((i1, i2))
        style_dict[(i1, i2)] = ':'#'dashed'
        arrow_dict[(i1, i2)] = 0.001*magnification

        v = points[i1] - points[i2]
        d = np.linalg.norm(v)

        if d < 2: ## do not show bp type for too small edges
            continue
        p = points[i2]+v/2 
        typ = item['LW'][1]+ item['LW'][2]
        

        # Compute points for each shape based on directional vector
        direc = v / d
        SCALAR = 0.8
        p2 = p + (-1 * direc * SCALAR)
        p1 = p + (direc * SCALAR)
        p=[p,p1,p2] # use p if hoog/hoog or sugar/sugar, otherwise p1 and p2


        if "." in typ: # DSSR couldn't determine properly
            continue
        if typ == "WW" and item['LW'][0] == 'c': #do not show cis watson crick pairs
            continue
        if item['LW'][0] == 'c':
            orient = 'k'
        else:
            orient = 'w'
        bp_markers.append([p, bp_map[typ], orient, item['LW'][0]+typ])
    return edges, bp_markers, bp_map

    # Get index of second nucleotide

    #     i1 = dssrids.index(item['nt1'])
    #     i2 = dssrids.index(item['nt2'])
    #     edges.append((i1, i2))
    #     style_dict[(i1, i2)] = ':'#'dashed'
    #     arrow_dict[(i1, i2)] = 0.001*magnification

    #     v = points[i1] - points[i2]
    #     d = np.linalg.norm(v)

    #     if d < 2: ## do not show bp type for too small edges
    #         continue
    #     p = points[i2]+v/2 
    #     typ = item['DSSR'][1]+ item['DSSR'][3]
        
    #     if "." in typ: # DSSR couldn't determine properly
    #         continue
    #     if typ == "WW": #do not show watson crick pairs
    #         continue
    #     if item['DSSR'][0] == 'c':
    #         orient = 'k'
    #     else:
    #         orient = 'w'
    #     bp_markers.append([p, bp_map[typ], orient, item['DSSR'][0]+typ])


def getBasePairingEdgesSaenger(dssrout, dssrids, points):
    #bp types: DSSR [ct][MWm][+-][MWm]
    dssr_bp_types = list(sre_yield.AllStrings('[012][0-9]'))
    l = [int(i) for i in dssr_bp_types]
    idx = np.argsort(l)
    dssr_bp_types = np.array(dssr_bp_types)[idx].tolist()
    # print(dssr_bp_types)
    
    #dssr_bp_types.remove("WW")
    
    #bp_marker_types = '[o^pdshP*]'
    #marker_bp_types = list(sre_yield.AllStrings(bp_marker_types))
    
    bp_map = {}
    for item in dssr_bp_types:
        bp_map[item] = item #'${}$'.format(item)
    
    
    magnification = max(1, min(len(dssrids)/40, 10))
    edges = []
    bp_markers = []

    for item in dssrout['pairs']:
        i1 = dssrids.index(item['nt1'])
        i2 = dssrids.index(item['nt2'])
        edges.append((i1, i2))
        style_dict[(i1, i2)] = ':'#'dashed'
        arrow_dict[(i1, i2)] = 0.001*magnification

        v = points[i1] - points[i2]
        d = np.linalg.norm(v)

        if d < 2: ## do not show bp type for too small edges
            continue
        p = points[i2]+v/2 
        typ = item['Saenger'].split("-")[0]
        
        if typ not in dssr_bp_types: # DSSR couldn't determine properly
            continue
        orient = 'w'
        bp_markers.append([p, bp_map[typ], orient, typ])

        
    return edges, bp_markers, bp_map


def getBackBoneEdges(ids, chids, dssrids, dssrout):

    magnification = max(1, min(len(ids)/40, 10))
    edges = []
    for i in range(len(dssrids) -1):
        #_,_,_, chid1 = process_resid(dssrout['nts'][i]["nt_id"])
        #_,_,_, chid_next = process_resid(dssrout['nts'][i+1]["nt_id"])
        chid1 = dssrout['nts'][i]["nt_id"].split(".")[2]#process_resid(dssrout['nts'][i]["nt_id"], model)
        chid_next = dssrout['nts'][i+1]["nt_id"].split(".")[2]#process_resid(dssrout['nts'][i+1]["nt_id"], model)
        
        try:
            if chid1 == chid_next:
                a = dssrids.index(dssrout['nts'][i]["nt_id"])
                b = dssrids.index(dssrout['nts'][i+1]["nt_id"])
                edges.append((a,b))
                style_dict[(a,b)] = 'solid'
                arrow_dict[(a,b)] = 10*np.sqrt(magnification)
        except:
            pass
    
    return edges

def getResNumPoints(points, ids, G, k=10, separation=1):
    from sklearn.neighbors import KDTree
    tree = KDTree(points)
    SCALAR = 2*separation
    rett = []
    labels = []
    for i in range(0,len(points), k):
        p = points[i]
        l = str(ids[i][1]) + ids[i][2].replace(" ","")
        out = tree.query_radius([p], r=10)[0][1:].tolist()
        if len(out) > 20:
            continue ## too congested
        out += [item[1] for item in G.edges(i)]
        neighbors = (np.array(points)[out])
        c = neighbors.mean(axis=0)
        try:
            ret = p + SCALAR*(p-c)/np.linalg.norm(p - c)
        except:
            continue
        rett.append(ret)
        labels.append(l)
    return rett, labels

def Plot(points, markers, ids, chids, dssrids, dssrout, prefix="", rotation=False, bp_type='DSSR',
        out_path=None, time_string="ac1", extra={'arrowsize':1, 'circlesize':1,
            'circle_labelsize':1, 'cols':['#FF9896', '#AEC7E8', '#90CC84', '#DBDB8D', '#FFFFFF'],
            'showNumberLabels': True, 'numberSeparation': 1, 'numberSize': 1, 'markerSize': 1
            }, mFIG_PATH='', mMEDIA_PATH=''):
    '''rotation is False if no rotation is wished, otherwise, one
    can a pass a value in radian e.g. np.pi , np.pi/2, np.pi/3 etc. '''
    FIG_PATH = mFIG_PATH
    MEDIA_PATH = mMEDIA_PATH
    
    BP_EDGE_COLOR = 'royalblue' ## base-pairing edge color
    BB_EDGE_COLOR = 'k' ## back bone edge color
    BP_MARKER_COLOR = 'royalblue' ## base-pairing annotation marker
    
    edge_color_dict = {}

    marker_size = float(extra['markerSize']) #default 1
    dssrids = list(dssrids) # for npz
    rotation_string = "" # used to append to file name
    if not rotation:
        pass
    else:
        rotation_string = str(rotation)
        rotation = np.radians(float(rotation))
        centroid = np.mean(points, axis=0)
        V = points - centroid
        theta = -1*rotation
        rot = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
        V_ = np.dot(rot, V.T).T
        points = centroid + V_

    magnification = max(1, min(len(points)/40,10))
    G = nx.DiGraph()
    cold = {'A': extra['cols'][0],
    'C': extra['cols'][1],    'G': extra['cols'][2],
    'U': extra['cols'][3],    'DA': extra['cols'][0],
    'DC': extra['cols'][1],    'DG': extra['cols'][2],
    'DT': extra['cols'][3],    }

    colors = []
    labels = {}
    for i, (point, marker) in enumerate(zip(points, markers)):
        G.add_node(i, pos=point)
        marker = marker.replace("$","")
        
        if marker in cold.keys():
            labels[i] = marker[-1]
        else:
            #print(chem_components['0MC'])
            if marker in chem_components.keys():
                parent = chem_components[marker].tolist()
                labels[i] = parent[-1].lower()
            else:
                labels[i] = 'X'
            log.append("Non standard residue {} ({}, chain {}) assigned label {}.".format(marker, ids[i], chids[i], labels[i]))
        try:
            colors.append(cold[marker])
        except:
            colors.append(extra['cols'][4])
    

    fig, ax = plt.subplots(1, figsize=(8*magnification, 6*magnification))
    
    
    #### draw edges #######
    # if no pairs do not get base pairing edges:
    pairings = None
    bp_markers = None
    bp_map = None
    
    edgecolors = []
    if 'pairs' in dssrout.keys():
        if bp_type == "dssr":
            pairings, bp_markers, bp_map = getBasePairingEdges(dssrout, dssrids, points)
        elif bp_type == "saenger":
            pairings, bp_markers, bp_map = getBasePairingEdgesSaenger(dssrout, dssrids, points)
        elif bp_type == "dssrLw":
            pairings, bp_markers, bp_map = getBasePairingEdgesDssrLw(dssrout, dssrids, points)
        elif bp_type == "none":
            pairings, bp_markers, bp_map = getBasePairingEdges(dssrout, dssrids, points)


        for item in pairings:
            # print(points[item[0]], points[item[1]])
            G.add_edge(item[0],item[1])
            edge_color_dict[(item[0],item[1])]= BP_EDGE_COLOR

    edges = getBackBoneEdges(ids, chids, dssrids, dssrout)
    for item in edges:
        G.add_edge(item[0],item[1])
        edge_color_dict[(item[0],item[1])]= BB_EDGE_COLOR

    d = deepcopy(G.edges)
    for item in d:
        if(points[item[0]][0] == points[item[1]][0] and points[item[0]][1] == points[item[1]][1]):
            #G.nodes[item[0]]['pos'] = (G.nodes[item[0]]['pos'] + np.random.random(2)*5)
            G.remove_edge(item[0],item[1])
    style = [style_dict[item] for item in G.edges]
    arrow = [arrow_dict[item]*extra['arrowsize'] for item in G.edges]
    nx.draw(G, nx.get_node_attributes(G, 'pos'), with_labels=False,
            node_size=160*magnification*extra['circlesize'],
            edgelist=[],
            node_color=colors, edgecolors='#000000')
    nx.draw_networkx_labels(G, nx.get_node_attributes(G, 'pos'), labels,
            font_size=(10+(magnification))*extra['circle_labelsize'], font_family='monospace',
            verticalalignment='center')
    nx.draw_networkx_edges(G, nx.get_node_attributes(G, 'pos'), edgelist=G.edges, style=style,
            edge_color=[edge_color_dict[item] for item in G.edges],
            arrowsize=arrow, width=1*magnification, arrowstyle='->')
    
    # If user does not want base pair annotations, turn these off by setting bp_type to something else

    #if bp_type == "dssr" and 'pairs' in dssrout.keys():
    #    for item in bp_markers:
    #        if item[2] == "k":
    #            item[2] = BP_MARKER_COLOR
    #        plt.scatter(item[0][0], item[0][1], marker=item[1], color=item[2], s = 80*magnification*marker_size,
    #                linewidth=1*magnification*marker_size, edgecolor=BP_MARKER_COLOR, label=item[3] )
    
    if bp_type == "saenger" and 'pairs' in dssrout.keys():
        for item in bp_markers:
            plt.text(item[0][0], item[0][1], item[1], color=BP_MARKER_COLOR, fontsize=10*np.sqrt(magnification)*marker_size)
    
    elif (bp_type == "dssrLw" or bp_type == "dssr") and 'pairs' in dssrout.keys():
        for item in bp_markers:
            if item[2] == "k":
                item[2] = BP_MARKER_COLOR
            if(item[1] == "ss" or item[1] == ">>" or item[1] == 'oo'): # just need one shape for these
                plt.scatter(item[0][0][0], item[0][0][1], marker=getCustomMarker(0, item), color=item[2], s = 80*magnification*marker_size,
                    linewidth=1*magnification*marker_size, edgecolor=BP_MARKER_COLOR, label=item[3] )
            else:
                # first shape!
                plt.scatter(item[0][1][0], item[0][1][1], marker=getCustomMarker(0, item), color=item[2], s = 80*magnification*marker_size,
                    linewidth=1*magnification*marker_size, edgecolor=BP_MARKER_COLOR, label=item[3] )
                
                # second shape!
                plt.scatter(item[0][2][0], item[0][2][1], marker=getCustomMarker(1, item), color=item[2], s = 80*magnification*marker_size,
                    linewidth=1*magnification*marker_size, edgecolor=BP_MARKER_COLOR, label=item[3] )        

    '''
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))True
    plt.legend(by_label.values(), by_label.keys(), title="non-WC bps")
    '''
    if bool(int(extra['showNumberLabels'])): # if user wants to show number labels
        poses, texts = getResNumPoints(points, ids, G, separation=float(extra['numberSeparation']))
        for i in range(len(poses)):
            plt.text(poses[i][0],poses[i][1], texts[i], color='saddlebrown',  fontsize=(10+(magnification))*float(extra['numberSize']))
    
    
    bounds = np.max(points, axis=0), np.min(points, axis=0)
    
    special_points = np.array([[bounds[0][0], bounds[1][1]], [bounds[1][0], bounds[0][1]]])
    center  =np.mean(special_points, axis=0)
    #d = np.linalg.norm((special_points[0] - special_points[1]))/2.828
    dir = (special_points[0] - center)
    dir2 = (special_points[1] - center)
    toplot = special_points[0] + (10)*dir/np.linalg.norm(dir)
    toplot2 = special_points[1] + (9)*dir2/np.linalg.norm(dir2)
    
    #toplot2 = center + (d+15)*dir2/np.linalg.norm(dir)

    plt.scatter(toplot[0], toplot[1], color="w",s=0)
    plt.scatter(toplot2[0], toplot2[1], color="w",s=0)

    plt.tight_layout()
    plt.gca().set_aspect('equal')
    #plt.savefig('{}/{}/{}{}{}.png'.format(MEDIA_PATH,FIG_PATH,prefix,time_string, int(extra['counter'])))
    #plt.savefig('{}/{}/{}{}{}.svg'.format(MEDIA_PATH,FIG_PATH,prefix,time_string, int(extra['counter'])))


    plt.close()


    # SAVE JSON of pertinent information to call regenerate labels!
    # maybe also print time string

    # for item in bp_map.keys():
    #     plt.scatter(0,0,marker=bp_map[item], color='k', label = 'c'+item, linewidth=1,
    #             edgecolor='k')
    #     plt.scatter(0,0,marker=bp_map[item], color='w', label = 't'+item, edgecolor='k',
    #             linewidth=1)
    # plt.legend()
    # plt.savefig('{}/{}/legend.svg'.format(MEDIA_PATH,FIG_PATH))

    # return '{}/{}{}{}.png'.format(FIG_PATH,prefix,time_string,rotation_string)
    return '{}/{}{}{}.svg'.format(FIG_PATH,prefix,time_string,int(extra['counter'])),'{}/{}{}{}.png'.format(FIG_PATH,prefix,time_string,int(extra['counter'])), "\n".join(log)


