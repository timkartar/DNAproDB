from Bio.PDB import MMCIFParser, PDBParser
import Bio
import json
import re, os, sys
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from .get_helix_coords import get_helix_coords, process_resid
#from plot import Plot
from math import cos, sin
import re 
from sklearn.neighbors import KDTree

DSSR_PATH = ''
CIF_PATH=''
FIG_PATH=''

dssrout = None
tree=None
conditional_bulging = True
Model = None
def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def circularLayout(n, m, d, theta, factor=False):
    poses = []
    rot = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
    for i in range(n):
        d = np.dot(rot, d)

        if factor:
            poses.append(m+d*np.sqrt(n))
        else:
            poses.append(m+d*1.5)

    return poses

def perp(v):
    v = v.tolist() + [0]
    z = [0,0,1]
    p = np.cross(v,z)[:2]
    p = p/(np.linalg.norm(p)+0.00000001)
    return p

def updateLoopPoints(start_pos, end_pos, val, helix_coords, factor=False):
    m = (start_pos + end_pos)/2
    poses = []
    v = (end_pos - start_pos)
    
    p = perp(v)
    ## generate points circularly using m,v,p
    n = len(val)
    #m = m + p*np.linalg.norm(v)*0.2
    d = start_pos - m
    
    
    theta = np.pi/(n+1)
    poses = circularLayout(n,m,d,theta, factor)
    neg_poses = circularLayout(n,m,d,-1*theta, factor)
    helix_coords_m = np.mean(helix_coords, axis=0)
    poses_m = np.mean(poses, axis=0)
    neg_poses_m = np.mean(neg_poses, axis=0)
    
    dis = np.linalg.norm(poses_m - helix_coords_m)
    n_dis = np.linalg.norm(neg_poses_m - helix_coords_m)
    
    if len(poses) == 0:
        return poses
    
    l = tree.query_radius(poses, r=5, count_only=True).sum()
    n_l = tree.query_radius(neg_poses, r=5, count_only=True).sum()
    
    if l < n_l:
        return poses
    else:
        return neg_poses

def generate_coords(helix_coords, helix_ids, dic, helix_dssrids, dssrout):
    pairs = []
    for item in dssrout['pairs']:
        pairs.append((item['nt1'],item['nt2']))

    
    positions = []
    markers = []
    ids = []
    chids = []
    dssrids = []
    for key, val in dic.items():
        start = key[0]
        end = key[1]
        if(start == None or end ==None):
            continue
        start_pos = helix_coords[key[0]]
        end_pos = helix_coords[key[1]]

        t_poses = []
        t_markers = []
        t_ids = []
        t_chids = []
        t_dssrids = []

        l = len(val)
        
        def check_multichain(val):
            chains = [item[2] for item in val]
            chains = list(set(chains))
            if len(chains) > 1:
                return True
            else:
                return False
        multi_chain = check_multichain(val)
        c1 = []
        c2 = []
        chains = []
        for item in val:
            if item[2] not in chains:
                chains.append(item[2])
        #chains = list(chains)
        #print(chains)
        for item in val:
            if item[2] == chains[0]:
                c1.append(item)
            else:
                c2.append(item)
        #print(c1, c2)

        for i in range(len(val)):
            item = val[i]
            t_ids.append(item[0]) 
            t_markers.append("${}$".format(item[1]))
            v = (end_pos - start_pos)
            pos = start_pos + v*(i+1)/(l+1)
            
            if multi_chain:
                if item in c1:
                    try:
                        c1_p11 = helix_coords[start - 2]
                    except:
                        c1_p11 = helix_coords[start + 2]
                    try:
                        c1_p12 = helix_coords[start + 2]
                    except:
                        c1_p12 = helix_coords[start - 2]
                    c1_p2 = helix_coords[start]
                    if np.linalg.norm(c1_p2 - c1_p11) < np.linalg.norm(c1_p2 - c1_p12):
                        c1_p1 = c1_p11
                    else:
                        c1_p1 = c1_p12

                    pos = c1_p2 + (c1.index(item)+1)*(c1_p2 - c1_p1)
                    if pos in helix_coords:
                        pos = c1_p2 - (c1.index(item)+1)*(c1_p2 - c1_p1)
                if item in c2:
                    c2_p1 = helix_coords[end]
                    try:
                        c2_p21 = helix_coords[end-2]
                    except:
                        c2_p21 = helix_coords[end+2]
                    try:    
                        c2_p22 = helix_coords[end+2]
                    except:
                        c2_p22 = helix_coords[end-2]

                    if np.linalg.norm(c2_p1 - c2_p21) < np.linalg.norm(c2_p1 - c2_p22):
                        c2_p2 = c2_p21
                    else:
                        c2_p2 = c2_p22
                    
                    #for i in range(len(c2)):
                    pos = c2_p1 + (len(c2) - (c2.index(item)))*(c2_p1 - c2_p2)
                    if pos in helix_coords:
                        pos = c2_p1 - (len(c2) - (c2.index(item)))*(c2_p1 - c2_p2)
                
            if pos in helix_coords:
                p = perp(v)
                tpos = pos + p*2 + v/(2*(l+1)) ## perpeindicular and out shift for overlap
                neg_tpos = pos - p*2 + v/(2*(l+1)) ## perpeindicular and out shift for overlap
                l = tree.query_radius([tpos], r=5, count_only=True)
                n_l = tree.query_radius([neg_tpos], r=5, count_only=True)
                #if dis > n_dis:
                if l < n_l:
                    pos = tpos
                else:
                    pos = neg_tpos

           
            t_poses.append(pos)
            t_chids.append(item[2])
            t_dssrids.append(item[3])
        v = helix_coords[start] - helix_coords[end]
        
        if not multi_chain or (multi_chain and np.linalg.norm(v) < 5):
            if not conditional_bulging:
                t_poses = updateLoopPoints(start_pos, end_pos, val, helix_coords, factor=False)
            elif (helix_dssrids[start], helix_dssrids[end]) in pairs or (helix_dssrids[end], helix_dssrids[start]) in pairs:
                t_poses = updateLoopPoints(start_pos, end_pos, val, helix_coords)
            elif np.linalg.norm(v) < 3: #threshold for bulging
                t_poses = updateLoopPoints(start_pos, end_pos, val, helix_coords, factor=True) 
            elif np.linalg.norm(v)/(len(val) +  0.0001) < 1.5: #threshold for bulging
                t_poses = updateLoopPoints(start_pos, end_pos, val, helix_coords, factor=True)

        #else:
        #    t_poses = updateLoopPoints(start_pos, end_pos, val, helix_coords) 
        positions+=(t_poses)
        markers+=(t_markers)
        ids+=(t_ids)
        chids+=(t_chids)
        dssrids+=(t_dssrids)
    
    return positions, markers, ids, chids, dssrids

def get_coords(nts, helix_ids, helix_coords, dssrids, dssrout):
    dic = {} #keys are like (start, end) values are like [(nt_id, rest1),...]

    start = None
    end = None
    l = []
    prev = False
    covered = []
    starters = []
    
    for item in dssrout['nts']:
        spl1, nt_id, rest1, chid =  process_resid(item['nt_id'], Model)  

        if nt_id not in helix_ids or item['nt_id'] not in dssrids:
            starters.append((nt_id, rest1, chid, item['nt_id']))
        else:
            break

    enders = []
    for item in dssrout['nts'][::-1]:
        spl1, nt_id, rest1, chid =  process_resid(item['nt_id'], Model)
        if nt_id not in helix_ids or item['nt_id'] not in dssrids:      
            enders.append((nt_id, rest1, chid, item['nt_id']))
        else:
            idx = dssrids.index(item['nt_id'])
            break
    
    dic[(0,0)] = starters
    dic[(idx, idx)] = enders
    
    curr_chain = None
    
    for item in dssrout['nts']:
        spl1, nt_id, rest1, chid =  process_resid(item['nt_id'], Model)

        if nt_id in helix_ids and item['nt_id'] in dssrids:
            if prev == False:
                start = dssrids.index(item['nt_id'])
            else:
                end = dssrids.index(item['nt_id'])
                prev = False
                dic[(start, end)] = l
                start = dssrids.index(item['nt_id'])
                covered += l
                l = []

        else:
            prev = True
            l.append((nt_id, rest1, chid, item['nt_id'])) 
    
    return generate_coords(helix_coords, helix_ids, dic, dssrids, dssrout)        

def orderData(points, markers, ids, chids, dssrids):
    unique_chains = np.unique(chids)
    d = {}
    for item in unique_chains:
        d[item] = []

    for i in range(len(ids)):
        d[chids[i]].append(ids[i][1])
    
    sorted_nice = []
    for k in d.keys():
        d[k] = np.sort(d[k])
        d[k] = ["{}{}".format(i,k) for i in d[k]]
        sorted_nice += d[k]

    resnumbers = []
    icodes = {}
    for i in range(len(ids)):
        resnumbers.append("{}{}".format(str(ids[i][1]), chids[i]))
        #full_resnumbers.append("{}{}".format(str(ids[i][1]) + ids[i][2].replace(' ',''), chids[i]))

    argsorted = []
    
    done_indices = []

    def find_all(l, k):
        ret = []
        for i in range(len(l)):
            if l[i] == k:
                ret.append(i)
        return ret

    for item in sorted_nice:
        indices = find_all(resnumbers, item)
        #idx = resnumbers.index(item)
        for idx in indices:
            if idx not in done_indices:
                break
            
        done_indices.append(idx)
        argsorted.append(idx)
    points = points[argsorted,:]
    markers = np.array(markers)[argsorted].tolist()
    chids = np.array(chids)[argsorted].tolist()
    dssrids = np.array(dssrids)[argsorted].tolist()
    #sys.exit()
    ids = np.array(ids)[argsorted].tolist()
    return points, markers, ids, chids, dssrids, d

def getTails(helix_dssrids, dssrids, chids, points):
    starters = {}
    enders = {}
    for item in np.unique(chids):
        starters[item] = []
        enders[item] = []
    
    starting = True
    for i in range(len(dssrids)-1):
        chid = chids[i]
        if dssrids[i] not in helix_dssrids:
            if(starting):
                starters[chid].append(i)
        else:
            starting = False
            if chids[i+1] != chid:
                starting=True

    
    ending = True
    rev_dssrids = dssrids[::-1]
    rev_chids = chids[::-1]
    for i in range(len(dssrids)-1):
        chid = rev_chids[i]
        if rev_dssrids[i] not in helix_dssrids:
            if(ending):
                enders[chid].append(len(dssrids) -1 - i)
        else:
            ending = False
            if rev_chids[i+1] != chid:
                ending=True
    
    #print([(item, starters[item]) for item in starters])
    #print([(item, enders[item]) for item in enders])
    for k in starters.keys():
        if len(starters[k]) == 0:
            continue
        ip1 = starters[k][-1] + 1
        ip2 = starters[k][-1] + 2
        v = points[ip1] - points[ip2]
        n = len(starters[k])
        for i in range(len(starters[k])):
            points[starters[k][n-i-1]] = points[ip1] + v*(i+1)
            #print(points[ip1],points[ip2], points[starters[k][n-i-1]])

    for k in enders.keys():
        if len(enders[k]) == 0:
            continue
        enders[k] = enders[k][::-1]
        

        ip1 = enders[k][0] - 1
        ip2 = enders[k][0] - 2

        v = points[ip1] - points[ip2]
        for i in range(len(enders[k])):
            points[enders[k][i]] = points[ip1] + v*(i+1)
    
    return starters, enders, points


#if __name__ == "__main__":
"""
cond_bulging: True = attempt to condense non-base-pairing nucleotides. False = always bulge them out as loop
prefix: everything before .cif or .pdb in the file
cif_file: cif or pdb file path
"""
def rnascape(prefix, structure, dssr, cond_bulging=True, mDSSR_PATH='', mFIG_PATH='' ):
    #DSSR_PATH = mDSSR_PATH
    FIG_PATH = mFIG_PATH
    global tree, dssrout, conditional_bulging, Model
    dssrout=dssr
    conditional_bulging = cond_bulging
    #prefix = sys.argv[1]
    #if os.path.exists("{}/{}.png".format(FIG_PATH, prefix)):
    #    sys.exit()

    # support both CIF and PDB files
    #if cif_file.endswith(".cif"):
    #    parser = MMCIFParser()
    #elif cif_file.endswith(".pdb"):
    #    parser = PDBParser()
    #model = parser.get_structure(prefix,"./vn/{}-assembly1.cif".format(prefix))[0]
    #model = parser.get_structure(prefix,cif_file)[0]
    model=structure #[0]
    Model = model
    
    figpath=''

    #with open(json_file,"r") as f:
    #    dssrout = json.load(f)

    helices = get_helix_coords(dssrout, model)
    if helices == None:
        from .get_helix_coords import get_cetroid
        coords = []
        markers = []
        ids = []
        chids = []
        dssrids = []
        for item in dssrout['nts']:
            spl1, nt_id, rest, chid = process_resid(item['nt_id'], Model)
            ntc = get_cetroid(model[chid][nt_id])
            coords.append(ntc)
            markers.append('${}$'.format(rest))
            ids.append(nt_id)
            chids.append(chid)
            dssrids.append(item['nt_id'])
        pca = PCA(n_components=2)
        points = pca.fit_transform(np.array(coords))
    #   V  Don’t know whether to include this line  V
    #   points, markers, ids, chids, dssrids, dic = orderData(points, markers, ids, chids, dssrids)
        

    else:
        helix_points, helix_ids, helix_markers, helix_chids, helix_dssrids = get_helix_coords(dssrout, model)
        tree=KDTree(helix_points)
        '''
        l = [helix_points, helix_ids, helix_markers, helix_chids, helix_dssrids]
        for item in l:
            print(len(item))

        for i in range(len(helix_markers)):
            plt.scatter(helix_points[i,0], helix_points[i,1], marker=helix_markers[i], edgecolors='none', color='black', s=200)
        plt.show()
        '''
        '''
        print(helix_dssrids,
                helix_points,
                helix_chids)
        '''
        
        rest_positions, rest_markers, rest_ids, rest_chids, rest_dssrids = get_coords(dssrout,
                helix_ids, helix_points, helix_dssrids, dssrout)
        
        points = np.array(helix_points.tolist() + rest_positions)
        
        markers = helix_markers + rest_markers
        
        ids = helix_ids + rest_ids
        
        chids = helix_chids + rest_chids
        dssrids = helix_dssrids + rest_dssrids
        
        points, markers, ids, chids, dssrids, dic = orderData(points, markers, ids, chids, dssrids)
        #points = updateLoopPoints(points, dssrids, dssrout)
        
        '''
        idx = (np.argsort(chids, kind='mergesort'))
        chids = np.array(chids)[idx].tolist()

        ids = np.array(ids)[idx].tolist()

        markers = np.array(markers)[idx].tolist()

        dssrids = np.array(dssrids)[idx].tolist()

        points= points[idx,:]
        '''

        starters, enders, points = getTails(helix_dssrids, dssrids, chids, points)
        #figpath = Plot(points, markers, ids, chids, dssrids, dssrout, prefix)
    return points, markers, ids, chids, dssrids, dssrout, prefix
