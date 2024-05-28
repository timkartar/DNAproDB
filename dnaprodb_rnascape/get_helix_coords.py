from Bio.PDB import MMCIFParser
import Bio
import json
import re, os, sys
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

bp_width = 2.5
TIGHT_HELIX_THRES=0.5
def get_cetroid(res):
    atoms = []
    for atom in res:
        atoms.append(atom.coord)
    return np.mean(atoms, axis=0)

def untillastnumber(l):
    revl = l[::-1]
    
    for idx in range(len(revl)):
        if not revl[idx].isdigit():
            break
    return revl[idx:][::-1]

def process_resid(nt1, model):
    spl1 = nt1.split(".")
    if len(spl1) < 6:
        if "^" in spl1[1]:
            tmp = spl1[1].split("^")[0]
            icode = spl1[1].split("^")[1]
            spl1[1] = tmp
        else:
            icode = ' ' 


        res1 = re.sub("\D",",",spl1[1]).split(",")[-1]
        

        rest1_n = untillastnumber(spl1[1])
        rest1 = spl1[1]
        
        if (rest1_n not in list('AUGC') + ['DA','DC','DG','DT']):
            het1 = 'H_' + rest1_n
        else:
            het1=' '
        chid = spl1[0]
    else:
        if spl1[5] == '':
            icode = ' '
        else:
            icode = spl1[5]
        rest1_n = spl1[3]
        rest1 = spl1[3] + spl1[4]
        res1 = spl1[4]
        if (rest1_n not in list('AUGC') + ['DA','DC','DG','DT']):
            het1 = 'H_' + rest1_n
        else:
            het1=' '
        chid = spl1[2]

    def check_len5_residue(spl1, id1, rest1, chid1):
        try:
            model[chid1][id1]
        except:
            for res in model[chid1]:
                if id1[1] == res.get_id()[1]:
                    spl1[3] = res.get_resname()
                    id1 = (res.get_id()[0], id1[1], id1[2])
                    rest1 = res.get_resname()
        return spl1, id1, rest1, chid1

    return check_len5_residue(spl1, (het1.strip("/"), int(res1), icode), rest1_n.strip("/"), chid)
    #return spl1, (het1.strip("/"), int(res1), icode), rest1_n.strip("/"), chid
    #return spl1, (het1.strip("/"), int(res1), icode), rest1.strip("/"), chid


def get_helix_coords(dssrout, model):
    helices = []
    if 'helices' not in dssrout:
        return None
    for helix in dssrout['helices']:
        helix_data = {
            'coords' : [],
            'markers' : [],
            'pairs' : [],
            'ids':[],
            'chids':[],
            'dssrids':[]
        }
        for pair in helix['pairs']:
            nt1 = pair['nt1']
            nt2 = pair['nt2']
            
            spl1, id1, rest1, chid1 = process_resid(nt1, model)
            spl2, id2, rest2, chid2 = process_resid(nt2, model)
            ntc1 = get_cetroid(model[chid1][id1])
            ntc2 = get_cetroid(model[chid2][id2])
            #    sys.exit()
            

            helix_data['coords'].append(ntc1)
            helix_data['coords'].append(ntc2)
            helix_data['markers'].append("${}$".format(rest1))
            helix_data['markers'].append("${}$".format(rest2))
            helix_data['pairs'].append((rest1, rest2))
            helix_data['ids'].append(id1)
            helix_data['ids'].append(id2)
            helix_data['chids'].append(chid1)
            helix_data['chids'].append(chid2)
            helix_data['dssrids'].append(nt1)
            helix_data['dssrids'].append(nt2)
             
        helices.append(helix_data)

    coords = []
    markers = []
    ids = []
    chids = []
    dssrids = []

    for item in helices:
        coords += item['coords']
        markers += item['markers']
        ids += item['ids']
        chids += item['chids']
        dssrids += item['dssrids']

    pca = PCA(n_components=2)        
    pcs = pca.fit_transform(np.array(coords))

    '''
    for i in range(len(markers)):
        plt.scatter(pcs[i,0], pcs[i,1], marker=markers[i], edgecolors='none', color='black', s=100)

    #plt.show()
    #plt.close()
    '''
    def transformHelix(helix, style="heur"):
        pc2 = pca.transform(helix['coords'])
        if style=="heur":
            first = (pc2[0] + pc2[1])/2
            last = (pc2[-1] + pc2[-2])/2
            n = len(pc2)//2
            def comp(first, last, n):
                pc1=[first]
                for i in range(n-1):
                    p = first + (i+1)*(last - first)/(n-1)
                    pc1.append(p)
                return pc1
            pc1 = comp(first, last, n)
            mean = pc2.mean(axis=0)
            ax_mean = np.array(pc1).mean(axis=0)
            if(np.linalg.norm(mean - ax_mean)) > 10: ##approximation is crude, need to divide into two
                m = n - n%2
                mid1 = (pc2[m] + pc2[m+1])/2  
                mid2 = (pc2[m+2] + pc2[m+3])/2 
                pc1 = comp(first, mid1, n//2) + comp(mid2, last, n-n//2)
            
            mean = pc2.mean(axis=0)
            ax_mean = np.array(pc1).mean(axis=0)
        return pc1

    helix_axes = []
    for item in helices:
        helix_axes += [transformHelix(item)]
    def processHelixAxes(helix_axes):
        axes = []

        toalign = []
        for axis in helix_axes:
            last = axis[-1]
            if toalign == []:
                toalign.append([axis])
                continue
            else:
                prev = toalign[-1][-1]
                prelast = prev[-1]
                predir = prev[-1] - prev[-2]
                predir = predir/np.linalg.norm(predir)
                first = axis[0]
                currdir = axis[1] - axis[0]
                currdir = currdir/np.linalg.norm(currdir)
                deviation = np.arccos(np.dot(predir,currdir))

                if np.linalg.norm(first - prelast) < 20: #(check if within certain distance)
                    if deviation < np.pi/5: # check if axis deviation less than pi/6 (30 degrees)
                        toalign[-1].append(axis)
                        # toalign.append([axis])
                    else:
                        toalign.append([axis])
                else:
                    toalign.append([axis])
            axes += axis
        aligned = []
        for item in toalign:
            if len(item) == 1:
                aligned+=item
                continue
            alll = []
            for ax in item:
                alll+=ax
            first = alll[0]
            last = alll[-1]
            n = len(alll)
            toappend = [first]
            for i in range(n-1):
                p = first + (i+1)*(last - first)/(n-1)
                toappend.append(p)
            aligned += [toappend]
        return aligned
    
    #TODO uncomment to merge helices
    helix_axes= processHelixAxes(helix_axes)
    
    ## extend tight helices ##
    #print(helix_axes[0])
    
    for item in helix_axes:
        v = item[1] - item[0]
        vl = (np.linalg.norm(v))
        if (vl / bp_width) < TIGHT_HELIX_THRES:
            n = len(item)
            for  i in range(1,n):
                item[i]  = item[i-1] + v * bp_width/vl
    
    #print(helix_axes[0])
    #sys.exit()
    #helix_axes = np.array(helix_axis)
    '''
    for helix_axis in helix_axes:
        for i in range(len(helix_axis)):
            d = np.array(helix_axis)
            plt.scatter(d[i,0], d[i,1], edgecolors='none', color='black', s=100)
    plt.show()
    plt.close()
    '''
    points = []

    for helix_axis in helix_axes:   
        for i in range(len(helix_axis)-1):
            local_axis = (helix_axis[i+1] - helix_axis[i])
            local_axis_3ded = np.array(local_axis.tolist() + [0])
            zaxis = np.array([0,0,1])
            local_perp = np.cross(local_axis_3ded, zaxis)[:2]
            local_perp = local_perp/np.linalg.norm(local_perp)
            p1 = helix_axis[i] + local_perp*bp_width
            p2 = helix_axis[i] - local_perp*bp_width
            points += [p2, p1]


        last_point1 = points[-2] + helix_axis[-1] - helix_axis[-2]
        last_point2 = points[-1] + helix_axis[-1] - helix_axis[-2]

        points += [last_point1, last_point2]
    points = np.array(points)
    
    '''
    for i in range(len(markers)):
        plt.scatter(points[i,0], points[i,1], marker=markers[i], edgecolors='none', color='black', s=100)
    plt.show()
    '''
    return points, ids, markers, chids, dssrids

if __name__ == "__main__":
    prefix = sys.argv[1]
    parser = MMCIFParser()
    model = parser.get_structure(prefix,"./vn/{}-assembly1.cif".format(prefix))[0]

    with open("./json/{}.json".format(prefix),"r") as f:
        dssrout = json.load(f)

    get_helix_coords(dssrout, model)



