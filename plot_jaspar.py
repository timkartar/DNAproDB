from pyjaspar import jaspardb
import numpy as np
from scipy.stats import entropy
import pandas as pd
import matplotlib.pyplot as plt

def makeLogo(gt, ax):
    import logomaker as lm
    
    e = entropy(gt,[0.25,0.25,0.25,0.25], axis=1, base=2)
    gt = gt * e[:,np.newaxis]
    gt_dict = {'A': [],
        'C':[],
        'G':[],
        'T':[]
        }
    for item in gt:
        gt_dict['A'].append(item[0])
        gt_dict['C'].append(item[1])
        gt_dict['G'].append(item[2])
        gt_dict['T'].append(item[3])

    gt_df = pd.DataFrame(gt_dict)
    colors = {'A': 'green',
        'C': 'blue',
        'G': '#FBA922',#'black',
        'T': 'red'
        }
    logo = lm.Logo(gt_df, color_scheme=colors, ax=ax, show_spines=False, alpha=0.7, fade_probabilities=False, font_name='DejaVu Sans')
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=15)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=15)
    ax.set_xlabel("Positions", fontsize=15)
    ax.set_ylabel("bits", fontsize=15)
    return logo

jdb_obj = jaspardb(release='JASPAR2024')
motifs = jdb_obj.fetch_motifs(all_versions=False)

l = []
unique_acs = set()
for motif in motifs:
    #print(motif.matrix_id)
    #print(motif.pwm)
    motif_string = str(motif).split("\n")

    for line in motif_string:
        if "Accession" in line:
            uniprots = (line.split("\t")[1].replace("[","").replace("]",""))

    if uniprots == "''":
        continue
    for item in uniprots.replace("'","").replace(" ","").split(","):
        #l.append(motif.matrix_id + "\t" + uniprots.replace("'","").replace(" ",""))
        l.append(item + "\t" + motif.matrix_id)
    for item in uniprots.replace("'","").replace(" ","").split(","):
        unique_acs.add(item)
    fig, ax = plt.subplots()
    pwm = np.array(list(motif.pwm.values())).T 
    makeLogo(pwm, ax)
    plt.savefig("./JASPAR/logos/{}.svg".format(motif.matrix_id))
    plt.close()

out = open("./JASPAR/uniprot_to_jaspar_map.tsv","w")
out.write("\n".join(l))
out.close()

