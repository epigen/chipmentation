
"""

This does:

- Get sets of genomic regions
- Get the nucleotide sequences from those regions
- Scan nucleotide strings for a PWM and score every basepair (using MOODS)
- Sum/Average overal profiles per position
- Plot

"""

import numpy as np
import pandas as pd
import cPickle as pickle
import pybedtools
import MOODS
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("whitegrid")

# Get genome fasta file
fasta = pybedtools.BedTool("/media/afr/cemm-backup/reference/hg19/forBowtie2/hg19.fa")


# Get PWM
PWM = [
    [0.83593967562, 0.9612721605805, 0.924068101918, 1.06685938079, 0.797339497119, 1.18221399883, 0.8166032303705, 0.731357810088, 0.8009085680465, 0.6413220540835, 1.24263470188, 1.27537587614, 1.093537759345, 1.270194775255, 0.607808747652, 1.1045452182715, 0.830270118282, 0.9203172509935, 0.7788696601795, 0.820042458443, 1.071536938715],
    [1.058289907845, 0.9590219891645, 1.307638114905, 0.9474308687955, 1.625234242265, 0.703349925951, 1.126706911395, 1.15491741682, 1.29497649068, 1.446820024825, 0.669382028924, 0.673327304706, 0.85726490824, 0.842842178544, 1.703477693555, 0.883047249808, 0.911871662974, 1.071061400505, 1.12081896184, 1.356389737925, 1.075154356595],
    [1.075154356595, 1.356389737925, 1.12081896184, 1.071061400505, 0.911871662974, 0.883047249808, 1.703477693555, 0.842842178544, 0.85726490824, 0.673327304706, 0.669787971405, 1.446820024825, 1.29497649068, 1.15491741682, 1.126706911395, 0.703349925951, 1.625234242265, 0.9474308687955, 1.307638114905, 0.9590219891645, 1.058289907845],
    [1.071536938715, 0.820042458443, 0.7788696601795, 0.9203172509935, 0.830270118282, 1.1045452182715, 0.607808747652, 1.270194775255, 1.093537759345, 1.27537587614, 1.21553439588, 0.6413220540835, 0.8009085680465, 0.731357810088, 0.8166032303705, 1.18221399883, 0.797339497119, 1.06685938079, 0.924068101918, 0.9612721605805, 0.83593967562]
]

# Calculate matches to PWM
D = dict()
for region in ["ctcf", "gata1", "pu1", "rest", "tss"]:
    if region != "tss":
        bed = pybedtools.BedTool(
            "/media/afr/cemm-backup/chipmentation/data/peaks/" +
            "K562_10M_CM_{0}_nan_nan_0_0_hg19/K562_10M_CM_{0}_nan_nan_0_0_hg19_peaks.motifCentered.bed".format(region.upper()))
    else:
        bed = pybedtools.BedTool("/media/afr/cemm-backup/hg19.cage_peak_coord_robust.TATA_Annotated.bed")

    # Get nucleotides
    seq = bed.sequence(s=True, fi=fasta)
    seqs = open(seq.seqfn).read().split("\n")[1::2]  # get fasta sequences

    # Match strings with PWM
    scores = list()
    for sequence in seqs:
        result = MOODS.search(sequence, [PWM], 30)
        scores.append([j for i, j in result[0]])
    D[region] = pd.DataFrame(scores)

# average
d = map(lambda x: x.apply(np.mean, axis=0), D.values())
d = d.drop(0)
d = d.T
d.columns = ['rest', 'ctcf', 'pu1', 'gata1']
pickle.dump(d, open("tf_bias.pickle", "wb"), protocol=pickle.HIGHEST_PROTOCOL)

# Plot
d.plot()
