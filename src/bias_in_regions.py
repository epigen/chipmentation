
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
DF = pd.DataFrame(columns=["sample", "signal", "average", "x", "type"])

for region in ["K562_10M_CM_CTCF_nan_nan_0_0_hg19",
               "K562_10M_CM_GATA1_nan_nan_0_0_hg19",
               "K562_10M_CM_PU1_nan_nan_0_0_hg19",
               "K562_10M_CM_REST_nan_nan_0_0_hg19",
               "tss"]:
    print(region)
    if region != "tss":
        bed = pybedtools.BedTool(
            "/media/afr/cemm-backup/chipmentation/data/peaks/" +
            "{0}/{0}_peaks.motifCentered.bed".format(region))
    else:
        bed = pybedtools.BedTool("/media/afr/cemm-backup/hg19.cage_peak_coord_robust.TATA_Annotated.bed")
        bed = bed.slop(b=200, genome="hg19")

    # Get nucleotides
    seq = bed.sequence(s=True, fi=fasta)
    seqs = open(seq.seqfn).read().split("\n")[1::2]  # get fasta sequences

    # Match strings with PWM
    scores = np.array(range(len(seqs[0]) - 20))
    for sequence in seqs:
        result = MOODS.search(sequence, [PWM], 30)
        scores += np.array([j for i, j in result[0]])
    D[region] = scores / float(len(seqs))
    if region != "tss":
        DFF = pd.DataFrame([scores / float(len(seqs)), range(1980)]).T
    else:
        DFF = pd.DataFrame([scores / float(len(seqs)), range(380)]).T
    DFF.columns = ["average", "x"]
    DFF["sample"] = region
    DFF["signal"] = "pwm"
    DFF["type"] = "raw"
    DF = DF.append(DFF).reset_index(drop=True)

# save
pickle.dump(D, open("/media/afr/cemm-backup/chipmentation/tf_bias.pickle", "wb"), protocol=pickle.HIGHEST_PROTOCOL)
D = pickle.load(open("/media/afr/cemm-backup/chipmentation/tf_bias.pickle", "r"))

pickle.dump(DF, open("/media/afr/cemm-backup/chipmentation/tf_bias.df.pickle", "wb"), protocol=pickle.HIGHEST_PROTOCOL)
DF = pickle.load(open("/media/afr/cemm-backup/chipmentation/tf_bias.df.pickle", "r"))

# Plot
for region, data in D.items():
    plt.plot(data[1000 - 380:1000 + 380], label=region)  # plot center 800 bp
plt.legend()
plt.show()


# PWMs example
from Bio import Seq, Alphabet, motifs
a = Alphabet.IUPAC.IUPACUnambiguousDNA
seqs = [
    Seq.Seq("ACTG", alphabet=a),
    Seq.Seq("ACTG", alphabet=a),
    Seq.Seq("ACTG", alphabet=a),
    Seq.Seq("ACTG", alphabet=a),
    Seq.Seq("ACTC", alphabet=a),
    Seq.Seq("ACTT", alphabet=a),
    Seq.Seq("ACTC", alphabet=a),
    Seq.Seq("ACTT", alphabet=a),
    Seq.Seq("ACTA", alphabet=a)
]
m = motifs.create(seqs, alphabet=Alphabet.IUPAC.IUPACUnambiguousDNA)
m.pwm
m.weblogo("logo.png")
