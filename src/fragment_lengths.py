import os
import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

cmd = """
getFragments() {
    NAME=$1
    samtools view /fhgfs/groups/lab_bock/shared/projects/chipmentation/data/mapped/${NAME}.trimmed.bowtie2.bam | \
    python /home/arendeiro/fragmentLengths.py $NAME
}
export -f getFragments

SAMPLES=(K562_10M_ATAC_H3K4ME1_nan_PE_1_1_hg19
    K562_10M_ATAC_PU1_nan_PE_1_1_hg19
    K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19
    K562_10M_CM_PU1_nan_PE_1_1_hg19
    K562_500K_ATAC_H3K4ME3_nan_01ULTN5_PE_1_1_hg19
    K562_500K_CM_H3K4ME3_nan_02ULTN5_PE_1_1_hg19
    K562_500K_CM_H3K4ME3_nan_05ULTN5_PE_1_1_hg19
    K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19
    K562_500K_CM_H3K4ME3_nan_5ULTN5_PE_1_1_hg19
    K562_500K_CM_IGG_nan_1ULTN5_PE_1_1_hg19)
parallel getFragments ::: ${SAMPLES[*]}
"""

os.system(cmd)

samples = [
    "K562_10M_ATAC_H3K4ME1_nan_PE_1_1_hg19", "K562_10M_ATAC_PU1_nan_PE_1_1_hg19",
    "K562_10M_CM_H3K4ME1_nan_PE_1_1_hg19", "K562_10M_CM_PU1_nan_PE_1_1_hg19",
    "K562_500K_ATAC_H3K4ME3_nan_01ULTN5_PE_1_1_hg19", "K562_500K_CM_H3K4ME3_nan_02ULTN5_PE_1_1_hg19",
    "K562_500K_CM_H3K4ME3_nan_05ULTN5_PE_1_1_hg19", "K562_500K_CM_H3K4ME3_nan_1ULTN5_PE_1_1_hg19",
    "K562_500K_CM_H3K4ME3_nan_5ULTN5_PE_1_1_hg19", "K562_500K_CM_IGG_nan_1ULTN5_PE_1_1_hg19"
]

for i in range(len(samples)):
    name = samples[i]
    dist = pickle.load(open(name + ".pickle", "r"))
    dist.pop(0)

    plt.plot(dist.keys(), dist.values(), '-', label=name)

    plt.legend(loc='upper right', shadow=False)
    plt.xlim(0, 300)
    plt.savefig(name + "_fragmentDistribution.pdf")
    plt.close()

    # Get percentiles
    import numpy as np
    l = [size for size in dist.keys() for _ in range(dist[size])]
    print(np.percentile(l, [20, 40, 60, 80]))
