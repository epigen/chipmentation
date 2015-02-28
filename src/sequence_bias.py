#!/usr/env python

"""
A more deep analysis of ChIPmentation data.

"""

from argparse import ArgumentParser
from collections import Counter
from matplotlib import pyplot as plt
import cPickle as pickle
import HTSeq
import numpy as np
import os
import sys
import random
import re
import seaborn as sns
import string
from Bio import SeqIO
import itertools

np.set_printoptions(linewidth=200)


def main():
    # Parse command-line arguments
    parser = ArgumentParser()

    # Global options
    # positional arguments
    parser.add_argument(dest='data_dir', type=str, help='Directory to save data to.')
    parser.add_argument(dest='results_dir', type=str, help='Directory to save data to.')
    parser.add_argument(dest='plots_dir', type=str, help='Directory to save plots to.')
    parser.add_argument('bam_files', nargs='*', help='Bam files')
    # optional arguments
    parser.add_argument('--duplicates', dest='duplicates', action='store_true')
    parser.add_argument('-k', dest='k', type=int, default=5)
    parser.add_argument('--genome', dest='genome', type=str, default='hg19')
    args = parser.parse_args()
    args = parser.parse_args(["projects/chipmentation/results/periodicity-data",
                              "projects/chipmentation/results/periodicity",
                              "projects/chipmentation/results/plots/periodicity",
                              "data/human/chipmentation/mapped/merged/DNase_UWashington_K562_mergedReplicates.bam",
                              "data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_CM.bam",
                              "data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_ChIP.bam",
                              "data/human/chipmentation/mapped/merged/H3K27me3_K562_10k500k_CM.bam",
                              "data/human/chipmentation/mapped/merged/H3K27me3_K562_500k_ChIP.bam",
                              "data/human/chipmentation/mapped/merged/PU1_K562_10mio_CM.bam",
                              "data/human/chipmentation/mapped/merged/CTCF_K562_10mio_CM.bam",
                              "/fhgfs/groups/lab_bock/arendeiro/projects/atac-seq/data/mapped/ASP14_50k_ATAC-seq_nan_nan_DOX_ATAC10-8_0_0.trimmed.bowtie2.shifted.dups.bam",
                              "/fhgfs/groups/lab_bock/arendeiro/projects/atac-seq/data/mapped/ASP14_50k_ATAC-seq_nan_nan_untreated_ATAC10-7_0_0.trimmed.bowtie2.shifted.dups.bam"
                              ]
                             )

    # TODO: mkdirs
    args.data_dir = os.path.abspath(args.data_dir)
    args.results_dir = os.path.abspath(args.results_dir)
    args.plots_dir = os.path.abspath(args.plots_dir)

    counts = dict()

    # Get sample names
    samples = {re.sub("\.bam", "", os.path.basename(sampleFile)): os.path.abspath(sampleFile) for sampleFile in args.bam_files}

    # for each sample, count kmers
    for sampleName, sampleFile in samples.items():
        print(sampleName)
        bam = HTSeq.BAM_Reader(sampleFile)
        counts[sampleName] = countKmers(bam, args.k)

    pickle.dump(counts, open(os.path.join(args.results_dir, "kmerCounts.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    # Get genome background
    kmers = ["".join(i) for i in itertools.product('ACTGN', repeat=5)]
    countsGenome = Counter()

    for seq_record in SeqIO.parse("/fhgfs/prod/ngs_resources/genomes/hg19/hg19.fa", "fasta"):
        print(seq_record.name)
        for kmer in kmers:
            countsGenome[kmer] = seq_record.seq.count(kmer)

    pickle.dump(countsGenome, open(os.path.join(args.results_dir, "kmerCountsGenome.pickle"), "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    # plot


def countKmers(bam, k, start=0, duplicates=False):
    """
    Count kmers in reads in Bam file.

    bam=HTSeq.BAM_Reader object.
    k=int.
    duplicates=bool.
    """
    # TODO:
    # implement strand-wise
    kmers = Counter()

    # For all alignments
    # Get k bases from 5 prime end
    # Add to Counter
    for aln in bam:
        # check if duplicate
        if not duplicates and aln.pcr_or_optical_duplicate:
            continue
        # get sequence kmer
        kmer = aln.read.seq[start:k + start]
        # append to counter
        kmers[kmer] += 1
    return kmers


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(0)
