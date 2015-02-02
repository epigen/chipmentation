#!/usr/env python

from argparse import ArgumentParser
import os
import HTSeq
import cPickle as pickle
import multiprocessing
import parmap
from collections import Counter
import itertools


def distances(feature, bam, fragment_size, duplicates=True, strand_wise=True, permutate=False):
    """
    Gets pairwise distance between reads in a single interval. Returns dict with distance:count.
    If permutate=True, it will randomize the reads in each interval along it.

    feature=HTSeq.GenomicInterval object.
    bam=HTSeq.BAM_Reader object.
    fragment_size=int.
    duplicates=bool.
    strand_wise=bool.
    permutate=bool.
    """
    dists = Counter()
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX']
    
    # FIXME:
    #if feature.chrom not in chroms:
    #    continue
    # Fetch all alignments in feature window
    alnsInWindow = bam[feature]
    if permutate:
        # randomize each alignment's position in window
        alns = list()
        for aln in alnsInWindow:
            aln.iv.start_d = random.randrange(feature.start, feature.end)
            alns.append(aln)
        alnsInWindow = alns

    # Measure distance between reads in window, pairwisely
    for aln1, aln2 in itertools.combinations(alnsInWindow, 2):
        # check if duplicate
        if not duplicates and (aln1.pcr_or_optical_duplicate or aln2.pcr_or_optical_duplicate):
            continue
        # check if in same strand
        if not strand_wise and aln1.iv.strand != aln2.iv.strand:
            continue
        # adjust fragment to size
        aln1.iv.length = fragment_size
        aln2.iv.length = fragment_size

        # get position relative
        dist = abs(aln1.iv.start_d - aln2.iv.start_d)
        # add +1 to dict
        if strand_wise:
            if aln1.iv.strand == "+":
                dists[dist] += 1
            elif aln1.iv.strand == "-":
                dists[-dist] += 1
        else:
            dists[dist] += 1
    return dists


if __name__ == '__main__':
    parser = ArgumentParser(
        description = 'count_distances_parallel.py',
        usage       = 'python count_distances_parallel.py <directory> input_pickle '
    )

    ### Global options
    # positional arguments
    parser.add_argument(dest='input_pickle', type=str, help='Pickle file to load.')
    parser.add_argument(dest='output_pickle', type=str, help='Pickle file to save to.')
    parser.add_argument(dest='bam_file', type=str, help = 'Bam file.')
    
    # optional arguments
    parser.add_argument('--strand-wise', dest='strand_wise', action='store_true')
    parser.add_argument('--duplicates', dest='duplicates', action='store_true')
    parser.add_argument('--permute', dest='permute', action='store_true')
    parser.add_argument('--fragment-size', dest='fragment_size', type=int, default=1)

    args = parser.parse_args()

    ### Read pickle with windows
    windows = pickle.load(open(args.input_pickle, 'r')) # dict-like list of tuples
    # convert list of tuples to list
    windows = [tup[1] for tup in windows]

    ### Make bam object
    bam = HTSeq.BAM_Reader(os.path.abspath(args.bam_file))

    ### Process in parallel and serialize result
    # parallel process and reduce to Counter
    dists = reduce(lambda x, y: x + y, parmap.map(distances, windows, bam, args.fragment_size, duplicates=args.duplicates, strand_wise=args.strand_wise, permutate=args.permute))
    
    ### Serialize
    pickle.dump(dists, open(args.output_pickle, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
