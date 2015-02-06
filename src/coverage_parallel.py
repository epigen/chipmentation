#!/usr/env python

from argparse import ArgumentParser
import os
import HTSeq
import cPickle as pickle
import multiprocessing
import parmap
import itertools
import numpy as np


def coverageInWindow(feature, bam, fragmentsize, orientation, duplicates, strand_wise, permutate):
    """
    Gets read coverage in a single genomic interval. Returns 1D np.array if strand_wise=False, 2D if True.

    feature=HTSeq.GenomicInterval object.
    bam=HTSeq.BAM_Reader object - Must be sorted and indexed with .bai file!
    fragmentsize=int.
    stranded=bool.
    duplicates=bool.
    """
    # Loop through TSSs, get coverage, append to dict
    # FIXME:
    # select only features in chroms
    chroms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX']
    
    # Initialize empty array for this feature
    if not strand_wise:
        profile = np.zeros(feature.length, dtype=np.float64)
    else:
        profile = np.zeros((2, feature.length), dtype=np.float64)

    # Check if feature is in bam index 
    if feature.chrom not in chroms or feature.chrom == "chrM":
        return profile

    # Fetch all alignments in feature window
    alnsInWindow = bam[feature]
    if permutate:
        # randomize each alignment's position in window
        alns = list()
        for aln in alnsInWindow:
            aln.iv.start_d = random.randrange(feature.start, feature.end)
            alns.append(aln)
        alnsInWindow = alns

    for aln in alnsInWindow:
        # check if duplicate
        if not duplicates and aln.pcr_or_optical_duplicate:
            continue
        aln.iv.length = fragmentsize # adjust to size

        # get position in relative to window
        if orientation:
            if feature.strand == "+" or feature.strand == ".":
                start_in_window = aln.iv.start - feature.start - 1
                end_in_window   = aln.iv.end   - feature.start - 1
            else:
                start_in_window = feature.length - abs(feature.start - aln.iv.end) - 1
                end_in_window   = feature.length - abs(feature.start - aln.iv.start) - 1
        else:
            start_in_window = aln.iv.start - feature.start - 1
            end_in_window   = aln.iv.end   - feature.start - 1
        
        # check fragment is within window; this is because of fragmentsize adjustment
        if start_in_window <= 0 or end_in_window > feature.length:
            continue

        # add +1 to all positions overlapped by read within window
        if not strand_wise:
            profile[start_in_window : end_in_window] += 1
        else:
            if aln.iv.strand == "+":
                profile[0][start_in_window : end_in_window] += 1
            else:
                profile[1][start_in_window : end_in_window] += 1
    return profile


if __name__ == '__main__':
    parser = ArgumentParser(
        description = 'coverage_parallel.py',
        usage       = 'python coverage_parallel.py [options] input_pickle output_pickle bam_file'
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
    parser.add_argument('--orientation', dest='orientation', action='store_true')
    parser.add_argument('--fragment-size', dest='fragment_size', type=int, default=1)

    args = parser.parse_args()

    ### Read pickle with windows
    windows = pickle.load(open(args.input_pickle, 'r')) # dict-like list of tuples

    ### Make bam object
    bam = HTSeq.BAM_Reader(os.path.abspath(args.bam_file))

    ### Process in parallel and serialize result
    # parallel process and reduce to Counter
    coverage = dict(
        zip(
            [tup[0] for tup in windows],
            parmap.map(
                coverageInWindow,
                [tup[1] for tup in windows], #similar to windows.values()
                bam,
                args.fragment_size,
                args.orientation,
                args.duplicates,
                args.strand_wise,
                args.permute
            )
        )
    )
    
    ### Serialize
    pickle.dump(coverage, open(args.output_pickle, "wb"), protocol=pickle.HIGHEST_PROTOCOL)
