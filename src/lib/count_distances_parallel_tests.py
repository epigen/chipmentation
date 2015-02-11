
bam = HTSeq.BAM_Reader(os.path.abspath("/home/arendeiro/data/human/chipmentation/mapped/merged/H3K4me3_K562_500k_CM.bam"))
windows = [HTSeq.GenomicInterval("chr1", 100000,110000), HTSeq.GenomicInterval("chr1", 200000,210000)]

dists = reduce(
    lambda x, y: x + y,
    parmap.map(distances, windows, bam, 1, True, True, False)
)

d1 = distances(windows[0], bam, 1,
    duplicates=True,
    strand_wise=True,
    permutate=False
)


d2 = distances(windows[1], bam, 1,
    duplicates=True,
    strand_wise=True,
    permutate=False
)

assert dists == d1 + d2
