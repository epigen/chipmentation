ChIPmentation
=====

Schmidl 2015 'ChIPmentation' project
-----
In its first version (prior to submission) most preprocessing was writen in `bash` and downstream analysis was done with `R`. `src/main.sh` was the main script which called some of the code in `src/pipelines,scripts,lib` etc...

The current version uses my Python repository [`chipseq-pipelines`](https://github.com/afrendeiro/chipseq-pipelines) to manage preprocessing of files much more efficiently using only the sample annotation sheet ([`chipmentation_annotation.csv`](chipmentation_annotation.csv)).

All downstram analysis is performed with scripts in `src` using libraries such as HTSeq and pybedtools for efficiency and clarity.


# TODO list (response to reviewers)

### Experimental procedures
- Transposase titration on H3K4me3 ChIPmentation (up to saturation)
- ChIPmentation on TFs with less cells (*e.g.* 100k, 200k, 500k, 1-2mio)
- Tagmentation of H3K4me3 ChIPed DNA on K562
- ATAC-seq on K562
- Paired end sequencing of ChIPmentation (H3K4me3 and PU.1)

### Analysis
- Address transposition event sequence biases at the TATA box with tagmented ChIP DNA (see [relevant document](https://raw.githubusercontent.com/ComputationalEpigenetics/chipmentation/d2172920013bb20dd19101c4cc7795caa0113c18/results/tn5_bias/README.md))
- Address transposition event sequence biases genome-wide (read frequency by k-mers genome-wide?)
- Show ChIPmentation does produce subnucleosomal fragments using paired end data
- Show ChIPmentation is more likely to generate only one (pair) of read from each IPed fragment at least for TFs
- Explain differences of footprints by showing ChIP, ATAC-seq and DNAse data (ChIPmentation should be a mixture of all)
- Compare number of peaks with footprints between DNAase, ATAC-seq and ChIPmentation
- Detect footprints of CM TFs on histone data (H3K4me1, H3K27Ac on TFs)

### Manuscript
- Better tracks to display the data
- Reinforce the idea of novelty by showing the addditional information captured (as reviwer #3 appreciates) in a better way
- Rephrase that ChIPmentation does allow the usage of less cells, but that is dependant on the abundance of the factor and antibody efficiency
- Mention that we're able to recover much smaller fragments due to no adapter dimers (show bioanalyzer profiles)
- Mention we remove the transposase by washing the ChIP with SDS (as previously demonstrated)
- Show library size distribution created with different amounts of transposase
