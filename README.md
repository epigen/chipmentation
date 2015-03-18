ChIPmentation
=====

Schmidl 2015 'ChIPmentation' project
-----
In its first version (prior to submission) most preprocessing was writen in `bash` and downstream analysis was done with `R`. `src/main.sh` was the main script which called some of the code in `src/pipelines,scripts,lib` etc...

The current version uses my Python repository [`chipseq-pipelines`](https://github.com/afrendeiro/chipseq-pipelines) to manage preprocessing of files much more efficiently using only the sample annotation sheet ([`chipmentation.sample_annotation.csv`](chipmentation.sample_annotation.csv)).

All downstream analysis is performed with scripts in `src` using libraries such as HTSeq and pybedtools for efficiency and clarity.

### Data
See annotation sheet here: [CSV](chipmentation.sample_annotation.csv)

See tracks here: [UCSC](http://genome.ucsc.edu/cgi-bin/hgTracks?org=human&hgt.customText=http://www.biomedical-sequencing.at/bocklab/arendeiro/chipmentation/bigWig/trackHub_hg19.curated.txt)

### Sample names
Sample names are a concatenation of the annotation fields: `cellLine, numberCells, technique, ip, patient, treatment, biologicalReplicate, technicalReplicate, genome`. 

Samples where all of the above are the same except for `technicalReplicate` (same sample sequenced twice independently), are merged in a new sample labeled with `technicalReplicate=0`, representing one biological replicate which was sequenced more than once. 

Further, Samples where all of the above are the same except for `technicalReplicate` AND `biologicalReplicate` , are merged in a new sample labeled with `technicalReplicate=0`, `biologicalReplicate=0`, representing a concatenation of all biological replicates of this type.

For downstream analysis I've been using mostly the later case: `technicalReplicate=0`, `biologicalReplicate=0`.

This creates a [***sheet with new samples*** containing these merged ones](https://github.com/ComputationalEpigenetics/chipmentation/blob/master/chipmentation.replicates.annotation_sheet.csv), and paired control samples, which I use for downstream.

### Project structure
As defined in [`chipseq-pipelines`](https://github.com/afrendeiro/chipseq-pipelines):

For all my projects:

`projectsroot`=/fhgfs/groups/lab_bock/shared/projects/
`htmlroot`=/fhgfs/groups/lab_bock/public_html/arendeiro/

```
projectsroot
|__ chipmentation
    |__ runs
    |__ data
    |   |__ fastq
    |   |__ fastqc
    |   |__ raw
    |   |__ mapped
    |   |__ coverage
    |   |__ peaks
    |   |__ motifs
    |__ results
         |__ plots

htmlroot
|__ chipmentation
    |__ bigwig
```
JSON description [here](https://github.com/ComputationalEpigenetics/chipmentation/blob/master/projectPaths.json).

Probably the most relevant folder is `/chipmentation/data/mapped`. There are 2 bam files for each sample (+ .bai indices):
- ChIP-seq samples:
    - `.trimmed.bowtie2.dups.bam` - samples trimmed, aligned and with marked duplicates
    - `.trimmed.bowtie2.nodups.bam` - samples trimmed, aligned and with duplicates removed
- ChIPmenmtation, ATAC-seq or ChIPed-tagmented naked DNA samples:
    - `.trimmed.bowtie2.shifted.dups.bam` - samples trimmed, aligned, reads shifted and with marked duplicates
    - `.trimmed.bowtie2.shifted.nodups.bam` - samples trimmed, aligned, reads shifted and with duplicates removed

I will document the other folders in due time.

# TODO list (response to reviewers)
### Experimental procedures
- Transposase titration on H3K4me3 ChIPmentation (up to saturation)
- ChIPmentation on TFs with less cells (*e.g.* 100k, 200k, 500k, 1-2mio)
- Tagmentation of H3K4me3 ChIPed DNA on K562
- ATAC-seq on K562
- Paired end sequencing of ChIPmentation (H3K4me3 and PU.1)

### Analysis
- Address transposition event sequence biases at the TATA box with tagmented ChIP DNA (see [relevant document](https://github.com/ComputationalEpigenetics/chipmentation/blob/master/results/tn5_bias/README.md))
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
- Show library size distribution created with different amounts of transposase (see [relevant document](https://github.com/ComputationalEpigenetics/chipmentation/tree/master/results/fragment_size))
