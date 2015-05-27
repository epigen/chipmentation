ChIPmentation
=====

Schmidl 2015 'ChIPmentation' project
-----
In its first version (prior to submission) most preprocessing was writen in `bash` and downstream analysis was done with `R`. `src/main.sh` was the main script which called some of the code in `src/pipelines,scripts,lib` etc...

The current version uses my Python repository [`chipseq-pipelines`](https://github.com/afrendeiro/chipseq-pipelines) to manage preprocessing of files much more efficiently using only the sample annotation sheet ([`chipmentation.sample_annotation.csv`](metadata/chipmentation.sample_annotation.csv)).

All downstream analysis is performed with scripts in `src` using libraries such as HTSeq and pybedtools for efficiency and clarity.

### Data
See annotation sheet here: [CSV](metadata/chipmentation.replicates.annotation_sheet.csv)

[ChIpmentation TrackHub at UCSC](http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&hubClear=http://biomedical-sequencing.at/bocklab/papers/schmidl2015/hub.txt&position=chr16:88485000-88885000)

### Sample names
Sample names are a concatenation of the annotation fields: `cellLine, numberCells, technique, ip, patient, treatment, biologicalReplicate, technicalReplicate, genome`. 

Samples where all of the above are the same except for `technicalReplicate` (same sample sequenced twice independently), are merged in a new sample labeled with `technicalReplicate=0`, representing one biological replicate which was sequenced more than once. 

Further, Samples where all of the above are the same except for `technicalReplicate` AND `biologicalReplicate` , are merged in a new sample labeled with `technicalReplicate=0`, `biologicalReplicate=0`, representing a concatenation of all biological replicates of this type.

For downstream analysis I've been using mostly the later case: `technicalReplicate=0`, `biologicalReplicate=0`.

This creates a [***sheet with new samples*** containing these merged ones](metadata/chipmentation.replicates.annotation_sheet.csv), and paired control samples, which I use for downstream.

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

JSON description [here](metadata/projectPaths.json).

Probably the most relevant folder is `/chipmentation/data/mapped`. There are 2 bam files for each sample (+ .bai indices):
- ChIP-seq samples:
    - `.trimmed.bowtie2.dups.bam` - samples trimmed, aligned and with marked duplicates
    - `.trimmed.bowtie2.nodups.bam` - samples trimmed, aligned and with duplicates removed
- ChIPmenmtation, ATAC-seq or ChIPed-tagmented naked DNA samples:
    - `.trimmed.bowtie2.shifted.dups.bam` - samples trimmed, aligned, reads shifted and with marked duplicates
    - `.trimmed.bowtie2.shifted.nodups.bam` - samples trimmed, aligned, reads shifted and with duplicates removed

I will document the other folders in due time.
