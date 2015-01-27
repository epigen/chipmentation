ChIPmentation
=====

Schmidl 2015 'ChIPmentation' project
-----
In its first version (prior to submission) most preprocessing was writen in `bash` and downstream analysis was done with `R`. `src/main.sh` was the main script which called some of the code in `src/pipelines,scripts,lib` etc...

The current version uses my Python repository [`chipseq-pipelines`](https://github.com/afrendeiro/chipseq-pipelines) to manage preprocessing of files much more efficiently using only the sample annotation sheet ([`chipmentation_annotation.csv`](chipmentation_annotation.csv)).

All downstram analysis is performed with scripts in `src` using libraries such as HTSeq and pybedtools for efficiency and clarity.
