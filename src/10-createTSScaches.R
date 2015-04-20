project.init2("chipmentation")

# Here's the new fast way to do this. Might want to cache this.
utility("funcGenomeSignals.R")

for (i in 1:nrow(psa)) {
	xcfile = psa[i, xcfile]
	if (!file.exists(xcfile)) {
		warning("DNE: ", xcfile);
		next();
	}
	message(psa[i, sampleName]);
	simpleCache(psa[i, sampleName], { x = bwSummaryOverBed("/media/nsheffield/red6/chipmentation/data/shiftedExactCuts/K562_100K_CM_CTCF_nan_nan_0_0_hg19.bigWig", "~/github/chipmentation/cagePeaks400bp.bed", nElem=400); xm = as.matrix(x[, -1, with=FALSE]); xm }, cacheSubDir="tss400bp", noload=TRUE, recreate=TRUE)
}



# And loading...
# 
get(psa[i, sampleName])
simpleCache(psa[i, sampleName], cacheSubDir="tss400bp", reload=TRUE)
