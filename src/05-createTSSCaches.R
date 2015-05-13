# Creates caches for signal surrounding TSSs.
project.init2("chipmentation")
#sParam = getSlurmParams()

eload(loadPSA())
eload(loadCageTSS())

# Create caches surrounding CAGE TSSs:

# For each experiment of interest
for (i in 1:nrow(SV$msa)) {
	# Grab the experiment exact cuts file
	xcfile = SV$msa[i, xcfile]
	if (!file.exists(xcfile)) {
		warning("DNE: ", xcfile);
		next();
	}
	message(SV$msa[i, sampleName]);
	SV$msa[i, xcfile]
	simpleCache(SV$msa[i, sampleName], { 
		x = bwSummaryOverBed(SV$msa[i, xcfile], dat$tss, nElem=400);
		xm = as.matrix(x[, -1, with=FALSE]);
		xm = flipRows(xm, negStrand);
		xm 
	}, cacheSubDir="signal/tss", noload=TRUE, recreate=FALSE, nofail=TRUE)
}

# Capped and binary caches:
for (i in 1:nrow(SV$msa)) {
	xcfile = SV$msa[i, xcfile]
	if (!file.exists(xcfile)) {
		warning("DNE: ", xcfile);
		next();
	}
	message(SV$msa[i, sampleName]);
	simpleCache(SV$msa[i, sampleName], cacheSubDir="signal/tss", assignToVariable="M", recreate=FALSE, nofail=TRUE, reload=TRUE)
	# set max to 1
	
	simpleCache(paste0(SV$msa[i, sampleName], "_cap"), { capData(M, .999); }, cacheSubDir="signal/tss_cap", noload=TRUE, recreate=FALSE, nofail=TRUE)
	simpleCache(paste0(SV$msa[i, sampleName], "_bin"), { M[M>1] = 1; M; }, cacheSubDir="signal/tss_bin", noload=TRUE, recreate=FALSE, nofail=TRUE)
	simpleCache(paste0(SV$msa[i, sampleName], "_cap5"), { M[M>5] = 5; M; }, cacheSubDir="signal/tss_cap5", noload=TRUE, recreate=FALSE, nofail=TRUE)
}



# Test for bad ones:

psa[, xcfile]

for (i in 1:nrow(psa)) {
	message(psa[i, sampleName], "\t" , 
file.exists(paste0(getOption("RCACHE.DIR"), "signal/tss/", psa[i, sampleName],".RData"))
	)
}


# I didn't use SLURM as it completed fast enough in just a loop, but here's
# some code in case I need to revisit this in the future:


preamble = substitute ({
	project.init("chipmentation", "projects/chipmentation")
	eload(loadPSA())
	eload(loadPeaks())
	factors = c("ctcf", "pu1", "gata1", "rest")
	iFactor=1
	i=48
} )
sParam = getSlurmParams(preamble=preamble, sourceProjectInit=FALSE)
sParam
dir.create("slurm")

simpleCache(SV$msa[i, sampleName], { x = bwSummaryOverBed(	SV$msa[i, xcfile], bedFile, nElem=pointsToExtract); xm = as.matrix(x[, -1, with=FALSE]); xm = flipRows(xm, negStrand); xm }, cacheSubDir=paste0("signal/", factorName), noload=TRUE, recreate=FALSE, nofail=TRUE, slurmParams=sParam)







