# Goes through each experiment, and for a bunch of sets of peaks
# (from different factors), extracts the exact cuts signal, which I will use
# later for summarizing stuff.

project.init2("chipmentation")
eload(loadPSA())
eload(loadPeaks())

# Create caches surrounding ChIP Peaks:
i=24; iFactor=4

# For each experiment of interest
for (i in 1:nrow(SV$msa)) {
	# Grab the experiment exact cuts file
	xcfile = SV$msa[i, xcfile]
	if (!file.exists(xcfile)) {
		warning("DNE: ", xcfile);
		next();
	}
	message(i, " of ", nrow(SV$msa), ": ", SV$msa[i, sampleName]);

	# For each factor (which has a set of peaks)
	for (iFactor in 1:length(SV$factors)) {
		factorName = SV$factors[iFactor]
		cacheSubDir = paste0("signal/", factorName)
		message(factorName, "\t", appendLF=FALSE)	
	
		# Grab the set of peaks
		bedFile = SV[[paste0(factorName, "BedFile")]]
		negStrand = which(SV[[factorName]]$V6=="-")
		pointsToExtract = 2000

		# Extract scores for those peaks
		simpleCache(SV$msa[i, sampleName], { x = bwSummaryOverBed(	SV$msa[i, xcfile], bedFile, nElem=pointsToExtract); xm = as.matrix(x[, -1, with=FALSE]); xm = flipRows(xm, negStrand); xm }, cacheSubDir=cacheSubDir, assignToVar="M", recreate=FALSE, nofail=TRUE)

		# Create capped, binned, and capped5 versions:
		simpleCache(paste0(SV$msa[i, sampleName], "_cap"), { capData(M, .999); }, cacheSubDir=paste0(cacheSubDir, "_cap"), noload=TRUE, recreate=FALSE, nofail=TRUE)
		simpleCache(paste0(SV$msa[i, sampleName], "_bin"), { M[M>1] = 1; M; }, cacheSubDir=paste0(cacheSubDir,"_bin"), noload=TRUE, recreate=FALSE, nofail=TRUE)
		simpleCache(paste0(SV$msa[i, sampleName], "_cap5"), { M[M>5] = 5; M; }, cacheSubDir=paste0(cacheSubDir,"_cap5"), noload=TRUE, recreate=FALSE, nofail=TRUE)

	} # loop through factors
} # loop through experiments

# Load complete caches like this:
# get(psa[i, sampleName])
# simpleCache(psa[i, sampleName], cacheSubDir="tss400bp", reload=TRUE)

