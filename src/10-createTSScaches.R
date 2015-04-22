# Creates caches for signal surrounding TSSs.
project.init2("chipmentation")
#sParam = getSlurmParams()

psa=loadPSA()
tss = loadCageTSS();
negStrand = which(tss$V6=="-")



# Create caches:
for (i in 1:nrow(msa)) {
	xcfile = msa[i, xcfile]
	if (!file.exists(xcfile)) {
		warning("DNE: ", xcfile);
		next();
	}
	message(msa[i, sampleName]);
	msa[i, xcfile]
	simpleCache(msa[i, sampleName], { x = bwSummaryOverBed(	msa[i, xcfile], dat$cage, nElem=400); xm = as.matrix(x[, -1, with=FALSE]); xm = flipRows(xm, negStrand); }, cacheSubDir="tss400bp", noload=TRUE, recreate=FALSE, nofail=TRUE)
}

# Capped and binary caches:
for (i in 1:nrow(msa)) {
	xcfile = msa[i, xcfile]
	if (!file.exists(xcfile)) {
		warning("DNE: ", xcfile);
		next();
	}
	message(msa[i, sampleName]);
	simpleCache(msa[i, sampleName], cacheSubDir="tss400bp", assignToVariable="M", recreate=FALSE, nofail=TRUE, reload=TRUE)
	# set max to 1
	
	simpleCache(paste0(msa[i, sampleName], "_cap"), { capData(M, .999); }, cacheSubDir="tss400bp_capped", noload=TRUE, recreate=FALSE, nofail=TRUE)
	simpleCache(paste0(msa[i, sampleName], "_bin"), { M[M>1] = 1; M; }, cacheSubDir="tss400bp_bin", noload=TRUE, recreate=FALSE, nofail=TRUE)
}



# Test for bad ones:

psa[, xcfile]

for (i in 1:nrow(psa)) {
	message(psa[i, sampleName], "\t" , 
file.exists(paste0(getOption("RCACHE.DIR"), "tss400bp/", psa[i, sampleName],".RData"))
	)
}

# Load complete caches like this:
# 
get(psa[i, sampleName])
simpleCache(psa[i, sampleName], cacheSubDir="tss400bp", reload=TRUE)


