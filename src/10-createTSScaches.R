# Creates caches for signal surrounding TSSs.

project.init2("chipmentation")

utility("funcGenomeSignals.R")

psa=loadPSA()

# Create caches:
for (i in 1:nrow(psa)) {
	xcfile = psa[i, xcfile]
	if (!file.exists(xcfile)) {
		warning("DNE: ", xcfile);
		next();
	}
	message(psa[i, sampleName]);
	psa[i, xcfile]
	simpleCache(psa[i, sampleName], { x = bwSummaryOverBed(	psa[i, xcfile], dat$cage, nElem=400); xm = as.matrix(x[, -1, with=FALSE]); xm }, cacheSubDir="tss400bp", noload=TRUE, recreate=FALSE, nofail=TRUE)
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


