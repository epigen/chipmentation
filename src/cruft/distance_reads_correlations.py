from matplotlib import pyplot

peaks = [
	[ # peak
		1,2,3,4 # windows
	],
	[
		2,3,4,5
	]	
]


correlations = []

for peak in peaks:
	for window in range(0, len(peak)):
		distReads = {}
		for bp in range(0, len(peak)):
			distance = abs(bp - window)
			reads = int(peak[bp])
			if distance not in distReads.keys():
				distReads[distance] = reads
			else:
				distReads[distance] = distReads[distance] + reads
		# calculate correlation (pseudocode)
		correlation = cor(distReads.keys(), distReads.values())
		correlations.append((window, correlation))

# plot correlations vs distances
pyplot.plot(correlations)

