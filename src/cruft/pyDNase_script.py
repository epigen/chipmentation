
# prepare bed track
awk -v OFS='\t' '{print $1, $2, $3, $4, $5, "+"}' /home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.bed > tmp
cp tmp /home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.bed

ipython

import matplotlib.pyplot as plt
import pyDNase
from pyDNase.footprinting import wellington


bam = "/home/arendeiro/data/human/chipmentation/mapped/merged/PU1_K562_10mio_CM.bam"
bed = "/home/arendeiro/data/human/chipmentation/peaks/PU1_K562_10mio_CM_peaks/PU1_K562_10mio_CM_peaks.bed"


#Load test data
reads = pyDNase.BAMHandler(bam) # pyDNase.example_reads()
regions = pyDNase.GenomicIntervalSet(bed) # pyDNase.example_regions()

#Plot cuts data
for i in range(1, 21):
	plt.subplot(10, 2, i)
	r = random.randint(0, len(regions))
	plt.plot(reads[regions[r]]["+"],c="red")
	plt.plot(-reads[regions[r]]["-"],c="blue")

#Footprint and plot the results
plt.subplot(2, 2, 1)
plt.plot(reads[regions[0]]["+"],c="red")
plt.plot(-reads[regions[0]]["-"],c="blue")
plt.ylabel('reads')
plt.subplot(2, 2, 3)
footprinter = wellington(regions[0], reads)
plt.plot(footprinter.scores,c="black")
plt.ylabel('footprint p-value')

plt.subplot(2, 2, 2)
plt.plot(reads[regions[1]]["+"],c="red")
plt.plot(-reads[regions[1]]["-"],c="blue")
plt.ylabel('reads')
plt.subplot(2, 2, 4)
footprinter = wellington(regions[1], reads)
plt.plot(footprinter.scores,c="black")
plt.ylabel('footprint p-value')

plt.show()



# Footprint all keep ones with cuttoff below -150
footprints = []
for region in range(len(regions)):
	footprinter = wellington(regions[region], reads)
	footprints.append(footprinter.footprints(withCutoff = -20))

with open("PU1_K562_10mio_CM.footprints.bed","w") as bedout:
	for foot in footprints:
		bedout.write(str(foot))

# bigwig - incomplete
print "fixedStep\tchrom=" + str(footprinter.interval.chromosome) + "\t start="+ str(footprinter.interval.startbp) +"\tstep=1"
for i in footprinter.scores:
	print i


