import sys
import csv

colours = {
    "IGG": "153,153,153", "INPUT": "153,153,153",  # grey
    "H3K36ME1": "158,36,0", "H3K36ME2": "158,36,0", "H3K36ME3": "158,36,0",  # orange
    "H3K4ME3": "0,158,115",  # bluish green
    "H3K4ME1": "230,159,0", "H3K14ac": "230,159,0",  # yellow
    "H3K27ME1": "0,114,178", "H3K27ME2": "0,114,178", "H3K27ME3": "0,114,178",  # blue
    "H3K9ME1": "86,180,233", "H3K9ME2": "86,180,233", "H3K9ME3": "86,180,233",  # sky blue
    "H3AC": "213,94,41", "H3K9AC": "213,94,41", "H3K27AC": "213,94,41", "H3K56AC": "213,94,41", "H3K56AC": "213,94,41",  # vermillion
    "H3K79ME1": "204,121,167", "H3K79ME2": "204,121,167", "H3K79ME3": "204,121,167",  # reddish purple
    "ATAC": "0,158,115", "nan": "0,158,115",
    "DNASE": "0,158,115",
    "CTCF": "83,66,2",
    "PU1": "110,2,44",
    "GATA1": "155,5,5",
    "GATA2": "81,3,3",
    "REST": "37,2,109",
    "CJUN": "42,3,81",
    "FLI1": "81,81,3"
}

wr = sys.stdout
wr.write("browser position chr21:28,049,584-38,023,583\n")

for line in csv.reader(iter(sys.stdin.readline, ''), delimiter='\t'):
    bigwig = line[0].split("/")[-1]
    name = bigwig.split(".")[0]

    ip = name.split("_")[3]
    colour = colours[ip]

    command = """track type=bigWig name='{0}' description='{0}' """.format(name)
    command += """height=32 visibility=full maxHeightPixels=32:32:25 bigDataUrl={0} color={1}\n""".format(line[0], colour)

    wr.write(command)
wr.close()
