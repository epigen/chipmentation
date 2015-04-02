import csv
import pandas as pd

df = pd.read_csv("/home/arendeiro/lab_bock/public_html/arendeiro/chipmentation/bigWig/trackHub_hg19.curated.txt", sep="\t")

names = df.apply(lambda x: x[0].split('\'')[1], axis=1)

attrs = list(names.apply(lambda x: x.split('_')))

order = pd.DataFrame(sorted(enumerate(attrs), key=lambda x: (x[1][3], x[1][1], x[1][2], x[1][0])))[0]

df = df.reindex(order)

df.to_csv("/home/arendeiro/lab_bock/public_html/arendeiro/chipmentation/bigWig/trackHub_hg19.curated.txt", index=False, quote=csv.QUOTE_NONE)
