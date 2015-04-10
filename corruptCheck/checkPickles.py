import os
import cPickle as pickle

pickles = os.listdir("/media/afr/cemm-backup/chipmentation/periodicity")

corrupt = []
for pic in pickles:
    print(pic)
    try:
        pickle.load(open(pic, "r"))
    except:
        corrupt.append(pic)
