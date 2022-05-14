
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('fname', metavar='fname', type=str, help='input filename')
parser.add_argument('size', metavar='N', type=int, help='picture width&height')
parser.add_argument('outfname', metavar='outfname', type=str,
    help='output filename (no extension)', default="out",  nargs='?')


args = parser.parse_args()
n = args.size
fname = args.fname
outfname = args.outfname
with open(fname, "r") as f:
    mylist = f.readlines()
    data = np.array(mylist, dtype=np.double).reshape((n,n))
    
    plt.imshow(data)
    plt.imsave(outfname + ".png",data) 
    print("done")