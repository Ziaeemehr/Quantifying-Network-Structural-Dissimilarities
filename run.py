import numpy as np
import rpy2.robjects as ro
from time import time


def Distances(g, h):
    r = ro.r
    r.source(r_file)
    p = r.D(g, h, 0.45, 0.45, 0.1)
    return list(p)[0]

start = time()
r_file = "d_distance.R"
D = Distances("data/g1.graphml", "data/g2.graphml")
print(D)

print("Done in {:g} seconds".format(time()-start))

start = time()
D = Distances("data/g1.graphml", "data/g2.graphml")
print("Done in {:g} seconds".format(time()-start))