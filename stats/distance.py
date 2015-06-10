import numpy as np
from numpy import *

#TODO fix crash if division by 0
def wOdum(data, alpha):
    length = data.__len__()
    numrows, numcols = shape(data)
    dists = np.zeros((numrows, numrows))

    for i in xrange(length):
        for j in xrange(length):
            dist = 0.0
            if i == j:
                dists[i,j] = dist
                dists[j,i] = dist
            else:
                num, den = 0.0, 0.0
                otus = data[i].__len__()
                for l in xrange(otus):
                    u = data[i,l]
                    v = data[j,l]
                    if u > 0 or v > 0:
                        num += abs(u-v)/(u+v)*(u+v)**alpha
                        den += (u+v)**alpha
                dist = num / den
                dists[i,j] = dist
                dists[j,i] = dist
    return dists


def MorisitaHorn(data):
    length = data.__len__()
    numrows, numcols = shape(data)
    dists = np.zeros((numrows, numrows))

    for i in xrange(length):
        n = data[i].sum()
        for j in xrange(length):
            m = data[j].sum()
            dist = 0.0
            if i == j:
                dists[i,j] = dist
                dists[j,i] = dist
            else:
                otus = data[i].__len__()
                num, den1, den2 = 0.0, 0.0, 0.0
                for l in xrange(otus):
                    u = data[i,l] / float(n)
                    v = data[j,l] / float(m)
                    num += u * v
                    den1 += u
                    den2 += v
                dist = 1 - 2 * (num / (den1 + den2))
                dists[i,j] = dist
                dists[j,i] = dist
    return dists
