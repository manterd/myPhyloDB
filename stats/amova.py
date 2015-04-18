import numpy as np
from numpy.random.mtrand import permutation
from database.utils import ordered_set
import multiprocessing as mp
import time
import math
import os


def amova(dm, trtList, perms):
    bigN, ncols = np.shape(dm)
    mySet = ordered_set(trtList)
    trts = float(len(mySet))

    lt = np.tril(dm, k=-1)
    subtotal = 0.0
    for x in np.nditer(lt):
        subtotal += x**2
    SStotal = float(subtotal)/float(bigN)

    SSwithin = SSW(mySet, trtList, lt)
    SSamong = SStotal - SSwithin
    MSamong = SSamong / (trts - 1)
    MSwithin = SSwithin / (bigN - trts)
    Fvalue = MSamong / MSwithin

    output = mp.Queue()
    numcore = mp.cpu_count()-1 or 1

    processes = [mp.Process(target=pval, args=(perms, numcore, trtList, mySet, lt, SStotal, trts, bigN, Fvalue, output, x)) for x in range(numcore)]
    for p in processes:
        p.start()
    for p in processes:
        p.join()

    results = [output.get() for p in processes]
    Ptotal = 0.0
    Ntotal = 0.0
    for j in results:
        Ptotal = np.core.umath.add(Ptotal, j["frac"])
        Ntotal = np.core.umath.add(Ntotal, j["iter"])

    pvalue = Ptotal/Ntotal

    result = ''
    result += 'Method name:             ' + str('AMOVA') + '\n'
    result += 'SSamong:                 ' + ("%.3f" % SSamong) + '\n'
    result += 'SSwithin:                ' + ("%.3f" % SSwithin) + '\n'
    result += 'SStotal:                 ' + ("%.3f" % SStotal) + '\n'
    result += 'MSamong:                 ' + ("%.3f" % MSamong) + '\n'
    result += 'MSwithin:                ' + ("%.3f" % MSwithin) + '\n'
    result += 'F-value:                 ' + ("%.3f" % Fvalue) + '\n'
    if pvalue < 0.001:
        result += 'P-value:                 ' + '<0.001\n'
    else:
        result += 'P-value:                 ' + ("%.3f" % pvalue) + '\n'
    result += 'Permutations:            ' + str(perms) + '\n'
    result += 'Samples:                 ' + str(bigN) + '\n'
    result += 'Treatments:              ' + str(trts) + '\n'
    return result


def SSW(mySet, trtList, dm):
    SSwithin = 0.0
    for i in mySet:
        subtotal = 0.0
        groupList = []
        for j in [j for j, x in enumerate(trtList) if x == i]:
            groupList.append(j)

        group = dm[groupList][:, groupList]
        for x in np.nditer(group):
            subtotal += x**2

        littleN = float(len(groupList))
        SSwithin += subtotal/littleN
    return SSwithin


def pval(perms, cores, trtList, mySet, lt, SStotal, trts, bigN, Fvalue, output, x):
    frac = 0
    seed = int(time.time()) + x
    np.random.mtrand.seed(seed)
    iter = (perms*1.0)/cores
    for i in range(int(math.floor(iter*x)), int(math.floor(iter*(x+1)-1))):
        rList = permutation(trtList)
        rSSwithin = SSW(mySet, rList, lt)
        rSSamong = SStotal - rSSwithin
        rMSamong = rSSamong / (trts - 1)
        rMSwithin = rSSwithin / (bigN - trts)
        rFvalue = rMSamong / rMSwithin
        if rFvalue >= Fvalue:
            frac += 1
    final = {}
    final["frac"] = frac
    final["iter"] = iter
    output.put(final)