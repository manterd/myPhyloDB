import math
import numpy as np
from numpy.random.mtrand import permutation
from database.utils import ordered_set
import multiprocessing as mp
import time
import math
import os


def homova(dm, trtList, perms):
    bigN, ncols = np.shape(dm)
    mySet = ordered_set(trtList)
    trts = len(mySet)

    lt = np.tril(dm, k=-1)
    Bart, allSSWi = Bartlett(mySet, trtList, lt)

    output = mp.Queue()
    if os.name == 'nt':
        numcore = 1
    else:
        numcore = mp.cpu_count()-1 or 1

    processes = [mp.Process(target=pval, args=(perms, numcore, trtList, mySet, lt, Bart, output, x)) for x in range(numcore)]
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
    result += 'Method name:             ' + str('HOMOVA') + '\n'
    result += 'Bartlett:                ' + ("%.3f" % Bart) + '\n'
    result += 'SSwithin/(N-1):                ' + str(allSSWi) + '\n'
    if pvalue < 0.001:
        result += 'P-value:                 ' + '<0.001\n'
    else:
        result += 'P-value:                 ' + ("%.3f" % pvalue) + '\n'
    result += 'Permutations:            ' + str(perms) + '\n'
    result += 'Samples:                 ' + str(bigN) + '\n'
    result += 'Treatments:              ' + str(trts) + '\n'
    return result


def Bartlett(mySet, trtList, dm):
    N = float(len(trtList))
    P = float(len(mySet))

    SSW, num2, den2 = 0.0, 0.0, 0.0
    allSSWi = {}
    for i in mySet:
        sub = 0.0
        groupList = []
        for j in [j for j, x in enumerate(trtList) if x == i]:
            groupList.append(j)

        group = dm[groupList][:, groupList]
        for x in np.nditer(group):
            sub += x**2

        littleN = float(len(groupList))
        SSWi = sub/littleN
        SSW += SSWi
        allSSWi[str(i)] = ("%.3f" % SSWi)

        num2 += (littleN - 1) * math.log(SSWi / (littleN - 1))
        den2 += 1/(littleN - 1) - 1/(N - P)

    num1 = (N - P) * math.log(SSW / (N - P))
    Bart = (num1 - num2) / (1 + (1/(3*(P - 1))) * den2)
    return Bart, allSSWi


def pval(perms, numcore, trtList, mySet, lt, Bart, output, x):
    frac = 0
    seed = int(time.time()) + x
    np.random.mtrand.seed(seed)
    iter = (perms*1.0)/numcore
    for i in range(int(math.floor(iter*x)), int(math.floor(iter*(x+1)))):
        rList = permutation(trtList)
        rBart = Bartlett(mySet, rList, lt)
        if rBart >= Bart:
            frac += 1
    final = {}
    final["frac"] = frac
    final["iter"] = iter
    output.put(final)