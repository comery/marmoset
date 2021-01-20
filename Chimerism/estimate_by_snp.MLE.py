#!/usr/bin/env python3

import sys
import numpy as np
from scipy.stats import binom
from scipy.optimize import minimize

if len(sys.argv) != 2:
    sys.exit('python3 %s <genotype_genomicDBI_gather.g.vcf.tab.snp.mat.filt>' % (sys.argv[0]))

inFile = sys.argv[1]

genotype_prop = {'00': 36/256, '01': 24/256, '02': 4/256, '10': 24/256, '11': 80/256, '12': 24/256, '20': 4/256, '21': 24/256, '22': 36/256}

def prop(rho, m1, m2):
    return 1/2*(rho*m1 + (1-rho)*m2)

def lik(rho):
    L = np.zeros(len(d1_array))
    '''
    for i in [0,1,2]:
        for j in [0,1,2]:
            p = prop(rho, i, j)
            C = 1 - binom.pmf(0, n_array, p) - binom.pmf(n_array, n_array, p)
            P = binom.pmf(d1_array, n_array, p)/C
            L += P*genotype_prop[str(i)+str(j)]
    '''
    #36/256 * binom.pmf(d1_array, n_array, prop(rho, 0, 0)) / (1 - binom.pmf(0, n_array, prop(rho, 0, 0)) - binom.pmf(n_array, n_array, prop(rho, 0, 0))) + \
    L = 24/256 * binom.pmf(d1_array, n_array, prop(rho, 0, 1)) / (1 - binom.pmf(0, n_array, prop(rho, 0, 1)) - binom.pmf(n_array, n_array, prop(rho, 0, 1))) + \
        4/256  * binom.pmf(d1_array, n_array, prop(rho, 0, 2)) / (1 - binom.pmf(0, n_array, prop(rho, 0, 2)) - binom.pmf(n_array, n_array, prop(rho, 0, 2))) + \
        24/256 * binom.pmf(d1_array, n_array, prop(rho, 1, 0)) / (1 - binom.pmf(0, n_array, prop(rho, 1, 0)) - binom.pmf(n_array, n_array, prop(rho, 1, 0))) + \
        80/256 * binom.pmf(d1_array, n_array, prop(rho, 1, 1)) / (1 - binom.pmf(0, n_array, prop(rho, 1, 1)) - binom.pmf(n_array, n_array, prop(rho, 1, 1))) + \
        24/256 * binom.pmf(d1_array, n_array, prop(rho, 1, 2)) / (1 - binom.pmf(0, n_array, prop(rho, 1, 2)) - binom.pmf(n_array, n_array, prop(rho, 1, 2))) + \
        4/256  * binom.pmf(d1_array, n_array, prop(rho, 2, 0)) / (1 - binom.pmf(0, n_array, prop(rho, 2, 0)) - binom.pmf(n_array, n_array, prop(rho, 2, 0))) + \
        24/256 * binom.pmf(d1_array, n_array, prop(rho, 2, 1)) / (1 - binom.pmf(0, n_array, prop(rho, 2, 1)) - binom.pmf(n_array, n_array, prop(rho, 2, 1)))
    #36/256 * binom.pmf(d1_array, n_array, prop(rho, 2, 2)) / (1 - binom.pmf(0, n_array, prop(rho, 2, 2)) - binom.pmf(n_array, n_array, prop(rho, 2, 2)))

    logLL = np.log(L).sum()
    print('##rho: %f, p: %f, logLL: %f' % (rho, prop(rho, 2, 0), logLL))
    return -logLL

d1_list = []
d2_list = []

with open(inFile) as f:
    for line in f:
        line = line.rstrip()
        if "CHROM" in line:
            continue
        tmp = line.split('\t')
        ref = tmp[3]
        alt = tmp[4].split(',')
        deps = tmp[7].split(',')
        deps_new = []
        for d in deps:
            if int(d) > 0:
                deps_new.append(int(d))
        if len(deps_new) >= 2:
            d1 = deps_new[0]
            d2 = deps_new[1]
            d1_list.append(d1)
            d2_list.append(d2)

d1_array = np.array(d1_list)
d2_array = np.array(d2_list)
n_array = d1_array + d2_array

lik_model = minimize(lik, np.array([0.2]), method='L-BFGS-B', bounds = ((0.0+1e-6, 0.5-1e-6),))
print(lik_model)

