#!/usr/bin/env python3

from collections import Counter
import pnaMotifs as pm
from pysam import FastaFile


fasta = "masked.fasta"
seq_obj = FastaFile(fasta)

dArr = []
alpha = b'ACGT'
bgCount = Counter(alpha)
for ref in seq_obj.references:
    ba = bytearray(seq_obj.fetch(ref), 'ascii')
    bgCount.update(ba)
    dArr.append(ba)

total = sum([bgCount[x] for x in bgCount if x != b'N'])

bg = {x: bgCount[x]/total for x in bgCount}
bg.update({78:.001})
for k in range(4, 13):
    pwm = pm.motifEMOOPS(dArr, k, bg)
    print('--------------')
    print("EM:")
    pm.printMotif(pwm, alpha, k)

    pwm = pm.motifGibbsOOPS(dArr, k, bg)
    print("Gibbs:")
    pm.printMotif(pwm, alpha, k)
