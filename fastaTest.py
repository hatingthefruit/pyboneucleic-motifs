#!/usr/bin/env python3

from collections import Counter
import pnaMotifs as pm
from Bio import SeqIO


fasta = "masked.fasta"

dArr = []
bgCount = Counter()

for record in SeqIO.parse(fasta, "fasta"):
    ba = bytearray(str(record.seq), 'utf-8')
    bgCount.update(ba)
    dArr.append(ba)

total = sum([bgCount[x] for x in bgCount if x != b'N'])

bg = {x: bgCount[x]/total for x in bgCount}
for k in range(4, 13):
    pwm = pm.motifEMOOPS(dArr, k, bg)
    print('--------------')
    print("EM:")
    pm.printMotif(pwm, bgCount.keys(), k)

    pwm = pm.motifGibbsOOPS(dArr, k, bg)
    print("Gibbs:")
    pm.printMotif(pwm, bgCount.keys(), k)
