#!/usr/bin/env python

from typing import Counter
import pnaMotifs as pm
import random

numData = random.randint(10, 40)
minLen = 15
maxLen = 80
dArr = []

alpha = b'ACGT'
k = 8
mStr = [random.choice(alpha) for x in range(k)]
print("Actual motif:", [chr(x) for x in mStr])
actM = bytearray(mStr)

for x in range(numData):
    currSize = random.randint(minLen, maxLen)
    currStr = bytearray(currSize)
    first = random.randint(0, currSize - k)
    currStr[first:first+k] = actM
    for y in range(first):
        currStr[y] = random.choice(alpha)
    for y in range(first+k, currSize):
        currStr[y] = random.choice(alpha)
    dArr.append(currStr)

pwm1 = pm.motifEM(dArr, k, {x: .25 for x in alpha})
pwm2 = pm.motifGibbsOOPS(dArr, k, {x: .25 for x in alpha})
print('--------------')

print("EM:")
pm.printMotif(pwm1, alpha, k)
print("Gibbs:")
pm.printMotif(pwm2, alpha, k)

dArr=[]

# With one replacement
for x in range(numData):
    currSize = random.randint(minLen, maxLen)
    currStr = bytearray(currSize)
    first = random.randint(0, currSize - k)
    currStr[first:first+k] = actM
    for i in range(4):
        currStr[first+random.randint(0, k-1)] = random.choice(alpha)
    for y in range(first):
        currStr[y] = random.choice(alpha)
    for y in range(first+k, currSize):
        currStr[y] = random.choice(alpha)
    dArr.append(currStr)

pwm = pm.motifEM(dArr, k, {x: .25 for x in alpha})
print('--------------')
print("EM:")
pm.printMotif(pwm, alpha, k)

pwm = pm.motifGibbsOOPS(dArr, k, {x: .25 for x in alpha})
print("Gibbs:")
pm.printMotif(pwm, alpha, k)