#!/bin/python

import pnaMotifs as pm
import random

# Todo: create some simple test data, run tests
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

pwm = pm.motifEMOOPS(dArr, k, {x: .25 for x in alpha})
print('--------------')

for y in alpha:
    print(chr(y), end=" ")
    for x in range(k):
        print("%f" % pwm[x][y], end=' ')
    print()

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

pwm = pm.motifEMOOPS(dArr, k, {x: .25 for x in alpha})
print('--------------')
for y in alpha:
    print(chr(y), end=" ")
    for x in range(k):
        print("%f" % pwm[x][y], end=' ')
    print()
