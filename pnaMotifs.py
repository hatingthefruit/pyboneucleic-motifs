from typing import ByteString, Dict, List, Match
import random
import math


def motifEMOOPS(sequences: List[ByteString], k: int, bgFreqs: Dict[int, float]):
    seqLens = [len(x) for x in sequences]  # The length of each input sequence
    z = [[0 for x in range(seqLens[y]-k)] for y in range(len(sequences))]
    numM = [0 for x in sequences]
    pwmCounts = [{x: 0 for x in bgFreqs} for y in range(k)]

    # Randomly generate initial motifs
    for i, seq in enumerate(sequences):
        for j in range(seqLens[i]-k):
            # Randomly choose if current subsequence is a motif
            if bool(random.getrandbits(1)):
                numM[i] += 1
                # Update count in each position for PWM
                for x in range(k):
                    pwmCounts[x][seq[j+x]] += 1
                break

    # Get the average sequence length and the average number of randomly assigned 'motifs' per sequence
    # to estimate the probability that a substring is a motif
    avgLen = sum(seqLens) / len(sequences)
    count = sum(numM)
    mProb = (count / len(sequences)) / avgLen

    # Generate random PWM values from randomly generated motifs
    pwm = [{x: pwmCounts[i][x]/count for x in pwmCounts[i]} for i in range(k)]

    
    rnd = 0
    logLO = 0
    logL = 0
    while rnd < 2 or logLO > logL:
        logLO = logL
        logL = 0
        pwmCounts = [{x: 0 for x in bgFreqs} for y in range(k)]
        # E-step
        count = 0
        for i, seq in enumerate(sequences):
            for j in range(seqLens[i]-k):
                currZ = mProb
                for x, char in enumerate(seq[j:j+k]):
                    currZ *= pwm[x][char]
                for x, char in enumerate(seq[j:j+k]):
                    pwmCounts[x][char] += currZ
                count += currZ
                # z[i][j] = P(M|substr in i at j) = P(W|M) * lmd
                z[i][j] = currZ
                if currZ != 0:
                    logL += math.log(currZ)

        # M-step
        pwm = [{x: pwmCounts[i][x]/count for x in pwmCounts[i]}
               for i in range(k)]
        mProb = (count / len(sequences)) / avgLen

        rnd += 1

    return pwm
