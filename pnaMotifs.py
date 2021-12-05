from typing import ByteString, Dict, Iterable, List, Union
import random
import math


"""
Expectation Maximization algorithm for motif finding. Given a list of sequences, a motif length, and background character frequencies, it returns a positional weight matrix representating a motif

For this implementation, we choose an OOPS model (Once Occurence Per String). 
"""

def motifEM(sequences: List[bytearray], k: int, bgFreqs: Dict[int, float]) -> List[Dict[int, float]]:
    seqLens = [len(x) for x in sequences]  # The length of each input sequence
    z = [[0 for x in range(seqLens[y]-k)] for y in range(len(sequences))]
    numM = [0 for x in sequences]
    pwmCounts = [{x: 0 for x in bgFreqs} for y in range(k)]

    # Randomly generate initial motifs
    for i, seq in enumerate(sequences):
        # Randomly pick a motif location
        j = random.randint(0, seqLens[i]-k-1)
        numM[i] += 1
        # Update count in each position for PWM
        for x in range(k):
            pwmCounts[x][seq[j+x]] += 1

    # Get the average sequence length and the average number of randomly assigned 'motifs' per sequence
    # to estimate the probability that a substring is a motif
    avgLen = sum(seqLens) / len(sequences)
    count = sum(numM)
    mProb = (count / len(sequences)) / avgLen

    # Generate random PWM values from randomly generated motifs
    pwm = [{x: pwmCounts[i][x]/count for x in pwmCounts[i]} for i in range(k)]
    # To prevent masked characters from interfering with the result, cap the possible probability for 'N' in any position and normalize all other probabilities in that position
    for x in range(k):
        if 78 in pwm[x] and pwm[x][78] > 1e-4:
            pwm[x][78] = 1e-4
            total = sum([pwm[x][y] for y in pwm[x]])
            pwm[x] = {y: pwm[x][y]/total for y in pwm[x]}

    rnd = 0  # Make sure we have at least one round
    logLO = 0  # Old logL
    logL = 0
    # Loop until log-likelihood is at a local optimum, or at least once
    while rnd < 1 or logLO > logL:
        logLO = logL
        logL = 0

        # Reset counts
        pwmCounts = [{x: 0 for x in bgFreqs} for y in range(k)]

        # E-step
        count = 0
        for i, seq in enumerate(sequences):
            for j in range(seqLens[i]-k):
                currZ = mProb
                # Calculate E value for the current location in Z
                for x, char in enumerate(seq[j:j+k]):
                    currZ *= pwm[x][char]
                # Update character counts for the current location in Z. Uxe currZ to take a weighted average later
                for x, char in enumerate(seq[j:j+k]):
                    if char != b'N':
                        pwmCounts[x][char] += currZ
                count += currZ
                # z[i][j] = P(M|substr in i at j) = P(W|M) * lmd
                z[i][j] = currZ
                if currZ != 0:
                    logL += math.log(currZ)

        # M-step
        pwm = [{x: pwmCounts[i][x]/count for x in pwmCounts[i]}
               for i in range(k)]
        # To prevent masked characters from interfering with the result, cap the possible probability for 'N' in any position and normalize all other probabilities in that position
        for x in range(k):
            if 78 in pwm[x] and pwm[x][78] > 1e-4:
                pwm[x][78] = 1e-4
                total = sum([pwm[x][y] for y in pwm[x]])
                pwm[x] = {y: pwm[x][y]/total for y in pwm[x]}
        mProb = (count / len(sequences)) / avgLen

        rnd += 1
    return pwm


"""
Gibbs sampling algorithm for motif finding. Given a list of sequences, a motif length, and background character frequencies, it returns a positional weight matrix representating a motif

For this implementation, we choose an OOPS model (Once Occurence Per String). For updating motif locations, we choose a deterministic (viterbi) method, or simply always picking the most likely location given the current model parameters.
"""


def motifGibbsOOPS(sequences: List[bytearray], k: int, bgFreqs: Dict[int, float]) -> List[Dict[int, float]]:
    seqLens = [len(x) for x in sequences]  # The length of each input sequence
    pwmCounts = [{x: 0 for x in bgFreqs} for y in range(k)]
    z = []

    # Randomly generate initial motif locations. The model for this function is OOPS, so we only generate one per sequence
    for i, seq in enumerate(sequences):
        j = random.randint(0, seqLens[i]-k)
        z.append(j)

    lastZ = []

    # Loop until there is not change in the motif locations
    while z != lastZ:
        # Reset pwmCounts
        pwmCounts = [{x: 0 for x in bgFreqs} for y in range(k)]
        lastZ = list(z)

        # Loop through all sequences
        for i, seq in enumerate(sequences):
            # Update pwmCounts for all sequences except the current on
            for j, int_seq in enumerate(sequences):
                currLoc = z[j]
                if j != i:
                    for x in range(k):
                        pwmCounts[x][int_seq[currLoc+x]] += 1

            # Use updated counts to update actual PWM
            pwm = [{x: pwmCounts[i][x]/(len(sequences) - 1)
                    for x in pwmCounts[i]} for i in range(k)]
            # To prevent masked characters from interfering with the result, cap the possible probability for 'N' in any position and normalize all other probabilities in that position
            for x in range(k):
                if 78 in pwm[x] and pwm[x][78] > 1e-4:
                    pwm[x][78] = 1e-4
                    total = sum([pwm[x][y] for y in pwm[x]])
                    pwm[x] = {y: pwm[x][y]/total for y in pwm[x]}
            # Find location such that Q_x/P_x (ratio of pattern to background probability) is maximized. Update z[i] to this new location, chosen using viterbi method (always the highest)
            max = 0
            maxPos = 0
            for j in range(seqLens[i] - k):
                score = 1
                bgScore = 1
                # Calculate Q_x and P_x concurrently
                for x in range(k):
                    score *= pwm[x][seq[j + x]]
                    bgScore *= bgFreqs[seq[j+x]]

                # Compare to maximimum ratio and update if needed
                ratio = score/bgScore
                if ratio > max:
                    max = ratio
                    maxPos = j
            z[i] = maxPos

    # Once loop has terminated calculate PWM one last time
    # First, reset pwmCounts
    pwmCounts = [{x: 0 for x in bgFreqs} for y in range(k)]

    # Then update counts at each location
    for i, seq in enumerate(sequences):
        currLoc = z[i]
        for x in range(k):
            pwmCounts[x][seq[currLoc+x]] += 1

    # Then calculate probabilities for each location
    pwm = [{x: pwmCounts[i][x]/len(sequences)
            for x in pwmCounts[i]} for i in range(k)]
    # To prevent masked characters from interfering with the result, cap the possible probability for 'N' in any position and normalize all other probabilities in that position
    for x in range(k):
        if 78 in pwm[x] and pwm[x][78] > 1e-4:
            pwm[x][78] = 1e-4
            total = sum([pwm[x][y] for y in pwm[x]])
            pwm[x] = {y: pwm[x][y]/total for y in pwm[x]}
    return pwm


"""
Helper function to make it easier to visualize the resulting PWM from the above algorithms
"""


def printMotif(pwm: List[Dict[int, float]], alpha: Union[Dict[int, float], Iterable[int]], k: int):
    for y in alpha:
        print(chr(y), end=" ")
        for x in range(k):
            print("%f" % pwm[x][y], end=' ')
        print()
