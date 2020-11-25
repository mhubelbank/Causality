import numpy as np

# Calculate P(B | Pa[1], ... Pa[K])
# Given:
# Pa[K] -- A vector of K Boolean values, one for each Pattern
# P(Pa | B) -- A vector of K values
# P(Pa | ~B) -- A vector of K values
# P(B) -- A scalar

def calcProb(Pa, Pa_B, Pa_nB, PB):
    K = len(Pa)
    # Calculate Likelihood = Product(k)(Pa[k] * Pa_B[k] + (1-Pa[k]) * (1-Ba_B[k]))
    aL = np.zeros((K,))
    for k in range(K):
        l = Pa[k] * Pa_B[k] + (1 - Pa[k]) * (1 - Pa_B[k])
        aL[k] = l
    L = np.prod(aL)
    # Now calculate W =  P(Pa) = Product(k)(Pa[k] * )
    aPa_B = np.zeros((K,))
    aPa_nB = np.zeros((K,))
    for k in range(K):
        #print('Pa = ', Pa, ', Pa_B = ', Pa_B)
        Pa_B1 = (Pa[k] * Pa_B[k] + (1 - Pa[k]) * (1 - Pa_B[k]))
        Pa_nB1 = (Pa[k] * Pa_nB[k] + (1 - Pa[k]) * (1 - Pa_nB[k]))
        aPa_B[k] = Pa_B1
        aPa_nB[k] = Pa_nB1
    W = np.prod(aPa_B) * PB + np.prod(aPa_nB) * (1 - PB)
    #print('L = ', L, ', PB = ', PB, ', W = ', W)
    result = L * PB / W
    return result


