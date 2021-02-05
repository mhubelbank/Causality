import numpy as np
from math import sqrt

class PDF:
    """ Construct the PDF from a list of bins, where a bin is defined as:
        (binNumber, min, max, prob).
        binNumber is the zero based number of the bin (0 - numBins-1)
        min and max are the limits of data values within the bin min <= value < max
        prob is the probability of a point falling within the current bin."""
    def __init__(self, binList, isDiscrete = False):
        self.bins = binList
        self.isDiscrete = isDiscrete
        self.binCount = len(self.bins)
        self.min = binList[0][1] # min of bin 0
        self.max = binList[self.binCount - 1][2] # max of last bin


    def P(self, value):
        """Return the probability of a given value."""
        outProb = 0.0 
        if value < self.min or value >= self.max:
            return outProb  # Outside the range.  Zero probability.
        for i in range(self.binCount):
            bin = self.bins[i]
            indx, start, end, prob = bin
            if value >= start and value < end:
                outProb = prob
                break
        return outProb
        
    def E(self):
        """Return the expected value (e.g. mean) of the disribution."""
        cum = 0
        for i in range(self.binCount):
            bin = self.bins[i]
            id, min, max, prob = bin
            if self.isDiscrete:
                value = min
            else:
                value = (min + max) / 2
            cum += prob * value

        exp = cum
        return exp

    def stDev(self):
        mean = self.E()
        cum = 0.0
        for i in range(self.binCount):
            bin = self.bins[i]
            id, min, max, prob = bin
            if self.isDiscrete:
                value = min
            else:
                value = (min + max) / 2
            
            cum += prob * (value - mean)**2
        var = cum
        std = sqrt(var)
        return std

    def ToHistogram(self):
        """Convert the pdf to a numpy array of probabilities [P(bin1), ..., P(binN)]"""
        return np.array([bin[3] for bin in self.bins])

    def getBin(self, binNum):
        return self.bins[binNum]

        