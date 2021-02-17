import numpy as np
from math import sqrt
import copy

epsilon = .001 # The maximum proability value considered to be zero

class PDF:
    """ Construct the PDF from a list of bins, where a bin is defined as:
        (binNumber, min, max, prob).
        binNumber is the zero based number of the bin (0 - numBins-1)
        min and max are the limits of data values within the bin min <= value < max
        prob is the probability of a point falling within the current bin."""
    def __init__(self, numSamples, binList, isDiscrete = False):
        self.N = numSamples
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

    def skew(self):
        mean = self.E()
        std = self.stDev()
        if std == 0:
            return 0.0
        cum = 0.0
        for i in range(self.binCount):
            bin = self.bins[i]
            id, min, max, prob = bin
            if self.isDiscrete:
                value = min
            else:
                value = (min + max) / 2.0
            cum += prob * (value-mean) / std
        return cum / self.N

    def kurtosis(self):
        mean = self.E()
        std = self.stDev()
        if std == 0:
            return 0.0
        cum = 0.0
        for i in range(self.binCount):
            bin = self.bins[i]
            id, min, max, prob = bin
            if self.isDiscrete:
                value = min
            else:
                value = (min + max) / 2
            cum += prob * ((value-mean) / std)**4
        return cum / self.N

    def ToHistogram(self):
        """Convert the pdf to a numpy array of probabilities [P(bin1), ..., P(binN)]"""
        return np.array([bin[3] for bin in self.bins])

    def SetHistogram(self, newHist):
        assert len(newHist) == len(self.bins), "PDF.SetHistogram: Cannot set histogram with different lenght than current distribution.  (new, original) = " + str((len(newHist), len(self.bins)))
        outBins = []
        for i in range(len(self.bins)):
            bin = self.bins[i]
            prob = newHist[i]
            newBin = bin[:-1] + (prob,)
            outBins.append(newBin)
        self.bins = outBins

    def getBin(self, binNum):
        return self.bins[binNum]

    def __add__(self, other):
        sHist = self.ToHistogram()
        oHist = other.ToHistogram()
        outHist = sHist + oHist
        out = copy.copy(self)
        out.SetHistogram(outHist)
        return out

    def __sub__(self, other):
        sHist = self.ToHistogram()
        oHist = other.ToHistogram()
        outHist = sHist - oHist
        out = copy.copy(self)
        out.SetHistogram(outHist)
        return out

    def __mult__(self, other):
        sHist = self.ToHistogram()
        oHist = other.ToHistogram()
        outHist = sHist * oHist
        out = copy.copy(self)
        out.SetHistogram(outHist)
        return out

    def compare(self, other):
        assert len(self.bins) == len(other.bins), "PDF.compare():  Bin sizes must match for each distribution " + str((len(self.bins, len(other.bins))))
        accum = 0.0
        for i in range(len(self.bins)):
            bin1 = self.bins[i]
            bin2 = other.bins[i]
            prob1 = bin1[3]
            prob2 = bin2[3]
            diff = abs(prob1 - prob2)
            accum += diff
        return accum / len(self.bins)

    def isNull(self):
        for val in self.ToHistogram():
            if val > epsilon:
                return False
        return True

    def __eq__(self, other):
        return (self - other).isNull()

