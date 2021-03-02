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

    def binValue(self, i):
        bin = self.bins[i]
        id, min, max, prob = bin
        if self.isDiscrete:
            value = min
        else:
            value = (min + max) / 2
        return value

    def minVal(self):
        return self.min
    
    def maxVal(self):
        return self.max

    def getBinForVal(self, value):
        if value < self.min:
            return 0
        elif value >= self.max:
            return len(self.bins)-1
        for i in range(self.binCount):
            bin = self.bins[i]
            indx, start, end, prob = bin
            if value >= start and value < end:
                return i
        print('no bin found', value, self.min, self.max)

    def P(self, valueSpec):
        """ Return the probability of a given value or range of values.
            valueSpec can be:
            number -- The probability of attaining the given value
                        Note: For real valued continuous variables, the probability
                        of attaining a given value is essentially zero.  In this case,
                        the probability of attaining a value within the discretized bin
                        associated with value is returned.  This is useful for comparison
                        purposes, but the returned probability is dependent on the discretization
                        density.
            (low, high) -- The probability of attaining a value in the range (low <= value < high)
                        high = None means infinity
                        low = None means negative infinity
        """
        outProb = 0.0
        if type(valueSpec) == type((0,)):
            # it's a range tuple
            assert len(valueSpec) == 2, 'pdf.P: valueSpec must be a single number or 2-tuple = ' + str(valueSpec)
            low, high = valueSpec
            if low is None:
                low = self.min - 1
            if low > self.max:
                return 0.0
            if high is None:
                high = self.max + 1
            if high <= self.min:
                return 0.0
            return self.Prange(low, high)
        value = valueSpec
        if value < self.min or value >= self.max:
            return outProb  # Outside the range.  Zero probability.
        for i in range(self.binCount):
            bin = self.bins[i]
            indx, start, end, prob = bin
            if value >= start and value < end:
                outProb = prob
                break
        return outProb

    def Prange(self, minVal, maxVal):
        """ Return the probability of x between 2 values
            i.e. P(minVal <= X < maxVal)
        """
        if minVal < self.min:
            minVal = self.min
        if maxVal > self.max:
            maxVal = self.max
        firstBin = self.getBinForVal(minVal)
        lastBin = self.getBinForVal(maxVal)
        cum = 0.0
        for i in range(firstBin, lastBin + 1):
            bin = self.bins[i]                
            indx, bmin, bmax, prob = bin
            if i == firstBin and i == lastBin:
                adjProb = (maxVal - minVal) / (bmax - bmin) * prob
            elif i == firstBin:
                adjProb = (bmax - minVal) / (bmax - bmin) * prob
            elif i == lastBin:
                adjProb = (maxVal - bmin) / (bmax - bmin) * prob
            else:
                adjProb = prob
            cum += adjProb
        return cum

    def E(self):
        """Return the expected value (e.g. mean) of the disribution."""
        cum = 0
        for i in range(self.binCount):
            bin = self.bins[i]
            id, min, max, prob = bin
            if self.isDiscrete:
                value = min
            else:
                value = (min + max) / 2.0
            value = self.binValue(i)
            cum += prob * value
        return cum

    mean = E
    
    def var(self):
        mean = self.E()
        cum = 0.0
        for i in range(self.binCount):
            bin = self.bins[i]
            id, min, max, prob = bin
            value = self.binValue(i)
            cum += prob * (value - mean)**2
        var = cum
        return var

    def stDev(self):
        var = self.var()
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
            value = self.binValue(i)
            cum += prob * ((value-mean) / std)**3
        return cum

    def kurtosis(self):
        """ Return the excess kurtosis of the distribution"""
        mean = self.E()
        std = self.stDev()
        if std == 0:
            return 0.0
        cum = 0.0
        for i in range(self.binCount):
            bin = self.bins[i]
            id, min, max, prob = bin
            value = self.binValue(i)
            cum += prob * ((value-mean) / std)**4
        
        return cum - 3

    def ToHistogram(self):
        """Convert the pdf to a numpy array of probabilities [P(bin1), ..., P(binN)]"""
        return np.array([bin[3] for bin in self.bins])

    def ToHistTuple(self):
        """Convert the pdf to a list of tuples of binVal, and P(bin) --  [(binVal1, P(bin1)), ..., (binValN, P(binN))]"""
        outTups = []
        for i in range(len(self.bins)):
            bin = self.bins[i]
            binProb = bin[3]
            outTups.append((self.binValue(i), binProb))
        return outTups

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

