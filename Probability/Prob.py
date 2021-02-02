# Module to compute field aggregates for a dataset
import numpy as np
import math
from Probability.pdf import PDF

class Sample:
    def __init__(self, ds, density = 1.0, discSpecs = None):
        """ Create a probability sample object.  Sample provides a mechanism for analyzing
            the probability space of a multivariate dataset.
            It can handle discrete as well as continuous variables.  Continuous probabilities
            are managed by discretization (binning) continuous variables into ranges.
            By default, the number of discretized bins for continuous variables is the square
            root of the number of samples.  This can be increased or decreased using the "density"
            parameter.
            Data(ds) is provided as a dictionary of variable name -> ValueList (List of variable values).
            ValueList should be the same length for all variables.  Each index in ValueList represents
            one sample.
            The density parameter is used to increase or decrease the number of bins used for continuous
            variables.  If density is 1 (default), then sqrt(N-samples) bins are used.  If density is
            set to D, then D * sqrt(N-samples) bins are used. 
            The discSpecs parameter is used to make recursive calls to the module while
            maintaining the discretization information, and should not be provided by the user.
        """
        self.ds = ds
        self.density = density
        self.fieldList = list(ds.keys())
        self._discreteVars = self._getDiscreteVars()
        self.fieldIndex = {}
        for i in range(len(self.fieldList)):
            key = self.fieldList[i]
            self.fieldIndex[key] = i
        # Convert to Numpy array
        npDat = []
        for field in self.fieldList:
            npDat.append(ds[field])
        self.aData = np.array(npDat)
        self.N = self.aData.shape[1]
        #print('N = ', self.N)
        self.fieldAggs = self.getAgg(ds)
        if discSpecs:
            self.discSpecs = self.fixupDiscSpecs(discSpecs)
            #print('self.aData = ', self.aData)
            self.discretized = self.aData
        else:
            self.discSpecs = self.calcDiscSpecs()
            self.discretized = self.discretize()

    def getAgg(self, ds):
        fieldList = list(ds.keys())

        #print('shape = ', aData.shape)
        mins = self.aData.min(1)
        maxs = self.aData.max(1)
        means = self.aData.mean(1)
        #print('means = ', means)
        stds = self.aData.std(1)
        outDict = {}
        for i in range(self.aData.shape[0]):
            fieldName = fieldList[i]
            aggs = (mins[i], maxs[i], means[i], stds[i])
            outDict[fieldName] = aggs
        return outDict

    def calcDiscSpecs(self):
        discSpecs = []
        for i in range(len(self.fieldList)):
            var = self.fieldList[i]
            minV = np.min(self.aData[i])
            maxV = np.max(self.aData[i])
            isDiscrete = var in self._discreteVars
            if isDiscrete:
                minV = int(minV)
                maxV = int(maxV)
                nBins = maxV - minV + 1
                binStarts = [v for v in range(minV, maxV + 1)]
                vals, counts = np.unique(self.aData[i], return_counts = True)
                hist = np.zeros((len(binStarts),))
                for j in range(len(vals)):
                    val = int(vals[j])
                    count = counts[j]
                    hist[val] = count
                edges = [v for v in range(minV, maxV + 2)]
            else:
                if self.N < 100:
                    nBins = 10
                else:
                    nBins = int(self.density * math.sqrt(self.N))
                hist, edges = np.histogram(self.aData[i], nBins)
            discSpecs.append((nBins, minV, maxV, edges, hist, isDiscrete))
        return discSpecs

    def fixupDiscSpecs(self, discSpecs):
        outSpecs = []
        for i in range(len(discSpecs)):
            discSpec = discSpecs[i]
            bins, min, max, edges, hist, isDiscrete = discSpec
            # Regenerate histogram.  The other data should use the original
            newHist, newEdges = np.histogram(self.aData[i], bins, (min, max))
            #assert edges == newEdges, 'Bad Edges'
            outSpecs.append((bins, min, max, edges, newHist, isDiscrete))
        #print('fixedDiscSpecs = ', outSpecs)
        return outSpecs

    def discretize(self):
        discretized = np.copy(self.aData)
        for i in range(len(self.fieldList)):
            field = self.aData[i]
            edges = self.discSpecs[i][3]
            dField = np.digitize(field, edges[:-1]) - 1
            discretized[i] = dField
        return discretized

    def toOriginalForm(self, discretized):
        data = {}
        for f in range(len(self.fieldList)):
            fieldName = self.fieldList[f]
            #print('discretized = ', discretized.shape, discretized)
            fieldVals = list(discretized[f, :])
            if len(fieldVals) == 0:
                return None
            data[fieldName] = fieldVals
        return data

    def getMidpoints(self, field):
        indx = self.fieldIndex[field]
        dSpec = self.discSpecs[indx]
        edges = dSpec[3]
        isDiscrete = dSpec[5]
        mids = []
        if isDiscrete:
            for i in range(len(edges) - 1):
                mids.append(edges[i])
        else:
            for i in range(len(edges) - 1):
                mids.append((edges[i] + edges[i+1]) / 2)
        return mids

    def getBucketVals(self, field):
        indx = self.fieldIndex[field]
        dSpec = self.discSpecs[indx]
        bucketCount = dSpec[0]
        return range(bucketCount)

    def filter(self, filtSpecs):
        #print('self.discretized = ', self.discretized)
        filtSpecs2 = []
        # Transform filtSpecs from (fieldName, value) to (fieldIndex, discretizedValue)
        for filt in filtSpecs:
            field, value = filt
            indx = self.fieldIndex[field]
            dSpec = self.discSpecs[indx]
            edges = list(dSpec[3])
            dValue = np.digitize(value, edges[:-1]) - 1
            filtSpecs2.append((indx, dValue))
            #print('dValue = ', dValue)
        remRecs = []
        for i in range(self.N):
            include = True
            for filt in filtSpecs2:
                fieldInd, dValue = filt
                if self.discretized[fieldInd, i] != dValue:
                    include = False
                    break
            if not include:
                remRecs.append(i)
        filtered = np.delete(self.aData, remRecs, 1)
        #print('filtered = ', filtered.shape,  filtered)
        # Now convert to orginal format, with only ecords that passed filter
        orgFilt = self.toOriginalForm(filtered)
        return orgFilt

    def binToVal(self, field, bin):
        indx = self.fieldIndex[field]
        dSpec = self.discSpecs[indx]
        min = dSpec[1]
        max = dSpec[2]
        val = (min + max) / 2
        return val

    def prob(self, rvName, value, givenSpec=None):
        """Return the probability that a variable is equal to a particular value
           given a set of conditions.  Returns a probability value 0 <= v <= 1.
           See distr (below) for a definition of the givenSpec.
        """
        d = self.distr(rvName, givenSpec)
        return d.P(value)

    P = prob

    def distr(self, rvName, givenSpecs = None):
        """Return a probability distribution as a PDF (see pdf.py) for the random variable
           indicated by rvName.
           If givenSpecs is provided, then will return the conditional distribution,
           otherwise will return the unconditional (i.e. marginal) distribution.
           This satisfies the following types of probability queries:
            - P(Y)
            - P(Y | X=x) -- Conditional probability
            - P(Y | X1=x1, ... ,Xk = xk) -- Multiple conditions
            - P(Y | X=x, Z) -- i.e. Conditionalize on Z
            - P(Y | X=x, Z1, ... Zk) -- Conditionalize on multiple variables
            Each spec in givenSpecs may contain a tuple or a single value:
            - (varName, value) indicates a condition
            - varName indicates conditionalizing on a variable
                (i.e. SUM(Y | X1=x1, ... , Xk = xk, Z = z) for all z in Z)
            If there is only one condition, it can be passed as s single tuple
            Examples:
            - distr('Y') -- The marginal probability of Y
            - distr('Y', [('X', 1)]) -- The probability of Y given X=1
            - distr('Y', ('X', 1)) -- The probability of Y given X=1 (same as above)
            - distr('Y', [('X1', 1), ('X2', 0)]) - The probability of Y given X1 = 1, and X2 = 0
            - distr('Y', [('X', 1), 'Z']) -- The probability of Y given X = 1, conditionalized on Z
        """
        isDiscrete = self.isDiscrete(rvName)
        if givenSpecs is None:
            # Marginal (unconditional) Probability
            indx = self.fieldIndex[rvName]
            dSpec = self.discSpecs[indx]
            hist = list(dSpec[4])
            edges = list(dSpec[3])
            minV = dSpec[1]
            maxV = dSpec[2]
            outHist = []
            for i in range(len(hist)):
                val = hist[i]
                outHist.append(val / self.N)
            pdfSpec = []
            for i in range(len(outHist)):
                start = edges[i]
                end = edges[i+1]
                pdfSpec.append((i, start, end, outHist[i]))
            outPDF = PDF(pdfSpec, isDiscrete=isDiscrete)
            return outPDF
        else:
            # Conditional Probability
            if type(givenSpecs) != type([]):
                givenSpecs = [givenSpecs]
            filtSpecs = [] # Filter specifications
            condSpecs = [] # Conditionalization specifications
            for givenSpec in givenSpecs:
                if type(givenSpec) == type((1,)):
                    # Is a conditional spec (varName, value)
                    # We filter the distribution based on the value
                    condVar, val = givenSpec
                    filtSpecs.append(givenSpec)
                else:
                    # Its a variable to conditionalize on.
                    condSpecs.append(givenSpec[0])
            newData = self.filter(filtSpecs)
            if not newData:
                # No data left after filtering.  Return a null distribution
                return None
            # Create a new probability object based on the filtered data
            filtSample = Sample(newData, density = self.density, discSpecs = self.discSpecs)
            if len(condSpecs) == 0:
                # No conditionalizing needed.  Return the pdf of the filtered data.
                outPDF = filtSample.distr(rvName)
            else:
                # Conditionalize on all indicated variables. I.e.,
                # SUM(P(filteredY | Z=z) * P(Z=z)) for all z in Z.
                accum = None
                conditionalizeOn = []
                for given in condSpecs:
                    conditionalizeOn.append(given)
                filtSpecs = self.jointCondSpecs(conditionalizeOn)
                for f in filtSpecs:
                    value = f[0]
                    given = conditionalizeOn[0]
                    probGV = self.prob(given, value) # P(Z=z)
                    if probGV == 0:
                         # Zero probability -- don't bother accumulating
                        continue
                    # probTGV is P(filteredY | Z=z) e.g., P(Y | X=1, Z=z)
                    probTGV = filtSample.distr(rvName, (given, value))
                    if not probTGV:
                        # Zero probability.  No need to accumulate
                        continue
                    probs = probTGV.ToHistogram() * probGV # Creates an array of probabilities
                    if accum is None:
                        accum = probs
                    else:
                        accum += probs
                # Now we start with a pdf of the original variable to establish the ranges, and
                # then replace the actual probabilities of each bin.  That way we maintain the
                # original bin structure. 
                template = filtSample.distr(rvName)
                outSpecs = []
                for i in range(len(accum)):
                    pdfBin = template.getBin(i)
                    newprob = accum[i]
                    newBin = pdfBin[:-1] + (newprob,)
                    outSpecs.append(newBin)
                outPDF = PDF(outSpecs, isDiscrete = isDiscrete)
            return outPDF

    PD = distr
    
    def jointValues(self, rvList):
        """ Return a list of the joint distribution values for a set of variables.
            I.e. [(rv1Val, rv2Val, ... , rvNVal)] for every combination of bin values.
        """
        nVars = len(rvList)
        rvName = rvList[0]
        vals = self.getMidpoints(rvName)
        if nVars == 1:
            return ([(val,) for val in vals])
        else:
            accum = []
            childVals = self.jointValues(rvList[1:]) # Recurse to get the child values
            for val in vals:
                accum += [(val,) + childVal for childVal in childVals]
            return accum

    def jointCondSpecs(self, rvList):
        condSpecList = []
        jointVals = self.jointValues(rvList)
        for item in jointVals:
            condSpecs = []
            for i in range(len(rvList)):
                rvName = rvList[i]
                val = item[i]
                spec = (rvName, val)
                condSpecs.append(spec)
            condSpecList.append(condSpecs)
        return condSpecList


    def pdfToProbArray(self, pdf):
        vals = []
        for bin in pdf:
            prob = bin[3]
            vals.append(prob)
        outArray = np.array(vals)
        return outArray

    def fieldStats(self, field):
        return self.fieldAggs[field]

    def _isDiscreteVar(self, rvName):
        vals = self.ds[rvName]
        isDiscrete = True
        for val in vals:
            if val != int(val):
                isDiscrete = False
                break
        return isDiscrete
    
    def _getDiscreteVars(self):
        discreteVars = []
        for var in self.fieldList:
            if self._isDiscreteVar(var):
                discreteVars.append(var)
        return discreteVars

    def isDiscrete(self, rvName):
        indx = self.fieldIndex[rvName]
        dSpec = self.discSpecs[indx]
        isDisc = dSpec[5]
        return isDisc
