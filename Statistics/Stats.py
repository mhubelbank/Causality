# Module to compute field aggregates for a dataset
import numpy as np
import math

class Stats:

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
            cardinality = len(np.unique(self.aData[i]))
            if cardinality < math.sqrt(self.N):
                nBuckets = cardinality
            elif self.N < 100:
                nBuckets = 10
            else:
                nBuckets = int(math.sqrt(self.N))
            hist, edges = np.histogram(self.aData[i], nBuckets)
            minV = np.min(self.aData[i])
            maxV = np.max(self.aData[i])

            discSpecs.append((nBuckets, minV, maxV, edges, hist))
        return discSpecs

    def fixupDiscSpecs(self, discSpecs):
        outSpecs = []
        for i in range(len(discSpecs)):
            discSpec = discSpecs[i]
            bins, min, max, edges, hist = discSpec
            # Regenerate histogram.  The other data should use the original
            newHist, newEdges = np.histogram(self.aData[i], bins, (min, max))
            #assert edges == newEdges, 'Bad Edges'
            outSpecs.append((bins, min, max, edges, newHist))
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
        mids = []
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

    def expected(self, pdf, discrete = False):
        cum = 0
        for i in range(len(pdf)):
            spec = pdf[i]
            bin, min, max, prob = spec
            if discrete:
                value = min
            else:
                value = (min + max) / 2
            cum += prob * value
        if discrete:
            exp = round(cum)
        else:
            exp = cum
        return exp
            
    def probDisc(self, target, bin):
        indx = self.fieldIndex[target]
        dSpec = self.discSpecs[indx]
        hist = list(dSpec[4])
        if bin > len(hist) - 1:
            bin = len(hist) - 1
        prob = hist[bin] / self.N
        return prob

    def prob(self, target, value):
        indx = self.fieldIndex[target]
        dSpec = self.discSpecs[indx]
        hist = list(dSpec[4])
        edges = list(dSpec[3])
        bin = np.digitize(value, edges[:-1]) - 1
        if bin > len(hist) - 1:
            bin = len(hist) - 1
        prob = hist[bin] / self.N
        return prob

    def condProb(self, target, targetVal, givenSpecs):
        if type(givenSpecs) != type([]):
            givenSpecs = [givenSpecs]
        filtSpecs = []
        condSpecs = []
        for givenSpec in givenSpecs:
            field, val = givenSpec
            # If val is None, we conditionalize.  If a specific value, we filter
            if val is not None:
                filtSpecs.append(givenSpec)
            else:
                condSpecs.append(givenSpec)
        newData = self.filter(filtSpecs)
        newStats = Stats(newData, self.discSpecs)
        #print('newPDF(B) = ', newStats.N, newStats.pdf('B'))
        return newStats.prob(target, targetVal)

    def pdf(self, target):
        indx = self.fieldIndex[target]
        dSpec = self.discSpecs[indx]
        hist = list(dSpec[4])
        edges = list(dSpec[3])
        minV = dSpec[1]
        maxV = dSpec[2]
        outHist = []
        for i in range(len(hist)):
            val = hist[i]
            outHist.append(val / self.N)
        outPDF = []
        for i in range(len(outHist)):
            start = edges[i]
            end = edges[i+1]
            outPDF.append((i, start, end, outHist[i]))
        return outPDF

    def pdfToProbArray(self, pdf):
        vals = []
        for bin in pdf:
            prob = bin[3]
            vals.append(prob)
        outArray = np.array(vals)
        return outArray

    def condPDF(self, target, givenSpecs):
        if type(givenSpecs) != type([]):
            givenSpecs = [givenSpecs]
        filtSpecs = []
        condSpecs = []
        for givenSpec in givenSpecs:
            field, val = givenSpec
            # If val is None, we conditionalize.  If a specific value, we filter
            if val is not None:
                filtSpecs.append(givenSpec)
            else:
                condSpecs.append(givenSpec[0])
        newData = self.filter(filtSpecs)
        if not newData:
            return []
        filtStats = Stats(newData, self.discSpecs)
        if len(condSpecs) == 0:
            outPDF = filtStats.pdf(target)
            #print('filtStats.discSpecs = ', filtStats.discSpecs)
        else:
            accum = None
            for given in condSpecs:
                values = self.getMidpoints(given)
                for value in values:
                    probGV = self.prob(given, value)
                    if probGV == 0:
                        continue
                    probTGV = filtStats.condPDF(target, (given, value))
                    if not probTGV:
                        continue
                    probs = self.pdfToProbArray(probTGV) * probGV
                    if accum is None:
                        accum = probs
                    else:
                        accum += probs
                break
            template = filtStats.pdf(target)
            outPDF = []
            for i in range(len(accum)):
                pdfBin = template[i]
                newprob = accum[i]
                newBin = pdfBin[:-1] + (newprob,)
                outPDF.append(newBin)
        return outPDF

    def fieldStats(self, field):
        return self.fieldAggs[field]

    def __init__(self, ds, discSpecs = None):
        self.ds = ds
        self.fieldList = list(ds.keys())
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
