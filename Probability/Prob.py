# Module to compute field aggregates for a dataset
import numpy as np
import math
from math import log, sqrt
from Probability import probCharts
from Probability.pdf import PDF

DEBUG = False

class ProbSpace:
    def __init__(self, ds, density = 1.0, power=1, discSpecs = None):
        """ Probability Space (i.e. Joint Probability Distribution) based on a a multivariate dataset
            of random variables provided.  'JPS' (Joint Probability Space) is an alias for ProbSpace.
            The Joint Probability Space is a multi-dimensional probability distribution that embeds all
            knowledge about the statistical relationships among the variables, and supports a powerful
            range of queries to expose that information.
            It can handle discrete as well as continuous variables.  Continuous probabilities
            are managed by dense discretization (binning) continuous variables into small ranges.
            By default, the number of discretized bins for continuous variables is the square
            root of the number of samples.  This can be increased or decreased using the "density"
            parameter (see below).  Discrete variables may be binary, categorical or numeric(integer).

            Data(ds) is provided as a dictionary of variable name -> ValueList (List of variable values).
            ValueList should be the same length for all variables.  Each index in ValueList represents
            one sample.

            The density parameter is used to increase or decrease the number of bins used for continuous
            variables.  If density is 1 (default), then sqrt(N-samples) bins are used.  If density is
            set to D, then D * sqrt(N-samples) bins are used.

            The power parameter is used for stochastic approximation of conditional probabilities.
            Power may range from 0 (test conditionality at mean only) up to 100 (test every comination
            of discretized variables).  Values > 0 test more and more points for conditional dependence,
            until at 100, all points are tested.  For linear relationships, power of 0 or 1 is sufficient,
            while for complex, discontinuous relationships, higher values may be necessary to achieve
            high precision.  Power allows a tradeoff between precision and run-time.  High values of
            power may result in unacceptably long run-times.  In practice, power <= 8 should suffice in
            most cases.

            The discSpecs parameter is used to make recursive calls to the module while
            maintaining the discretization information, and should not be provided by the user.

            ProbSpace includes a 'Plot' function that requires matplotlib, and produces a probability
            distribution plot for each variable.

            The main functions of ProbSpace are:
            - P(...) -- Returns the numerical probability of an event, given a set of conditions.
                P can return joint probabilities as well as univariate probabilities.
            - E(...) -- Returns the expected value (i.e. mean) of a variable, given a set of conditions.
            - distr(...) -- Returns a univariate probability distribution (see pdf.py) of a variable
                given a set of conditions.
            - dependence(...) -- Measures the dependence between two variables with optional conditioning
                on a set of other variables (i.e. conditional dependence.).
        """
        self.power = power
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
        self.probCache = {} # Probability cache
        self.distrCache = {} # Distribution cache
        self.fieldAggs = self.getAgg(ds)
        if discSpecs:
            self.discSpecs = self.fixupDiscSpecs(discSpecs)
            self.discretized = self.discretize()
        else:
            self.discSpecs = self.calcDiscSpecs()
            self.discretized = self.discretize()

    def getAgg(self, ds):
        fieldList = list(ds.keys())
        numObs = self.aData.shape[1]  # Number of observations
        if numObs > 0:
            mins = self.aData.min(1)
            maxs = self.aData.max(1)
            means = self.aData.mean(1)
            stds = self.aData.std(1)
        outDict = {}
        for i in range(self.aData.shape[0]):
            fieldName = fieldList[i]
            if numObs:
                aggs = (mins[i], maxs[i], means[i], stds[i])
            else:
                aggs = (0,0,0,0)
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
                    hist[j] = count
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
            fieldVals = list(discretized[f, :])
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

    def filteredSpace(self, filtSpec, power=None, density=None, discSpecs=None):
        """ Return a new ProbSpace object with just the data that passes the filter.
            The returned object represents the multivariate joint probability space
            of the filtered data.
        """
        if power is None:
            power = self.power
        if density is None:
            density = self.density
        #print('filtSpec = ', filtSpec)
        filtDat = self.filter(filtSpec)
        newPS = ProbSpace(filtDat, power = power, density = density, discSpecs = discSpecs)
        return newPS

    def filter(self, filtSpec):
        """ Filter the data based on a set of filterspecs: [filtSpec, ...]
            filtSpec := (varName, value) or (varName, lowValue, highValue)
            Returns a data set in the original {varName:[varData], ...}
            dictionary format.
        """
        filtdata = self.filterDat(filtSpec)
        # Now convert to orginal format, with only records that passed filter
        outData = self.toOriginalForm(filtdata)
        return outData

    def filterDat(self, filtSpec, adat = None):
        """ Filter the data in its array form and return a filtered array.
        """
        if adat is None:
            adat = self.aData
        # Transform filtSpecs from (varName, value) to (varIndex, value) or for
        # continuous variables: (varIndex, low, high)
        filtSpec2 = []
        for filt in filtSpec:
            if len(filt) == 2:
                var, val = filt
                indx = self.fieldIndex[var]
                dSpec = self.discSpecs[indx]
                filtSpec2.append((indx, val))
            else:
                field, low, high = filt
                indx = self.fieldIndex[field]
                dSpec = self.discSpecs[indx]
                edges = list(dSpec[3])
                varMin = dSpec[1]
                varMax = dSpec[2]
                if low is None:
                    low = varMin
                if high is None:
                    high = varMax
                filtSpec2.append((indx, low, high))

        remRecs = []
        for i in range(self.N):
            include = True
            for filt in filtSpec2:
                if len(filt) == 2:
                    fieldInd, targetVal = filt
                    if type(targetVal) == type((0,)):
                        targetVal = val[0]
                    val = adat[fieldInd, i]
                    if val != targetVal:
                        include = False
                        break
                else:
                    fieldInd, low, high = filt
                    val = adat[fieldInd, i]
                    if val < low or val >= high:
                        include = False
                        break
            if not include:
                remRecs.append(i)
        # Remove all the non included rows
        filtered = np.delete(adat, remRecs, 1)
        return filtered

    def binToVal(self, field, bin):
        indx = self.fieldIndex[field]
        dSpec = self.discSpecs[indx]
        min = dSpec[1]
        max = dSpec[2]
        val = (min + max) / 2
        return val

    def makeHashkey(self, targetSpec, givenSpec):
        if type(targetSpec) == type([]):
            targetSpec = tuple(targetSpec)
        if type(givenSpec) == type([]):
            givenSpec = tuple(givenSpec)
        hashKey = (targetSpec, givenSpec)
        return hashKey

    def E(self, target, givensSpec, power=None):
        """ Returns the expected value (i.e. mean) of the distribution
            of a single variable given a set of conditions.  This is
            a convenience function equivalent to:
                distr(target, givensSpec).E()

            targetSpec is a single variable name.
            givensSpec is a conditional specification (see distr below
            for format)
        """
        d = self.distr(target, givensSpec, power=power)
        return d.E()

    def prob(self, targetSpec, givenSpec=None):
        """ Return the probability of a variable or (set of variables)
            attaining a given value (or range of values) given a set
            of conditionalities on other variables.
            'P' is an alias for prob.
            The basic form is the probability of Y given X or the probability
            of targetSpec given givenSpec, where X and Y represent events or lists
            of simultaneous events:
                P(Y=y | X=x) = P(targetSpec | givenSpec)

            targetSpec (target specification) defines the result to be returned
            (i.e. the event or set of events whose probability is to be determined).
            A target specification may take one of several forms:
            - 2-tuple (varName, value) for the probability
                of attaining a single value
            - 3-tuple: (varName, minValue, maxValue) indicating an interval:
                [minValue, maxValue) (i.e. minValue <= value < maxValue).
                minValue or maxValue may be None.  A minValue of None implies 
                -infinity.  A maxValue of None implies infinity.
            - list of either of the above or any combination.  In this case, the joint 
                probability of the events is returned.
            
            givenSpec (optional) has the same format as targetSpec with equivalent meanings
            for the givens.  A given specification supports one additional flavor which
            is a single variable name.  This means that that variable should be
            "conditionalized" on.  
            The three flavors (varName, 2-tuple, 3-tuple)
            may be mixed within a givenSpec, presented as a list of givens.

            Examples:
            - prob(('A', 1)) -- The probability that variable A takes on the value 1.
            - prob([('A', 1), ('B', 2)]) -- The (joint) probability that A is 1 and B is 2.
            - prob(('A', .1, .5)) -- The probability that A is in the range [.1, .5).
            - prob(('A', .1, .5), [('B', 0, 1), ('C', -1, None)] -- The probability that
                    A is on interval [.1, .5) given that B is on interval [0,1) and
                    C is on interval [-1, infinity).
            - prob(('A', .1, .5), [('B', 0, 1), ('C', -1, None), 'D'] -- The probability that
                variable A is on interval [.1, .5) given that B is on interval [0, 1) and
                    C is on interval [-1, infinity), conditionalized on D.
    
            Conditionalization is taking the probability weighted sum of the results for every value
            of the conditionalizing variable or combination of conditionalizing variables.
            For example:
            - P(A=1 | B=2, C) is: sum over all (C=c) values( P(A=1 | B=2, C=c) * P(C=c))
        """
        if givenSpec is None:
            givenSpec = []
        if type(givenSpec) != type([]):
            givenSpec = [givenSpec]
        cacheKey = self.makeHashkey(targetSpec, givenSpec)
        if cacheKey in self.probCache.keys():
            return self.probCache[cacheKey]
        if type(targetSpec) == type([]):
            # Joint probability
            # Separate unbound (e.g. A) specifications from bound (e.g. A=1) specifications
            # Prob doesn't return unbound results.
            targetSpecU, targetSpec = self.separateSpecs(targetSpec)
            assert len(targetSpecU) == 0, 'prob.P: All target specifications must be bound (i.e. specified as tuples).  For unbound returns, use distr.)'
            result = self.jointProb(targetSpec, givenSpec)
            self.probCache[cacheKey] = result
            return result
        else:
            rvName = targetSpec[0]
            if len(targetSpec) == 2:
                valSpec = targetSpec[1]
            else:
                valSpec = targetSpec[1:]
            d = self.distr(rvName, givenSpec)
            result = d.P(valSpec)
            self.probCache[cacheKey] = result
            return result

    P = prob

    def distr(self, rvName, givenSpecs=None, power=None):
        """Return a univariate probability distribution as a PDF (see pdf.py) for the random variable
           indicated by rvName.
           If givenSpec is provided, then will return the conditional distribution,
           otherwise will return the unconditional (i.e. marginal) distribution.
           This satisfies the following types of probability queries:
            - P(Y) -- (marginal) Probability distribution of Y
            - P(Y | X=x) -- Conditional probability
            - P(Y | X1=x1, ... ,Xk = xk) -- Multiple conditions
            - P(Y | X=x, Z) -- i.e. Conditionalize on Z
            - P(Y | X=x, Z1, ... Zk) -- Conditionalize on multiple variables

            rvName is the name of the random variable whose distribution is requested.
            givenSpec (given specification) defines the conditions (givens) to
            be applied.
            A given specification may take one of several forms:
            - 2-tuple (varName, value) - Variable taking on a given value.
            - 3-tuple: (varName, minValue, maxValue) indicating an interval:
                [minValue, maxValue) (i.e. minValue <= value < maxValue).
                minValue or maxValue may be None.  A minValue of None implies 
                -infinity.  A maxValue of None implies infinity.
            - variable name: A variable to conditionalize on.
            - list of any of the above or any combination of above.

            Examples:
            - distr('Y') -- The (marginal) probability of Y
            - distr('Y', [('X', 1)]) -- The probability of Y given X=1.
            - distr('Y', [('X', 1, 2)]) -- The probability of Y given 1 <= X < 2.
            - distr('Y', ('X', 1)) -- The probability of Y given X=1 (same as above)
            - distr('Y', [('X1', 1), ('X2', 0)]) - The probability of Y given X1 = 1, and X2 = 0
            - distr('Y', [('X', 1), 'Z']) -- The probability of Y given X = 1, conditionalized on Z

            Conditionalization is taking the probability weighted sum of the results for every value
            of the conditionalizing variable or combination of conditionalizing variables.
            For example:
            - P(Y | X=1, Z) is: sum over all (Z=z) values( P(Y | X=1, Z=z) * P(Z=z))
        """
        if DEBUG:
            print('Prob.Sample: P(' , rvName, '|', givenSpecs , ')')
        if power is None:
            power = self.power
        if givenSpecs is None:
            givenSpecs = []
        if type(givenSpecs) != type([]):
            givenSpecs = [givenSpecs]
        cacheKey = self.makeHashkey(rvName, givenSpecs)
        if cacheKey in self.distrCache.keys():
             return self.distrCache[cacheKey]
        isDiscrete = self.isDiscrete(rvName)
        indx = self.fieldIndex[rvName]
        dSpec = self.discSpecs[indx]
        bins = dSpec[0]

        if not givenSpecs:
            # Marginal (unconditional) Probability
            bins = dSpec[0]
            hist = list(dSpec[4])
            if not hist:
                hist = [0] * bins
            edges = list(dSpec[3])
            minV = dSpec[1]
            maxV = dSpec[2]
            outHist = []
            for i in range(len(hist)):
                cnt = hist[i]
                if self.N > 0:
                    outHist.append(cnt / self.N)
                else:
                    outHist.append(0)
            pdfSpec = []
            for i in range(len(outHist)):
                start = edges[i]
                end = edges[i+1]
                pdfSpec.append((i, start, end, outHist[i]))
            outPDF = PDF(self.N, pdfSpec, isDiscrete=isDiscrete)
            self.distrCache[cacheKey] = outPDF
            return outPDF
        else:
            # Conditional Probability
            condSpecs, filtSpecs = self.separateSpecs(givenSpecs)
            if filtSpecs:
                # Create a new probability object based on the filtered data
                filtSpace = self.filteredSpace(filtSpecs, density = self.density, power = self.power, discSpecs=self.discSpecs)
            else:
                filtSpace = self
            if not condSpecs:
                # Nothing to conditionalize on.  We can just return the selected variable
                # from the filtered distribution
                outPDF = filtSpace.distr(rvName)
            else:
                # Conditionalize on all indicated variables. I.e.,
                # SUM(P(filteredY | Z=z) * P(Z=z)) for all z in Z.
                accum = np.zeros((bins,))
                conditionalizeOn = []
                for given in condSpecs:
                    conditionalizeOn.append(given)
                condFiltSpecs = self.getCondSpecs(conditionalizeOn, power=power)
                #print('condFiltSpecs = ', condFiltSpecs)
                #countRatio = float(self.N) / filtSample.N
                allProbs = 0.0 # The fraction of the probability space that has been tested.
                for cf in condFiltSpecs:
                    probZ = self.jointProb(cf) # P(Z=z)
                    if probZ == 0:
                        #print('probZ = 0')
                        # Zero probability -- don't bother accumulating
                        continue
                    # probYgZ is P(Y | Z=z) e.g., P(Y | X=1, Z=z)
                    probYgZ = filtSpace.distr(rvName, cf)
                    if probYgZ:
                        probs = probYgZ.ToHistogram() * probZ # Creates an array of probabilities
                        accum += probs
                    allProbs += probZ
                    #print('distr: cf = ', cf, ', probs = ', probYgZ.stats(), ', probZ = ', probZ)
                accum = accum / allProbs
                N = self.distr(rvName, filtSpecs).N
                # Now we start with a pdf of the original variable to establish the ranges, and
                # then replace the actual probabilities of each bin.  That way we maintain the
                # original bin structure. 
                template = self.distr(rvName)
                outSpecs = []
                for i in range(len(accum)):
                    pdfBin = template.getBin(i)
                    newprob = accum[i]
                    newBin = pdfBin[:-1] + (newprob,)
                    outSpecs.append(newBin)
                outPDF = PDF(N, outSpecs, isDiscrete = isDiscrete)
            self.distrCache[cacheKey] = outPDF
            return outPDF

    PD = distr

    def getCondSpecs(self, condVars, power=2, deltaAdjust=1):
        rawCS = self.getCondSpecs2(condVars, power = power, deltaAdjust = deltaAdjust)
        #print('rawCS = ', rawCS)
        outCS = []
        for spec in rawCS:
            #print('spec = ', spec)
            currPS = self
            outSpec = []
            # Adjust pseudo filters by the mean and std of the conditional
            for s in range(len(spec)):
                varSpec = spec[s]
                #print('varSpec = ', varSpec)
                var = varSpec[0]
                #print('var = ', var)
                if self.isDiscrete(var):
                    outSpec.append(varSpec)
                else:
                    min, max = varSpec[1:]
                    distr = currPS.distr(var)
                    mean = distr.E()
                    std = distr.stDev()
                    #print('mean, std = ', mean, std)
                    varSpec = (var, mean + min * std, mean + max * std)
                    outSpec.append(varSpec)
                # Use resulting conditional space of this variable for the sample
                # of the next one.  That way, we sample using the mean and std
                # of the variable in the conditionaed space of the previous vars.
                if s != len(spec) - 1:
                    currPS = currPS.filteredSpace([varSpec])
            #print('outSpec = ', outSpec)
            outCS.append(outSpec)
        #print('outCS = ', outCS)
        return outCS

    def getCondSpecs2(self, condVars, power=2, deltaAdjust=1):
        """ Produce a set of conditional specifications for stochastic
            conditionalization, given
            a set of variables to conditionalize on, and a desired power level.
            Power determines how many points to use to conditionalize on.
            Zero indicates conditionalize on the mean alone. 1 uses the mean
            and two other points (one on either side of the mean).
            2 Uses the mean plus 4 other points (2 on each side of the mean).
            Power values (p) less than 100 will test p * 2 + 1 values for each
            variable.range
            power value > 100 indicates that all values will be tested,
            which can be extremely processor intensive.
            Conditional specifications provide a list of lists of tuple:
            [[(varName1, value1_1), (varName2, value2_1), ... (varNameK, valueK_1)],
             [(varName1, value1_2), (varName2, value2_2), ... (varNameK, valueK_2)],
             ...
             [(varName1, value1_N), (varName2, value2_N), ... (varNameK, valueK_N)]]
            Where K is the number of conditional variables, and N is the total number
            of combinations = K**(2 * P + 1) for values of P < 100.
            """
        def getTestVals(self, rv, deltaAdjust):
            deltaBase = .05
            delta = deltaBase * deltaAdjust
            isDiscrete = self.isDiscrete(rvName)
            if isDiscrete or testPoints is None:
                # If is Discrete, return all values
                testVals = self.getMidpoints(rv)
                testVals = [testVal for testVal in testVals]
            else:
                # If continuous, sample values at testPoint distances from the mean
                testVals = []
                for tp in testPoints:
                    if tp == 0:
                        # For 0, just use the mean
                        testVals.append((- delta, + delta))
                    else:
                        # For nonzero, test points mean + tp and mean - tp
                        testVals.append((-tp - delta, -tp + delta))
                        testVals.append((tp - delta, tp + delta))
            return testVals
        levelSpecs0 = [.5, .75, .25, 1.0, 1.5, .25, .75, 1.25, 1.75, 2.0]
        maxLevel = 3 # Largest standard deviation to sample
        if power <= 10:
            levelSpecs = levelSpecs0
        elif power < 100:
            levelSpecs = np.range(1/power, maxLevel + 1/power, 1/power)
        else:
            levelSpecs = None
        if levelSpecs:
            testPoints = [0] + levelSpecs[:power]
        else:
            # TestPoints None means test all values
            testPoints = None

        # Find values for each variable based on testPoints
        nVars = len(condVars)
        rvName = condVars[0]
        if nVars == 1:
            # Only one var to do.  Find the values.
            vals = getTestVals(self, rvName, deltaAdjust=deltaAdjust)
            if self.isDiscrete(rvName):
                return [[(rvName, val)] for val in vals]
            else:
                return [[((rvName,) + val)] for val in vals]
        else:
            # We're not on the last var, so recurse and build up the total set
            accum = []
            vals = getTestVals(self, rvName, deltaAdjust=deltaAdjust)
            childVals = self.getCondSpecs(condVars[1:], power) # Recurse to get the child values
            for val in vals:
                if self.isDiscrete(rvName):
                    accum += [[(rvName, val)] + childVal for childVal in childVals]
                else:
                    accum += [[(rvName,) + val] + childVal for childVal in childVals]
            return accum
        # End of getCondSpecs
        
    def getCondSpecs_orig(self, condVars, power=2, deltaAdjust=1):
        """ Produce a set of conditional specifications for stochastic
            conditionalization, given
            a set of variables to conditionalize on, and a desired power level.
            Power determines how many points to use to conditionalize on.
            Zero indicates conditionalize on the mean alone. 1 uses the mean
            and two other points (one on either side of the mean).
            2 Uses the mean plus 4 other points (2 on each side of the mean).
            Power values (p) less than 100 will test p * 2 + 1 values for each
            variable.range
            power value > 100 indicates that all values will be tested,
            which can be extremely processor intensive.
            Conditional specifications provide a list of lists of tuple:
            [[(varName1, value1_1), (varName2, value2_1), ... (varNameK, valueK_1)],
             [(varName1, value1_2), (varName2, value2_2), ... (varNameK, valueK_2)],
             ...
             [(varName1, value1_N), (varName2, value2_N), ... (varNameK, valueK_N)]]
            Where K is the number of conditional variables, and N is the total number
            of combinations = K**(2 * P + 1) for values of P < 100.
            """
        def getTestVals(self, rv, deltaAdjust):
            deltaBase = .05
            delta = deltaBase * deltaAdjust
            isDiscrete = self.isDiscrete(rvName)
            if isDiscrete or testPoints is None:
                # If is Discrete, return all values
                testVals = self.getMidpoints(rv)
                testVals = [(testVal,) for testVal in testVals]
            else:
                # If continuous, sample values at testPoint distances from the mean
                testVals = []
                p = self.distr(rv)
                mean = p.E()
                std = p.stDev()
                for tp in testPoints:
                    if tp == 0:
                        # For 0, just use the mean
                        testVals.append((mean - delta * std, mean + delta * std))
                    else:
                        # For nonzero, test points mean + tp and mean - tp
                        testVals.append((mean - tp * std - delta * std, mean - tp * std + delta * std))
                        testVals.append((mean + tp * std - delta * std, mean + tp * std + delta * std))
            return testVals
        levelSpecs0 = [.5, .75, .25, 1.0, 1.5, .25, .75, 1.25, 1.75, 2.0]
        maxLevel = 3 # Largest standard deviation to sample
        if power <= 10:
            levelSpecs = levelSpecs0
        elif power < 100:
            levelSpecs = np.range(1/power, maxLevel + 1/power, 1/power)
        else:
            levelSpecs = None
        if levelSpecs:
            testPoints = [0] + levelSpecs[:power]
        else:
            # TestPoints None means test all values
            testPoints = None

        # Find values for each variable based on testPoints
        nVars = len(condVars)
        rvName = condVars[0]
        if nVars == 1:
            # Only one var to do.  Find the values.
            vals = getTestVals(self, rvName, deltaAdjust=deltaAdjust)
            if self.isDiscrete(rvName):
                return [[(rvName, val)] for val in vals]
            else:
                return [[((rvName,) + val)] for val in vals]
        else:
            # We're not on the last var, so recurse and build up the total set
            accum = []
            vals = getTestVals(self, rvName, deltaAdjust=deltaAdjust)
            childVals = self.getCondSpecs(condVars[1:], power) # Recurse to get the child values
            for val in vals:
                if self.isDiscrete(rvName):
                    accum += [[(rvName, val)] + childVal for childVal in childVals]
                else:
                    accum += [[(rvName,) + val] + childVal for childVal in childVals]
            return accum
        # End of getCondSpecs_orig

    def adjustSpec(self, spec, delta):
        outSpec = []
        for var in spec:
            if len(var) == 2:
                # Discrete.  Don't modify
                outSpec.append(var)
            else:
                varName, low, high = var
                mid = (low + high) / 2.0
                oldDelta = (high - low) / 2.0
                newDelta = delta * oldDelta
                newVar = (varName, mid - newDelta, mid + newDelta)
                #print('oldVar = ', var, 'newVar = ', newVar, mid, oldDelta, newDelta, deltaAdjust)
                outSpec.append(newVar)
        #print('old = ', spec, ', new =', outSpec, ', delta = ', deltaAdjust)
        return outSpec

    def dependence_new(self, rv1, rv2, givensSpec=[], power = None):
        """ givens is [given1, given2, ... , givenN]
        """
        if power is None:
            power = self.power
        if givensSpec is None:
            givensSpec = []
        if type(givensSpec) != type([]):
            givensSpec = [givensSpec]
        maxAttempts = 5
        minVals = 100
        maxVals = int(self.N**.5) * 5
        accum = 0.0
        accumProb = 0.0
        # Get all the combinations of rv1, rv2, and any givens
        # Depending on power, we test more combinations.  If level >= 100, we test all combos
        # For level = 0, we just test the mean.  For 1, we test the mean + 2 more values.
        # For level = 3, we test the mean + 6 more values.

        # Separate the givens into bound (e.g. B=1, 1 <= B < 2) and unbound (e.g., B) specifications.
        givensU, givensB = self.separateSpecs(givensSpec)
        condFiltSpecs = self.getCondSpecs(givensU + [rv2], power)
        prevGivens = None
        prevProb1 = None
        numTests = 0
        for spec in condFiltSpecs:
            # Compare P(Y | Z) with P(Y | X,Z)
            # givens is conditional on spec without rv2
            #print('spec = ', spec)
            deltaAdjust = 1.0
            for attempt in range(maxAttempts):
                givens = spec[:-1]
                if givens != prevGivens:
                    # Only recompute prob 1 when givensValues change
                    if givens:
                        prob1 = self.distr(rv1, givens + givensB)
                    else:
                        prob1 = self.distr(rv1, givensB)
                    testRanges0 = self.getCondSpecs([rv1], power=4)
                    testRanges = [tv[0][1:] for tv in testRanges0]
                    prevProb1 = prob1
                    prevGivens = givens
                else:
                    # Otherwise use the previously computed prob1
                    prob1 = prevProb1
                prob2 = self.distr(rv1, spec + givensB)
                #print('prob1.E() = ', prob1.E(), prob1.stDev(), prob1.N)
                #print('prob2.E() = ', prob2.E(), prob2.stDev(), prob2.N)
                #print('prob2.N = ', prob2.N)
                N = prob2.N
                oldDa = deltaAdjust
                #print('Attempt = ', attempt + 1, N, maxVals, ', spec = ', spec, ', deltaAdj = ', deltaAdjust)
                if attempt != maxAttempts - 1 and (N < minVals or N > maxVals):
                    if N == 0:
                        deltaAdjust *= 5
                    elif N < minVals:
                        deltaAdjust *= min([minVals * 1.1 / float(N), 2.2])
                    else:
                        # N > maxVals
                        deltaAdjust *= max([maxVals * .9 / float(N), .5])
                    spec = self.adjustSpec(spec, deltaAdjust / oldDa)
                    #print('Too few datapoints for: ', spec, prob1.N, prob2.N)
                    continue
                else:
                    break
            deps = []
            testsComplete = 0
            allProbs = 0.0
            if not self.isDiscrete(rv1):
                for testRange in testRanges:
                    p2Test = prob2.P(testRange)
                    #if p2Test == 0:
                    #    continue
                    p1Test = prob1.P(testRange)
                    err = ((p1Test - p2Test)**2 / (p1Test + p2Test) * 2 * p1Test) if p1Test + p2Test > 0 else 0.0
                    #err = (p1Test - p2Test)**2
                    #print('p1Test = ', p1Test, ', p2Test = ', p2Test, ', err = ', err)
                    allProbs += p1Test
                    testsComplete += 1
                    deps.append(err)
                # Rescale deps by the total probability measured
                #print('deps = ', deps)
                if allProbs > 0:
                    deps = [dep / allProbs for dep in deps]
            else:
                err = prob1.compare(prob2)
                deps.append(err)
            #dep = max(deps)
            dep = sum(deps)
            #print('dep = ', dep)
            # We accumulate any measured dependency multiplied by the probability of the conditional
            # clause.  This way, we weight the dependency by the frequency of the event.
            condProb = self.jointProb(spec)
            accum += dep * condProb
            accumProb += condProb # The total probability space assessed
            numTests += 1
            if DEBUG:
                print('spec = ', spec, ', givensB = ', givensB, ', dep = ', deps, prob2.N)
        if accumProb > 0.0:
            # Normalize the results for the probability space sampled by dividing by accumProb
            dependence = accum / accumProb
            return dependence
            H = 0.03
            L = 0.001
            if dependence < L:
                calDep = dependence / (2*L)
            else:
                calDep = (dependence - L) / (H-L) / 2 + .5
            #return dependence * (.5 / 0.000128)
            return calDep

        print('Cond distr too small: ', rv1, rv2, givens)
        return 0.0

    def dependence(self, rv1, rv2, givensSpec=[], power = None):
        """ givens is [given1, given2, ... , givenN]
        """
        if power is None:
            power = self.power
        if givensSpec is None:
            givensSpec = []
        if type(givensSpec) != type([]):
            givensSpec = [givensSpec]
        maxAttempts = 8
        minVals = 100
        maxVals = int(self.N**.5) * 5
        accum = 0.0
        accumProb = 0.0
        # Get all the combinations of rv1, rv2, and any givens
        # Depending on power, we test more combinations.  If level >= 100, we test all combos
        # For level = 0, we just test the mean.  For 1, we test the mean + 2 more values.
        # For level = 3, we test the mean + 6 more values.

        # Separate the givens into bound (e.g. B=1, 1 <= B < 2) and unbound (e.g., B) specifications.
        givensU, givensB = self.separateSpecs(givensSpec)
        condFiltSpecs = self.getCondSpecs(givensU + [rv2], power)
        prevGivens = None
        prevProb1 = None
        numTests = 0
        for spec in condFiltSpecs:
            allDiscrete = True
            for varSpec in spec:
                var = varSpec[0]
                if not self.isDiscrete(var):
                    allDiscrete = False
                    break
            # Compare P(Y | Z) with P(Y | X,Z)
            # givens is conditional on spec without rv2
            #print('spec = ', spec)
            deltaAdjust = 1.0
            if allDiscrete:
                attempts = 1
            else:
                attempts = maxAttempts
            for attempt in range(attempts):
                givens = spec[:-1]
                if givens != prevGivens:
                    # Only recompute prob 1 when givensValues change
                    if givens:
                        prob1 = self.distr(rv1, givens + givensB)
                    else:
                        prob1 = self.distr(rv1, givensB)
                    prevProb1 = prob1
                    prevGivens = givens
                else:
                    # Otherwise use the previously computed prob1
                    prob1 = prevProb1
                prob2 = self.distr(rv1, spec + givensB)
                N = prob2.N
                oldDa = deltaAdjust
                #print('Attempt = ', attempt + 1, N, maxVals, ', spec = ', spec, ', deltaAdj = ', deltaAdjust)
                if attempt != attempts - 1 and (N < minVals or N > maxVals):
                    if N == 0:
                        deltaAdjust *= 5
                    elif N < minVals:
                        deltaAdjust *= min([minVals * 1.1 / float(N), 1 + 1.2 / ((attempt+1)**.5)])
                    else:
                        # N > maxVals
                        deltaAdjust *= max([maxVals * .9 / float(N), 1 - .5 / ((attempt+1)**.5)])
                    spec = self.adjustSpec(spec, deltaAdjust / oldDa)
                    continue
                else:
                    break
            deps = []
            if prob2.N == 0:
                continue
            if not self.isDiscrete(rv1):
                mean1 = prob1.E()
                mean2 = prob2.E()
                dep1 = abs((mean1 - mean2))
                std1 = prob1.stDev()
                std2 = prob2.stDev()
                if std1 + std2 > 0:
                    dep2 = abs((std1 - std2) / (std1 + std2))
                    #dep2 = abs(std1 - std2)
                else:
                    dep2 = 0.0
                sk1 = prob1.skew()
                sk2 = prob2.skew()
                ku1 = prob1.kurtosis()
                ku2 = prob2.kurtosis()
                dep3 = (sk1 - sk2)**2
                dep4 = (ku1 - ku2)**2
                dep3 = 0
                dep4 = 0
            else:
                dep1 = prob1.compare(prob2)
                #print('prob1.stats, prob2.stats = ', prob1.stats(), prob2.stats())
                dep2 = dep3 = dep4 = 0.0
            dep = max([dep1, dep2, dep3, dep4])
            #print('dep = ', dep)
            # We accumulate any measured dependency multiplied by the probability of the conditional
            # clause.  This way, we weight the dependency by the frequency of the event.
            condProb = self.jointProb(spec)
            accum += dep * condProb
            accumProb += condProb # The total probability space assessed
            numTests += 1
            if DEBUG:
                print('spec = ', spec, ', givensB = ', givensB, ', dep = ', deps, prob2.N)
        if accumProb > 0.0:
            # Normalize the results for the probability space sampled by dividing by accumProb
            dependence = accum / accumProb
            #return dependence
            H = .95
            L = .05
            if dependence < L:
                calDep = dependence / (2*L)
            else:
                calDep = (dependence - L) / (H-L) / 2 + .5
            #return dependence * (.5 / 0.000128)
            return calDep

        print('Cond distr too small: ', rv1, rv2, givens)
        return 0.0

    def separateSpecs(self, specs):
        """ Separate bound and unbound variable specs,
            and return (unboundSpecs, boundSpecs).  While we're
            at it, fixup any exact value specs for continuous
            variables to make a small range of std deviations
            around the valeu.
        """
        delta = .05
        uSpecs = []
        bSpecs = []
        for spec in specs:
            if type(spec) == type((0,)):
                # It is a bound spec
                var = spec[0]
                if len(spec) == 2 and not self.isDiscrete(var):
                    # Continuous variable with a single value.  Change
                    # to a fixed small range around that value.
                    val = spec[1]
                    std = self.distr(var).stDev()
                    bSpecs.append((var, val - delta*std, val + delta*std))
                else:
                    bSpecs.append(spec)
            else:
                # Unbound
                uSpecs.append(spec)
        return uSpecs, bSpecs

    def dependence_old(self, rv1, rv2, givensSpec=None, power=None):
        """ Calculate the dependence between two random variables, optionally
            given a set of bound or unbound other variables.
        """
        if power is None:
            power = self.power
        if givensSpec is None:
            givensSpec = []
        if type(givensSpec) != type([]):
            givensSpec = [givensSpec]
        # Separate the givens into bound (e.g. B=1, 1 <= B < 2) and unbound (e.g., B) specifications.
        givensU, givensB = self.separateSpecs(givensSpec)
        accum = 0.0
        accumProb = 0.0
        # Get all the combinations of rv1, rv2, and any givens
        # Depending on level, we test more combinations.  If level >= 100, we test all combos
        # For level = 0, we just test the mean.  For 1, we test the mean + 2 more values.
        # For level = 3, we test the mean + 6 more values.
        condFiltSpecs = self.getCondSpecs(givensU + [rv2], power)
        prevGivens = None
        prevProb1 = None
        numTests = 0
        for spec in condFiltSpecs:
            # spec is rv2 + any other unbound conditions to test.
            # Compare P(Y | Z) with P(Y | X,Z)
            # givens is all unbound conditions without rv2 (i.e. X)
            givens = spec[:-1]
            if givens != prevGivens:
                # Only recompute prob 1 when givensValues change
                if givens:
                    prob1 = self.distr(rv1, givens + givensB)
                else:
                    prob1 = self.distr(rv1, givensB)
            else:
                # Otherwise use the previously computed prob1
                prob1 = prevProb1
            prob2 = self.distr(rv1, spec + givensB)
            if prob2.N < 5:
                continue
            prevGivens = givens
            prevProb1 = prob1
            mean1 = prob1.E()
            mean2 = prob2.E()
            dep1 = abs((mean1 - mean2))
            std1 = prob1.stDev()
            std2 = prob2.stDev()
            if std1 + std2 > 0:
                dep2 = abs((std1 - std2) / (std1 + std2))
                dep2 = abs(std1 - std2)
            else:
                dep2 = 0.0
            if not self.isDiscrete(rv1):
                sk1 = prob1.skew()
                sk2 = prob2.skew()
                ku1 = prob1.kurtosis()
                ku2 = prob2.kurtosis()
                dep3 = abs((sk1 - sk2))
                dep4 = abs((ku1 - ku2))
                dep3 = 0
                dep4 = 0
            else:
                dep1 = prob1.compare(prob2)
                dep2 = dep3 = dep4 = 0.0
            dep = max([dep1, dep2, dep3, dep4])
            # We accumulate any measured dependency multiplied by the probability of the conditional
            # clause.  This way, we weight the dependency by the frequency of the event.
            if DEBUG:
                print('spec = ', spec, ', givensB = ', givensB, ', dep = ', dep, dep1, dep2, dep3, dep4, prob2.N)
            condProb = self.jointProb(spec + givensB)
            accum += dep * condProb
            accumProb += condProb # The total probability space assessed
            numTests += 1
        if accumProb > 0.0:
            # Normalize the results for the probability space sampled by dividing by accumProb
            dependence = accum / accumProb
            return dependence
        else:
            print('Cond distr too small.  Accuracy may be impaired: ', rv1, rv2, givensSpec)
            return 0.0

    def independence(self, rv1, rv2, givensSpec=None, power=None):
        """
            Calculate the independence between two variables, and an optional set of givens.
            This is a heuristic inversion
            of the dependence calculation to match other independence measures which return
            the likelihood of the null hypothesis that the variables are dependent.
            A threshold of .1 is generally used.  Values below that are considered dependent.
            givens are formatted the same as for prob(...).
            TO DO: Calibrate to an exact p-value.
        """
        dep = self.dependence(rv1, rv2, givensSpec=givensSpec, power=power)
        # Bound it to [0, 1]
        dep = max([min([dep, 1]), 0])
        ind = 1 - dep
        return ind
        if dep > .15:
            # .15 is the experimental best threshold for separating dependence from independence.
            # Invert the dependence in interval (.15, 1.0] and map it onto [0, .1]
            ind = (1-dep) * (.1 / .85)
        else:
            # Invert the dependence in interval [0, .15] and map it onto (.1, 1]
            ind = 1 - dep * (.9 / .15)
        return ind

    def isIndependent(self, rv1, rv2, givensSpec=None, power=None):
        """ Determines if two variables are independent, optionally given a set of givens.
            Returns True if independent, otherwise False
        """
        ind = self.independence(rv1, rv2, givensSpec = givensSpec, power = power)
        # Use .1 (90% confidence as threshold.
        return ind > .5

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

    def jointProb(self, varSpecs, givenSpecs=None):
        """ Return the joint probability given a set of variables and their
            values.  varSpecs is of the form (varName, varVal).  We want
            to find the probability of all of the named variables having
            the designated value.
            Join Probability is calculated as e.g.:
            - P(A, B, C) = P(A | B,C) * P(B | C) * P(C)
        """
        if givenSpecs is None:
            givenSpecs = []
        accum = []
        nSpecs = len(varSpecs)
        for i in range(nSpecs):
            spec = varSpecs[i]
            if i == nSpecs - 1:
                accum.append(self.prob(spec, givenSpecs))
            else:
                nextSpecs = varSpecs[i+1:]
                accum.append(self.prob(spec, nextSpecs + givenSpecs))

        # Return the product of the accumulated probabilities
        allProbs = np.array(accum)
        jointProb = float(np.prod(allProbs))
        return jointProb

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

    def corrCoef(self, rv1, rv2):
        """Pearson Correlation Coefficient (rho)
        """
        indx1 = self.fieldIndex[rv1]
        indx2 = self.fieldIndex[rv2]
        dat1 = self.aData[indx1,:]
        dat2 = self.aData[indx2,:]
        mean1 = dat1.mean()
        mean2 = dat2.mean()
        num1 = 0.0
        denom1 = 0.0
        denom2 = 0.0
        for i in range(self.N):
            v1 = dat1[i]
            v2 = dat2[i]
            diff1 = v1 - mean1
            diff2 = v2 - mean2
            num1 += diff1 * diff2
            denom1 += diff1**2
            denom2 += diff2**2
        rho = num1 / (denom1**.5 * denom2**.5)
        return rho

    def dat2ps(self, adat):
        ds = self.toOriginalForm(adat)
        jds = ProbSpace(ds, density = self.density, discSpecs = self.discSpecs, power=self.power)
        return jds

    def distill(self, filters, minPoints):
        """
            Progressively filter the data set until at least minPoints
            records (observations) remain, or all filters have been used.
            This is used to filter in the case where filters may be overspecified.
        """
        adat = self.aData
        for f in range(len(filters)):
            filt = filters[f]
            filtDat = self.filterDat2([filt], adat = adat)
            sh = filtDat.shape
            datLen = sh[1]
            #if datLen < 1:
            #    print('eliminating: ', filt)
            #    break
            adat = filtDat
        return adat

    def filterDat2(self, filtSpec, adat = None):
        if adat is None:
            adat = self.aData
        remRecs = []
        N = adat.shape[1]
        for i in range(N):
            include = True
            for filt in filtSpec:
                if len(filt) == 2:
                    field, val = filt
                    fieldInd = self.fieldIndex[field]
                    if adat[fieldInd, i] != val:
                        include = False
                        break
                else:
                    field, valL, valH = filt
                    fieldInd = self.fieldIndex[field]
                    fieldVal = adat[fieldInd, i]
                    if fieldVal < valL or fieldVal >= valH:
                        include = False
                        break
            if not include:
                remRecs.append(i)
        # Remove all the non included rows
        filtered = np.delete(adat, remRecs, 1)
        return filtered

    def Predict(self, Y, X):
        """
            Y is a single variable name.  X is a dataset.
        """
        dists = self.PredictDist(Y, X)
        preds = [dist.E() for dist in dists]
        return preds

    def Classify(self, Y, X):
        """
            Y is a single variable name.  X is a dataset.
        """
        assert self.isDiscrete(Y), 'Prob.Classify: Target variable must be discrete.'
        dists = self.PredictDist(Y, X)
        preds = [dist.mode() for dist in dists]
        return preds

    def PredictDist(self, Y, X):
        """
            Y is a single variable name.  X is a dataset.
        """
        delta = .05 # small range of standard deviations around the variable
        minVals = 100
        maxVals = int(self.N**.5) * 5
        maxAttempts = 8
        outPreds = []
        # Make sure Y is not in X
        vars = list(X.keys())
        try:
            vars.remove(Y)
        except:
            pass
        # Sort the independent variables by dependence with Y
        deps = [(self.dependence(var, Y, power=3), var) for var in vars]
        deps.sort()
        deps.reverse()
        vars = []
        # Remove any independent independents
        for i in range(len(deps)):
            dep = deps[i]
            if dep[0] < .5:
                print('Prob.PredictDist: rejecting variables due to independence from target(p-value, var): ', deps[i:])
                break
            else:
                vars.append(dep[1])
        #print('vars = ', vars)
        # Vars now contains all the variables that are not indpendent from Y, sorted in
        # order of highest dependence.

        # Get the standard deviation for each variable
        sigmas = {}
        for var in vars:
            aggs = self.fieldAggs[var]
            sigma = aggs[3]
            sigmas[var] = sigma

        # Get the number of items to predict:
        numTests = len(X[vars[0]])
        allDiscrete = True
        for var in vars:
            if not self.isDiscrete(var):
                allDiscrete = False
                break
        if allDiscrete:
            attempts = 1
        else:
            attempts = maxAttempts
        for i in range(numTests):
            deltaAdjust = 1.0
            for attempt in range(attempts):
                delta2 = deltaAdjust * delta
                filts = []
                for var in vars:
                    sigma = sigmas[var]
                    #print('sigma = ', sigma)
                    val = X[var][i]
                    if self.isDiscrete(var):
                        filts.append((var, val))
                    else:
                        # For continuous, use a small range around the value
                        filts.append((var, val - delta2 * sigma, val + delta2 * sigma))
                # Determine the minimum and maximum number of points to use for filtering the distribution
                #print('#', i, ', filts = ', filts)
                adat = self.distill(filts, minVals)
                jps = self.dat2ps(adat)
                dist = jps.distr(Y)
                N = dist.N
                #print('N = ', N)
                if attempt != attempts - 1 and (N < minVals or N > maxVals):
                    if N == 0:
                        deltaAdjust *= 5
                    elif N < minVals:
                        deltaAdjust *= min([minVals * 1.1 / float(N), 1 + 1.2 / ((attempt+1)**.5)])
                    else:
                        # N > maxVals
                        deltaAdjust *= max([maxVals * .9 / float(N), 1 - .5 / ((attempt+1)**.5)])
                        continue
                else:
                    outPreds.append(dist)
                    break
        return outPreds

    def Plot(self):
        """ Plot the distribution of each variable in the joint probability space
            using matplotlib.
        """
        inf = 10**30
        plotDict = {}
        minX = inf
        maxX = -inf
        numPts = 1000
        pdfs = []
        for v in self.fieldList:
            d = self.distr(v)
            minval = d.minVal()
            maxval = d.maxVal()
            if maxval > maxX:
                maxX = maxval
            if minval < minX:
                minX = minval
            pdfs.append(d)
        xvals = []
        for i in range(numPts):
            rangex = maxX - minX
            incrx = rangex / numPts
            xvals.append(minX + i * incrx)
        for i in range(len(self.fieldList)):
            yvals = []
            var = self.fieldList[i]
            pdf = pdfs[i]
            for j in range(numPts):
                xval = xvals[j]
                if j == numPts - 1:
                    P = pdf.P((xval, maxX))
                else:
                    P = pdf.P((xval, xvals[j+1]))
                yvals.append(P)
            plotDict[var] = yvals
        plotDict['_x_'] = xvals
        probCharts.plot(plotDict)

