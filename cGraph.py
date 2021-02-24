import networkx
import math
import independence
from Probability import Prob
from Probability import pdf
import numpy as np
import direction

VERBOSE = 1
# rvList is a list of random variable (rv) objects
# data is a dictionary keyed by random variable name, containing a list of observed values.
#   Each variable's list should be the same length
class cGraph:
    def __init__(self, rvList, data={}):
        self.g = networkx.DiGraph()
        self.rvDict = {}
        for rv in rvList:
            if rv.name in self.rvDict.keys():
                raise 'Duplicate variable name = ' + rv.name
            self.rvDict[rv.name] = rv
        self.rvList = list(self.rvDict.keys())
        self.rvList.sort()
        self.g.add_nodes_from(self.rvDict.keys())
        edges = []
        for rv in rvList:
            for pa in rv.parentNames:
                edges.append((pa, rv.name))
        self.g.add_edges_from(edges)
        self.data = data
        edges = self.g.edges()
        self.edgeDict = {}
        for edge in edges:
            s, d = edge
            if s in self.edgeDict.keys():
                self.edgeDict[s].append(edge)
            else:
                self.edgeDict[s] = [edge]
            if d in self.edgeDict.keys():
                self.edgeDict[d].append(edge)
            else:
                self.edgeDict[d] = [edge]
        # Create a probability sample object for later use
        self.prob = Prob.Sample(self.data, density = .2, precision=2)
        for rvName in self.rvList:
            rv = self.rvDict[rvName]
            if rv.isObserved:
                self.vPrint(rvName)

    def vPrint(self, varName):
        """ Print statistical information about a variable"""
        d = self.prob.distr(varName)
        print('stats(', varName,':', 'm1 =', d.E(), ', m2 =', d.stDev(), ', m3 =',  d.skew(), ', m4 =', d.kurtosis())
        
    def isExogenous(self, varName):
        rv = self.rvDict[varName]
        return not rv.parentNames

    def printGraph(self):
        print('Nodes:', self.g.nodes())
        print('Edges:', self.g.edges())


    def getAdjacencies(self, node):
        return self.edgeDict[node]

    def getParents(self, node):
        parents = []
        adj = self.getAdjacencies(node)
        for a in adj:
            if a[1] == node:
                parents.append(a[0])
        return parents

    def isChild(self, parentNode, childNode):
        if parentNode in self.getParents(childNode):
            return True
        return False

    def isDescendant(self, ancestor, descendant):
        parents = self.getParents(descendant)
        if ancestor in parents:
            return True
        else:
            for parent in parents:
                 if self.isDescendant(ancestor, parent):
                    return True
        return False

    def isAdjacent(self, node1, node2):
        adj = self.getAdjacencies(node1)
        for a in adj:
            if a[0] == node2 or a[1] == node2:
                return True
        return False

    def combinations(self, inSet):
        c = []
        
        for i in range(len(inSet)):
            u = inSet[i]
            for v in inSet[i+1:]:
                c.append((u, v))
        return c

    def makeDependency(self, u, v, w, isDep):
        if u < v:
            d = (u, v, w, isDep)
        else:
            d = (v, u, w, isDep)
        return d

    def calcNDependencies(self, order, n=0):
        def combinations(n, r):
            return math.factorial(n)/ (math.factorial(r) * math.factorial(n - r))
        if n == 0:
            n = len(g.nodes())
        nDeps = n * (n-1) / 2
        nCDeps = 0
        for o in range(order):
            r = o + 1
            if r > n - 2:
                break
            nCDeps += nDeps * combinations(n-2, r)
        return (nDeps, nCDeps, nDeps + nCDeps)

    def getCombinations(self, nodes=None, order=3, minOrder = 1):
        print ('order = ', order)
        from itertools import combinations
        if nodes is None:
            nodes = self.g.nodes()
        allCombos = []
        for o in range(minOrder, order+1):
            combos = combinations(nodes, o)
            allCombos += combos
        return allCombos


    def computeDependencies(self, order):
        deps = []
        nodes = list(self.g.nodes())
        nodes.sort()
        #print('nodes = ', nodes, ', order = ', order)
        cNodes = self.getCombinations(nodes, order)
        for i in range(len(nodes)):
            node1 = nodes[i]
            if not self.rvDict[node1].isObserved:
                continue
            for j in range(i, len(nodes)):
                node2 = nodes[j]
                if node1 == node2 or not self.rvDict[node2].isObserved:
                    continue
                isAdjacent = self.isAdjacent(node1, node2)
                isSeparated = not isAdjacent and networkx.d_separated(self.g, {node1}, {node2}, {})
                dep = self.makeDependency(node1, node2, None, not isSeparated)
                deps.append(dep)
                for c in cNodes:
                    if node1 in c or node2 in c:
                        continue
                    # Verify that every member of c is observed.  If not, we skip this combo.
                    allObserved = True
                    for m in c:
                        if not self.rvDict[m].isObserved:
                            allObserved = False
                            break
                    if not allObserved:
                        continue
                    isSeparated = not isAdjacent and networkx.d_separated(self.g, {node1}, {node2}, set(c))
                    dep = self.makeDependency(node1, node2, c, not isSeparated)
                    deps.append(dep)
        return deps
    
    def formatDependency(self, dep):
        # dep is:  from, to, given, isDependent
        u, v, w, isDep = dep
        if isDep:
            rel = 'is not independent from'
        else:
            rel = 'is independent from'
        if w is None:
            given = ''
        else:
            given = 'given ' + str(w)
        out = u + ' ' + rel + ' ' + v + ' ' + given
        return out

    def printDependencies(self, deps):
        print('Implied Dependencies:\n')
        for d in deps:
            print(self.formatDependency(d))

    # Test the model for consistency with a set of data.
    # Format for data is {variableName: [variable value]}.
    # That is, a dictionary keyed by variable name, containing
    # a list of data values for that variable's series.
    # The lengths of all variable's lists should match.
    # That is, the number of samples for each variable must
    # be the same.
    #
    # Returns (confidence, numTotalTests, [numTestsPerType], [numErrsPerType], [errorDetails])
    # Where:
    #   - confidence is an estimate of the likelihood that the data generating process defined
    #       by the model produced the data being tested.  Ranges from 0.0 to 1.0.
    #   - numTotalTests is the number of independencies and dependencies implied by the model.
    #   - numTestsPerType is a list, for each error type, 0 - nTypes, of the number of tests that
    #       test for the given error type.
    #   - numErrsPerType is a list, for each error type, of the number of failed tests.
    #   - errorDetails is a list of failed tests, each with the following format:
    #       [(errType, x, y, z, isDep, errStr)]
    #       Where:
    #           errType = 0 (Exogenous variables not independent) or;
    #                    1 (Expected independence not observed) or; 
    #                   2 (Expected dependence not observed)
    #           x, y, z are each a list of variable names that
    #               comprise the statement x _||_ y | z.
    #               That is x is independent of y given z.
    #           isDep True if a dependence is expected.  False for 
    #               independence
    #           pval -- The p-val returned from the independence test
    #           errStr A human readable error string describing the error
    #
    def TestModel(self, data=None, order=3):
        numTestTypes = 3
        errors = []
        if data is None:
            data = self.data
        numTestsPerType = [0] * numTestTypes
        numErrsPerType = [0] * numTestTypes
        deps = self.computeDependencies(order)
        if VERBOSE:
            print('Testing Model for', len(deps), 'Independencies')
        for dep in deps:
            x, y, zlist, isDep = dep
            X = data[x]
            Y = data[y]
            #print('X = ', X)
            #print('Y = ', Y)
            Z = []
            if zlist is None:
                zlist = []
            if zlist:
                for z in zlist:
                    zdat = data[z]
                    Z.append(zdat)
                pval = independence.test([X], [Y], Z)
            else:
                pval = independence.test([X], [Y])
            print(x, y, zlist)
            errStr = None
            testType = -1
            if not Z and self.isExogenous(x) and self.isExogenous(y):
                testType = 0
            elif not isDep:
                testType = 1
            else:
                testType = 2
            numTestsPerType[testType] += 1
            if testType == 0 and pval < .1:
                errStr = 'Error (Type 0 -- Exogenous variables not independent) -- Expected: ' + self.formatDependency(dep) + ' but dependence was detected. P-val = ' + str(pval)
            elif testType == 2 and pval > .1:
                errStr = 'Warning (Type 2 -- Unexpected independence) -- Expected: ' +  self.formatDependency(dep) + ' but no dependence detected.  P-val = ' + str(pval)
                errType = 2
            elif testType == 1 and pval < .1:
                errStr = 'Error (Type 1 -- Unexpected dependence) -- Expected: ' + self.formatDependency(dep) + ' but dependence was detected. P-val = ' + str(pval)
                errType = 1
            if errStr:
                if VERBOSE:
                    print('***', errStr)
                    #5/0
                errors.append((testType, [x], [y], list(zlist), isDep, pval, errStr))
                numErrsPerType[testType] += 1
            elif VERBOSE:
                print('.',)
        confidence = 1.0
        failurePenaltyPerType = [1, 1, 1]
        errorRatios = [0.0] * numTestTypes
        for i in range(numTestTypes):
            nTests = numTestsPerType[i]
            nErrs = numErrsPerType[i]
            if nTests > 0:
                ratio = nErrs / nTests
                errorRatios[i] = ratio
                confidence -= ratio * failurePenaltyPerType[i] / numTestTypes
        confidence = max([confidence, 0.0])
        numTotalTests = len(deps)
        if VERBOSE:
            print('Model Testing Completed with', len(errors), 'error(s).  Confidence = ', round(confidence * 100, 1), '%')
        return (confidence, numTotalTests, numTestsPerType, numErrsPerType, errors)

    def testDirection(self):
        epsilon = .0001
        rvList = list(self.rvList)
        rvList.sort()
        edges = self.g.edges
        if VERBOSE:
            print('Testing', len(edges), 'directionalities:')
        errors = 0
        for edge in edges:
            x, y = edge
            R = direction.direction(self.data[x], self.data[y])
            if abs(R) > epsilon:
                determinate = True
                if R > 0:
                    cause = x
                    effect = y
                else:
                    cause = y
                    effect = x
            else:
                determinate = False
            shouldBe = x + ' -> ' + y
            marker = ''
            if VERBOSE:
                if determinate:
                    if cause + ' -> ' + effect != shouldBe:
                        errors += 1
                        marker = '***'
                    print(marker, cause, '->', effect, ', should be: ', shouldBe, ',R = ', R)
                else:
                    if shouldBe != 'Indeterminate':
                        marker = '***'
                        errors += 1
                    print(marker, x, '--', y, ', Indeterminate, should be: ', shouldBe, ',R = ', R)
        if VERBOSE:
            print('\nerrors = ', errors, '/', len(edges), '\n')


    def findExogenous(self):
        rvList = list(self.rvList)
        rvList.sort()
        accum = {}
        for v in rvList:
            accum[v] = 0.0
        numVars = len(rvList)
        for i in range(numVars):
            x = rvList[i]
            for j in range(i+1, numVars):
                y = rvList[j]
                if x == y:
                    continue
                R = direction.direction(self.data[x], self.data[y])

                if R > 0:
                    leastCausal = y
                else:
                    leastCausal = x
                accum[leastCausal] += abs(R)
        scores = [(accum[key], key) for key in accum.keys()]
        scores.sort()
        exos = []
        for tup in scores:
            var = tup[1]
            if not exos:
                exos.append(var)
            else:
                isExo = True
                for exo in exos:
                    pval = independence.test([self.data[var]], [self.data[exo]])
                    print('ind ', var, '-', exo, ' = ', pval)
                    if pval < .1:
                        isExo = False
                        break
                        
                if isExo:
                    exos.append(var)
                else:
                    break
                        
        if VERBOSE:
            print('\nExogenous = ', exos, '\n')
            #print('scores = ', scores, accum)

    def intervene(self, targetRV, doList, controlFor = []):
        """ Implements Intverventions (Level2 of Ladder of Causality)
            of the form P(Y | do(X1=x1)).  That is, the Probability
            of Y given that we set X1 to x1.  This is generalized
            to allow multiple interventions on different variables.
            doList is the set of interventions: [(varName1, val1), ..., (varNamek, valk)].
            We return a probability distribution that can be further queried,
            e.g., as to the probability of a value, or the expected value
            (see Probability/Prob.py and pdf.py)
        """
        # Filter out any interventions for which the target is not a descendant of the
        # intevening variable.  The effect of those interventions will alsways be zero.
        doListF = []
        for item in doList:
            rv, value = item
            if targetRV in networkx.descendants(self.g, rv):
                # It is a descendant.  Keep it.
                doListF.append(item)
        if not doListF:
            # No causal effects.  Return P(target)
            return self.prob.distr(targetRV)

        # Find all the backdoor paths and identify the minimum set of variables (Z) that
        # block all such paths without opening any new paths.
        blockingSet = self.findBackdoorBlockingSet(doListF[0][0], targetRV)
        # Now we compute the probability distribution of Y conditionalized on all of the blocking
        # variables.
        given = doList + blockingSet + controlFor
        distr = self.prob.distr(targetRV, given)
        # We return the probability distribution
        return distr

    def ACE(self, cause, effect):
        """ Average Causal Effect of cause on effect.
        """
        causeDistr = self.prob.distr(cause)
        causeMean = causeDistr.E()
        causeStd = causeDistr.stDev()
        effectAtMean = self.intervene(effect, [(cause, causeMean)]).E()
        tests = [-1, -.5, -.2, .2, .5, 1 ]
        testResults = []
        for test in tests:
            testBound = causeMean + (causeStd * test)
            diff = testBound - causeMean
            effectAtBound = self.intervene(effect, [(cause, testBound)]).E()
            ace = (effectAtBound - effectAtMean) / diff
            testResults.append(ace)
        #print('testResults = ', testResults)
        tr = np.array(testResults)
        final = float(np.mean(tr))
        #print('tr = ', list(tr))
        #print('ACE = ', effectAtMean, effectAtUpper, ace)
        return final


    def CDE(self, cause, effect):
        """ Controlled Direct Effect of cause on effect
        """
        causeDistr = self.prob.distr(cause)
        causeMean = causeDistr.E()
        causeStd = causeDistr.stDev()
        if not self.isChild(cause, effect):
            # Can't have a direct effect if cause is not a parent of effect
            return 0.0
        bdBlocking = self.findBackdoorBlockingSet(cause, effect)
        fdBlocking = self.findFrontdoorBlockingSet(cause, effect)
        given = bdBlocking + fdBlocking
        effectAtMean = distr = self.prob.distr(effect, [(cause, causeMean)] + given).E()
        tests = [-1, -.5, -.2, .2, .5, 1 ]
        testResults = []
        for test in tests:
            testBound = causeMean + (causeStd * test)
            diff = testBound - causeMean
            effectAtBound = distr = self.prob.distr(effect, [(cause, testBound)] + given).E()
            cde = (effectAtBound - effectAtMean) / diff
            testResults.append(cde)
        #print('testResults = ', testResults)
        tr = np.array(testResults)
        final = float(np.mean(tr))
        return final

    def findBackdoorBlockingSet(self, source, target):
        """ Find the minimal set of nodes that block all backdoor paths from source
            to target.
        """
        maxBlocking = 3
        bSet = []
        # find all paths from parents of source to target.
        parents = self.getParents(source)            
        #print('parents = ', parents)
        # Create a graph view that removes the links from the source to its parents
        def includeEdge(s, d):
            #print('source, dest = ', s, d)
            if d == source:
                return False
            return True

        pathNodes = {}
        vg = networkx.subgraph_view(self.g, filter_edge=includeEdge)
        for parent in parents:
            paths = networkx.all_simple_paths(vg, parent, target)
            #print('paths = ', [path for path in paths])
            for path in paths:
                #print('path = ', path)
                # Remove the last node of the path -- always the target
                intermediates = path[:-1]
                #print('int = ', intermediates)
                for i in intermediates:
                    if i not in pathNodes:
                        pathNodes[i] = 1
                    else:
                        pathNodes[i] += 1
        pathTups = []
        # First look for single node solutions
        for node in pathNodes.keys():
            cnt = pathNodes[node]
            outTup = (cnt, node)
            pathTups.append(outTup)
               
        # Sort the nodes in descending order of the number of paths containing it
        pathTups.sort()
        pathTups.reverse()
        combos = [(tup[1],) for tup in pathTups]
        #print('pathNodes = ', pathNodes.keys())
        # Now add any multiple field combinations.  Order is not significant here.
        multiCombos = self.getCombinations(pathNodes.keys(), maxBlocking, minOrder=2)
        combos += multiCombos
        #print('combos = ', combos)
        for nodeSet in combos:
            testSet = set(nodeSet)
            #print('testSet = ', list(testSet))
            if networkx.d_separated(self.g, set(parents), {target}, testSet):
                bSet = list(testSet)
                break

        print('BDblocking = ', bSet)
        return bSet

    def findFrontdoorBlockingSet(self, source, target):
        backdoorSet = self.findBackdoorBlockingSet(source, target)
        maxBlocking = 2
        bSet = []
        # Create a graph view that removes the direct link from the source to the destination
        def includeEdge(s, d):
            #print('source, dest = ', s, d)
            if s == source and d == target:
                return False
            return True

        pathNodes = {}
        vg = networkx.subgraph_view(self.g, filter_edge=includeEdge)
        # Use that view to find all indirect paths from source to dest
        paths0 = networkx.all_simple_paths(vg, source, target)
        paths = [path for path in paths0]
        print('paths = ', paths)
        if len(paths) == 0:
            # No indirect paths
            return []
        for path in paths:
            #print('path = ', path)
            # Remove the first and last node of the path -- always the source and target
            intermediates = path[1:-1]
            #print('int = ', intermediates)
            for i in intermediates:
                if i not in pathNodes:
                    pathNodes[i] = 1
                else:
                    pathNodes[i] += 1
        pathTups = []
        # First look for single node solutions
        for node in pathNodes.keys():
            cnt = pathNodes[node]
            outTup = (cnt, node)
            pathTups.append(outTup)
               
        # Sort the nodes in descending order of the number of paths containing it
        pathTups.sort()
        pathTups.reverse()
        combos = [(tup[1],) for tup in pathTups]
        #print('pathNodes = ', pathNodes.keys())
        # Now add any multiple field combinations.  Order is not significant here.
        multiCombos = self.getCombinations(pathNodes.keys(), maxBlocking, minOrder=2)
        combos += multiCombos
        print('combos = ', combos)
        for nodeSet in combos:
            testSet = set(list(nodeSet) + list(backdoorSet))
            #print('testSet = ', list(testSet))
            if networkx.d_separated(vg, {source}, {target}, testSet):
                bSet = list(nodeSet)
                break

        print('FDblocking = ', bSet)
        return bSet

