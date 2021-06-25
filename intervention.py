import networkx
import numpy as np

from cGraph import cGraph


class Intervention():
    def __init__(self, rvList, data):
        self.cGraph = cGraph(rvList, data)
        self.prob = self.cGraph.prob

    def intervene(self, targetRV, doList, controlFor = [], power=1):
        """ Implements Intverventions (Level 2 of Ladder of Causality)
            of the form P(Y | do(X1=x1),Z).  That is, the Probability
            of Y given that we set X1 to x1 and control for Z.  This is generalized
            to allow multiple interventions on different variables.
            doList is the set of interventions: [(varName1, val1), ..., (varNamek, valk)].
            We return a probability distribution that can be further queried,
            e.g., as to the probability of a value, or the expected value
            (see Probability/Prob.py and pdf.py)
        """
        # Filter out any interventions for which the target is not a descendant of the
        # intervening variable.  The effect of those interventions will always be zero.
        doListF = []
        for item in doList:
            rv, value = item
            if targetRV in networkx.descendants(self.g, rv):
                # It is a descendant.  Keep it.
                doListF.append(item)
        if not doListF:
            # No causal effects.  Return P(target)
            return self.prob.distr(targetRV, controlFor)

        # Find all the backdoor paths and identify the minimum set of variables (Z) that
        # block all such paths without opening any new paths.
        blockingSet = self.cGraph.findBackdoorBlockingSet(doListF[0][0], targetRV)
        # Now we compute the probability distribution of Y conditionalized on all of the blocking
        # variables.
        given = doList + blockingSet + controlFor
        distr = self.prob.distr(targetRV, given, power=power)
        # We return the probability distribution
        return distr

    def ACE(self, cause, effect, power=1):
        """ Average Causal Effect of cause on effect.
        """
        causeDistr = self.prob.distr(cause)
        causeMean = causeDistr.E()
        causeStd = causeDistr.stDev()
        tests = [.2, .5, 1 ]
        testResults = []
        for test in tests:
            lowBound = causeMean - (causeStd * test)
            highBound = causeMean + (causeStd * test)
            diff = highBound - lowBound
            effectAtLow = self.intervene(effect, [(cause, lowBound)], power=power).E()
            effectAtHigh = self.intervene(effect, [(cause, highBound)], power=power).E()
            ace = (effectAtHigh - effectAtLow) / diff
            testResults.append(ace)
        #print('testResults = ', testResults)
        tr = np.array(testResults)
        final = float(np.mean(tr))
        #print('ACE = ', effectAtMean, effectAtUpper, ace)
        return final

    # Total effect: measures the expected increase in effect as the cause changes from 0 to 1,
    # while the mediator is allowed to track the change in the cause naturally
    def TE(self, cause, effect):
        pass

    # Controlled direct effect: measures the expected increase in effect as the cause changes from 0 to 1,
    # while the mediator is set to a specified level uniformly over the population
    # Holds the mediator at a constant level for the entire population (pg 77)
    def CDE(self, cause, effect, mediator):
        pass
