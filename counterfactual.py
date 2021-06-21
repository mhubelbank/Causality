from cGraph import cGraph


class Counterfactual:
    # DEFINITIONS
    # Identifiable: A causal quantity E[Y(X)] is identifiable if we can compute it from a purely statistical quantity
    # E[Y|X], if there is no backdoor path from X->Y.

    # rv_list: A list of random variable (rv) objects; is a SEM if fully specified.
    # data: The evidence/"reality"; a dictionary keyed by random variable name mapped to its observed value.
    # u_space: Probability space for exogenous "Universe" variables U. Provides nondeterminism; "represent[s] our
    # uncertainty as to the identity of the subject under consideration or, when the subject is known, what other
    # characteristics that subject has that might have bearing on our problem."
    def __init__(self, rv_list, data, u_space=None):
        self.rv_list = rv_list
        self.data = data
        self.cGraph = cGraph(rv_list, data)
        self.is_sem = self.set_sem()

    def find_closest_world(self, record, cf, rv_list=None):
        """Finds the "closest world" to the given record in which the given counterfactual(s) hold.

        :param record: an observation that we wish to perform counterfactual queries on
        :param cf: the counterfactual(s) to apply
        :param rv_list: (optional) a list of RV targets to output
        :return: dictionary, with keyset filtered on rv_list if included
        """

        pass
    
    def deterministic(self, cf, rv_list=[]):
        """Computes a deterministic query at the unit-level.
        A SEM is required, as it represents the mechanism by which each variable obtains its value.
        Pertains to a single unit of the population in which we know the value of every relevant variables.

        :param cf: The counterfactual(s) to be computed.
            Ex. "What would Y be if X were 1?" -> cf = {"X": 1}
        :param rv_list: The desired consequence(s) of the counterfactual.
            If empty, all RVs in the model are computed using the numpy and the SEM's system of linear equations.
            Ex. "What "What would Y be if X were 1?" -> rv_list = ["Y"]
        :return:
        """
        self.check_sem()  # Need a fully specified SEM

        # Check invertibility -- have a method in the SEM class for this ?

        # 1. Use evidence (self.data) to determine the value of U.

        # 2. Modify the model, M, by removing the structural equations for the variables in X and replacing them with
        # the appropriate functions X=x, to obtain the modified model, M_x.

        # 3. Use the modified model M_x and the value of U to compute the value of Y (consequence of counterfactual).
        # Use numpy to solve system of linear equations?

        pass

    # Useful for non-invertible unit-level counterfactual queries.
    def probabilistic(self, cf, u_space=None, rv_list=[]):
        self.check_sem()
        self.check_u(u_space)
        pass

    def act(self, rv, val):
        self.SEM[rv] = val

    # Counterfactuals in linear models
    # Compute the effect of treatment on the treated
    # Allows us to answer hypothetical questions about individuals from population data.
    # ETT = The effect of treatment on the entire population
    # Assumes a linear system (see p106-107)
    def ETT(self, cause, effect, condition_on):
        pass

    # Excess risk ratio

    # def ERR(self):
    #     pass

    # Probability of necessity -- measures the degree to which a decision was necessary for the given outcome.
    def PN(self):
        pass

    # Probability of sufficiency -- measures the degree to which the action not taken would have been sufficient for
    # the desired outcome.
    def PS(self):
        pass

    # Probability of necessity AND sufficiency
    # Assumes monotonicity (consistency) between treatment and effect
    def PNS(self):
        pass

    def mediate(self):
        # Call all the below functions and print output using ratios at top of p122/123
        # Can we identify the natural direct and indirect effects?
        # Can the total effect be decomposed into the sum of the natural direct and indirect effects?
        # What fraction of the effect is necessarily due to a mediator?
        pass

    # Natural direct effect -- measures the expected increase in effect as the treatment changes from 0 to 1,
    # while the mediator is set to whatever value it would have attained for each individual prior to the change (T=0)
    def NDE(self, cause, effect, mediator):
        pass

    # Natural indirect effect -- measures the expected increase in effect when the treatment is held constant at T=0
    # and the mediator changes to whatever value it would have attained for each individual under T=1.
    # Captures the portion of the effect that can be explained by mediation alone, while disabling the capacity
    # of the effect to respond to the treatment.
    def NIE(self, cause, effect, mediator):
        # Calculate via formula AND TE - NDE
        pass

    # Determines whether the NDE and NIE can be identified.
    # We can identify these provided that there exists a set W of measured covariates such that:
    # 1. No member of W is a descendant of T
    # 2. W blocks all backdoor paths from M to Y (after removing T->M and T->Y)
    # 3. The W-specific effect of T on M is identifiable (using experiments or adjustments)
    # 4. The W-specific joint effect of {T,M} on Y is identifiable (using experiments or adjustments)
    def can_identify_natural_effects(self):
        pass

    # Measures the extent to which the effect of X on Y is explained by X's effect on the mediator Q.
    def mediation(self):
        pass

    # Determine whether this Counterfactual instance's rv_list contains a fully specified SEM.
    def set_sem(self):
        for rv in self.rv_list:
            if rv.forwardFunc is None:
                return False
        return True

    # Need a fully specified parametric model (SEM) to compute unit-level counterfactual queries.
    # data: Need to have data for every key in SEM and vice versa.
    def check_sem(self):
        if not self.is_sem:
            raise Exception("A fully specified SEM is needed for deterministic CF queries.")
        rv_names = [rv.name for rv in self.rv_list]
        if set(self.data.keys()) != set(rv_names):
            raise Exception("A complete observation is needed for deterministic CF queries.")

    def check_u(self, u_space=None):
        if u_space is None and self.u_space is None:
            raise Exception("A U-distribution is needed for probabilistic CF queries.")