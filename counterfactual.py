from intervention import Intervention
import math
import pandas as pd
import numpy as np
from sklearn import metrics
from sklearn import preprocessing
from scipy import spatial
from numpy import dot
from numpy.linalg import norm


class Counterfactual:
    # DEFINITIONS
    # Identifiable: A causal quantity E[Y(X)] is identifiable if we can compute it from a purely statistical quantity
    # E[Y|X], if there is no backdoor path from X->Y.

    # rv_list: A list of random variable (rv) objects; is a SEM if fully specified.
    # data: The evidence/"reality"; a dictionary keyed by random variable name mapped to observed values.
    def __init__(self, rv_list, data):
        self.rv_list = rv_list
        self.data = data
        self.int = Intervention(rv_list, data)
        self.cGraph = self.int.cGraph
        self.prob = self.cGraph.prob
        self.is_sem = self.set_sem()

    def cf_closest_worlds(self, unit, cf, target_rv, conditional=None, power=1):
        """Finds the "closest world" to the given unit in which the given counterfactual(s) hold.

        :param unit: an observation that we wish to perform counterfactual queries on; a dictionary
        :param cf: the counterfactual to apply; a tuple
        :param conditional:
        :param target_rv: an RV target to output; the consequent
        :param power:
        :return:
        """
        if unit.keys() != self.data.keys():
            raise Exception("Invalid unit observation entered.")

        if conditional is None:
            # then we condition on all variables except cf, target_rv
            conditional = [rv_name for rv_name in unit.keys()
                           if rv_name != cf[0] and rv_name != target_rv]

        # TODO sort conditional clause -- ascending order of independence from cf?

        df_worlds = self.find_closest_worlds(unit, cf, conditional)
        data_worlds = df_worlds.to_dict('list')  # convert the closest_worlds dataframe to a data dictionary

        print(data_worlds)

    def find_closest_worlds(self, unit, cf, conditional_rv):
        # Restrict the dataset based on the conditional clause -- find the k closest worlds (1%)
        num_obs = len(self.data[list(self.data)[0]])
        k = int(round(.01 * num_obs))

        df = pd.DataFrame.from_dict(self.data)

        df_obs = df[df[cf[0]] == cf[1]]  # Conditionalize on the counterfactual
        df_obs = df_obs[conditional_rv].copy()  # Filter columns on the given cond clauses
        k = min(k, len(df_obs))

        # Filter the unit observation to also include only conditional columns
        unit_cond = {}
        for rv_name, rv_val in unit.items():
            if rv_name in conditional_rv:
                unit_cond[rv_name] = rv_val

        # Split the dataset into continuous and discrete vars
        cols_disc = []
        cols_cont = []
        [cols_disc.append(col) if self.prob.isDiscrete(col) else cols_cont.append(col) for col in df_obs.columns]

        df_disc = df_obs[cols_disc].copy()
        unit_disc = {k: v for k, v in unit_cond.items() if k in cols_disc}
        unit_disc = np.array(list(unit_disc.values())).flatten()

        df_cont = df_obs[cols_cont].copy()
        unit_cont = {k: v for k, v in unit_cond.items() if k in cols_cont}
        unit_cont = np.array(list(unit_cont.values())).flatten()

        def jac_coeff(row_disc):
            row = row_disc.to_numpy().reshape(-1, 1).flatten()
            jac = metrics.jaccard_score(row, unit_disc, zero_division=1.0, average='weighted')
            print('JAC / row:', row, 'vs unit:', unit_disc, ':', jac)
            return jac

        def euc_dis(row_cont):
            row = row_cont.to_numpy().flatten()
            euc = np.linalg.norm(row - unit_cont)
            sim = 1 / math.exp(euc)  # src: https://tinyurl.com/5dfmyh3e
            print('EUC / row:', row, 'vs unit:', unit_cont, ':', sim)
            return sim
            # cos = 1 - spatial.distance.cosine(row, unit_cont)
            # print('COS / row:', row, 'vs unit:', unit_cont, ':', cos)
            # return cos

        euc_dists = df_cont.apply(euc_dis, axis=1)
        df_cont_dist = pd.DataFrame(data={"sim": euc_dists, "idx": euc_dists.index})

        jac_coeffs = df_disc.apply(jac_coeff, axis=1)
        df_disc_jac = pd.DataFrame(data={"match": jac_coeffs, "idx": jac_coeffs.index})

        df_combined = pd.merge(df_cont_dist, df_disc_jac, on='idx')
        averaged_sim = df_combined.apply(lambda row: round(row['sim'] * (len(cols_cont) / len(df_obs.columns)) +
                                                           row['match'] * (len(cols_disc) / len(df_obs.columns)), 2),
                                         axis=1)
        df_combined = pd.DataFrame(data={"sim": averaged_sim, "idx": averaged_sim.index})
        df_combined.sort_values(by="sim", ascending=False, inplace=True)

        print(df_combined)

        k_closest_obs = pd.DataFrame()
        for i in range(0, k):
            idx_original = df_combined.iloc[i]["idx"]
            k_closest_obs = k_closest_obs.append(df.loc[int(idx_original)])

        return k_closest_obs

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
    # u_space: Probability space for exogenous "Universe" variables U. Provides nondeterminism; "represent[s] our
    # uncertainty as to the identity of the subject under consideration or, when the subject is known, what other
    # characteristics that subject has that might have bearing on our problem."
    def probabilistic(self, cf, u_space=None, rv_list=[]):
        self.check_sem()
        self.check_u(u_space)
        pass

    def act(self, rv, val):
        self.data[rv] = val

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
        if u_space is None:
            raise Exception("A U-distribution is needed for probabilistic CF queries.")