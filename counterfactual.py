import math
from mpmath import mp
import pandas as pd
import numpy as np
from sklearn import metrics
from numpy.linalg import norm
import re

from intervention import Intervention


class Counterfactual:
    # DEFINITIONS
    # Identifiable: A causal quantity E[Y(X)] is identifiable if we can compute it from a purely statistical quantity
    # E[Y|X], if there is no backdoor path from X->Y.

    # rv_list: A list of random variable (rv) objects; is a SEM if fully specified.
    # data: The evidence/"reality"; a dictionary keyed by random variable name mapped to observed values.
    def __init__(self, rv_list, data, model=None, sem=None):
        self.rv_list = rv_list
        self.data = data
        self.int = Intervention(rv_list, data)
        self.cGraph = self.int.cGraph
        self.prob = self.cGraph.prob
        self.model = model
        self.sem = sem  # The SEM is a list of equations, structured the same as in the model files

    def cf(self, X, Y, cf, pdf=False, conditional=None):
        """Generates a distribution prediction of the target RV (Y) based on the counterfactual query and conditional
            clause, using the given unit observation (X) to find the "closest worlds" which will restrict the dataset

        :param X: an observation that we wish to perform counterfactual queries on; a list of tuples
                        NOTE: All exogenous variables must be specified.
        :param Y: an RV target to output (the consequent); a String (RV name)
        :param cf: the counterfactual to apply; a tuple
        :param conditional: conditional clause -- the dataset features to be used in finding the closest worlds; a list
        :return:
        """
        X = X.copy()
        x_dict = dict(X)
        y_obs = (Y, x_dict[Y])
        if Y in x_dict.keys():
            X.remove(y_obs)

        # vars will contain all the RVs that are not independent from Y, sorted in descending order of dependence with Y
        deps = [(self.prob.dependence(rv[0], Y, power=3), rv[0]) for rv in X]
        deps.sort(reverse=True)
        vars = []
        for i in range(len(deps)):
            dep = deps[i]
            if dep[0] < .5:
                print('Rejecting variables due to independence from target(p-value, var): ', deps[i:])
                break
            else:
                vars.append(dep[1])

        # Compute the SubSpace with filters based on given X with counterfactual clause replacement and Y removed
        filts = vars.copy()
        for index, rv in enumerate(filts):
            if rv == cf[0]:
                filts[index] = (rv, cf[1])
            else:
                filts[index] = (rv, x_dict[rv])
        fs = self.prob.SubSpace(filts, minPoints=10, maxPoints=0.001 * self.prob.N)
        # TODO include a check here for cf clause of each observation (delete obs if not cf)

        df = pd.DataFrame.from_dict(fs.ds)
        df_sim = self.find_similarity_scores(x_dict, Y, cf, df, list(reversed(vars)))

        out = '\nE_weighted(' + str(Y) + ') with cf ' + str(cf[0]) + '=' + str(cf[1]) + ': ' + str(df_sim[Y + "_wsim"].sum())
        if pdf:
            dist = fs.distr(Y)
            out += '\nPDF Results (Unweighted):'
            out += '\nE = ' + str(dist.E())
            out += '\nSt Dev = ' + str(dist.stDev())
            out += '\nSkew = ' + str(dist.skew())
            out += '\nKurtosis = ' + str(dist.kurtosis())

        return out

    def find_similarity_scores(self, x_dict, Y, cf, df, vars_sorted):
        score_cols = [col for col in df.columns if (col != cf[0] and col != Y)]
        cols_disc, cols_cont = [], []
        [cols_disc.append(col) if self.prob.isDiscrete(col) else cols_cont.append(col) for col in score_cols]
        df_disc, df_cont = df[cols_disc].copy(), df[cols_cont].copy()
        x_disc = np.array([v for k, v in x_dict.items() if k in cols_disc]).flatten()
        x_cont = np.array([v for k, v in x_dict.items() if k in cols_cont]).flatten()

        # Use the Jaccard coefficient for discrete variables
        def jac_coeff(row_disc):
            row = row_disc.to_numpy().reshape(-1, 1).flatten()
            jac = metrics.jaccard_score(row, x_disc, zero_division=1.0, average='weighted')
            # print('JAC / row:', row, 'vs unit:', x_disc, ':', jac)
            return jac

        # Use a Euclidean distance transformation for continuous variables
        def euc_dist(row_cont):
            row = row_cont.to_numpy().flatten()
            euc = np.linalg.norm(row - x_cont)
            sim = 1 / mp.exp(euc)  # src: https://tinyurl.com/5dfmyh3e
            # print('EUC / row:', row, 'vs unit:', x_cont, ':', round(sim, 5))
            return sim

        jac_coeffs = df_disc.apply(jac_coeff, axis=1)
        df_disc_jac = pd.DataFrame(data={"match": jac_coeffs, "idx": jac_coeffs.index})

        euc_dists = df_cont.apply(euc_dist, axis=1)
        df_cont_dist = pd.DataFrame(data={"dist": euc_dists, "idx": euc_dists.index})

        # We want the minimum non-zero distance to be ~= 0.01; 0.01 = min_dist**power
        df_cont_dist_gt_0 = df_cont_dist[abs(df_cont_dist['dist'] - 0.0) > 0.1 ** 100]
        power = mp.log(0.01) / mp.log(df_cont_dist_gt_0['dist'].min())
        df_cont_dist["dist"] = df_cont_dist["dist"] ** power

        # Compute a weighted-average similarity score from the discrete and continuous scores
        df_combined = pd.merge(df_cont_dist, df_disc_jac, on='idx')
        idx = list(df_combined['idx'])
        averaged_sim = df_combined.apply(lambda row: round(row['dist'] * (len(cols_cont) / len(df.columns)) +
                                                           row['match'] * (len(cols_disc) / len(df.columns)), 5),
                                         axis=1)
        df_combined = pd.DataFrame(data={"sim": averaged_sim, "idx": idx})
        df_combined.sort_values(by="sim", ascending=False, inplace=True)

        # Create a new dataframe with the SubSpace observations sorted by similarity score
        df_with_sim = pd.DataFrame()
        for i in range(len(df_combined)):
            idx_original = int(df_combined.iloc[i]["idx"])
            sim = df_combined.iloc[i]["sim"]
            df_with_sim = df_with_sim.append(df.loc[idx_original])
            df_with_sim.loc[idx_original, "sim"] = sim
        df_with_sim = df_with_sim[vars_sorted + [Y, 'sim']]

        # Weight the target variable predictions based on the similarity score to produce a weighted-avg expectation
        sim_sum = df_with_sim['sim'].sum()
        for i, row in df_with_sim.iterrows():
            df_with_sim.loc[i, Y + "_wsim"] = row[Y] * row['sim'] / sim_sum

        return df_with_sim

    def deterministic(self, unit, cf, target_rv):
        """Computes a deterministic query at the unit-level.
        A SEM is required, as it represents the mechanism by which each variable obtains its value.
        Pertains to a single unit of the population in which we know the value of every relevant variables.

        :param cf: the counterfactual(s) to be computed; a tuple
            Ex. "What would Y be if X were 1?" -> cf = ("X", 1)
        :param target_rv: an RV target to output (the consequent); a String (RV name)
            Ex. "What "What would Y be if X were 1?" -> target_rv = "Y"
            TODO If empty, all RVs in the model are computed using the numpy and the SEM's system of linear equations.
        :return:
        """
        self.make_sem()  # Need a fully specified SEM
        u_sem = self.make_u_sem()

        # TODO Check invertibility -- have a method in the SEM class for this ?

        # 1. Use evidence (self.data) to determine the value of U_v for each RV v in the SEM.

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

    # Returns a SEM modified to include U-vars (unit)
    # X = U_X
    # H = aX + U_H
    # Y = bX + cH + U_Y
    def make_u_sem(self):
        u_sem = {}
        for rv in self.rv_list:
            eq = rv.forwardFunc
            if len(rv.parentNames) == 0:  # If exogenous, U_X = X
                u_sem[rv.name] = rv.name
            else:  # If endogenous, U_X = X - F(X)
                rhs = eq[(re.search(' = ', eq)).span()[1]:]
                u_sem[rv.name] = rv.name + ' - ' + rhs
        return u_sem

    def solve_u_sem_with_cf(self, u_sem, unit, cf, target_rv):
        u_sem_solved = {}
        # Replace each RV in each equation with its observed value
        for rv, eq in u_sem.items():  # For each RV equation in the U-SEM
            def pad(s):
                return ' ' + s + ' '  # Pad to ensure that substrings containing this RV name are skipped

            for rv_term in u_sem.keys():  # For each RV term in this equation
                eq = eq.replace(pad(rv_term), pad(unit[rv_term]))  # Replace the RV term with its observed value
            u_sem_solved[rv] = exec(eq)

        # Now that we have all the U-terms, we can solve for the target RV using the given counterfactual
        eq_target = self.sem[target_rv] + ' + ' + u_sem_solved[target_rv]  # X = F(X) + U_X
        for rv_term in u_sem.keys():
            if rv_term == cf[0]:
                eq_target = eq_target.replace(pad(rv_term), pad(cf[1]))
            else:
                eq_target = eq_target.replace(pad(rv_term), pad(unit[rv_term]))

        return exec(eq_target)

    def check_u(self, u_space=None):
        if u_space is None:
            raise Exception("A U-distribution is needed for probabilistic CF queries.")

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
