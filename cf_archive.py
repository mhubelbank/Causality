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
    def __init__(self, rv_list, data):
        self.rv_list = rv_list
        self.data = data
        self.int = Intervention(rv_list, data)
        self.cGraph = self.int.cGraph
        self.prob = self.cGraph.prob
        self.sem = None  # The SEM is a list of equations, structured the same as in the model files

    def cf_closest_worlds(self, unit, cf, target_rv, conditional=None):
        """Generates a distribution prediction of the target RV based on the counterfactual query and conditional
            clause, using the given unit observation to find the "closest worlds" which will restrict the dataset

        :param unit: an observation that we wish to perform counterfactual queries on; a dictionary (TODO CHANGE TO LIST OF TUPLES)
                        NOTE: All exogenous variables must be specified.
        :param cf: the counterfactual to apply; a tuple
        :param conditional: conditional clause -- the dataset features to be used in finding the closest worlds; a list
        :param target_rv: an RV target to output (the consequent); a String (RV name)
        :return: a ProbSpace prediction for the target RV
        """
        # unit = self.make_unit(unit)
        # if unit.keys() != self.data.keys():
        #     raise Exception("Invalid unit observation entered; each RV in the dataset must be accounted for. "
        #                     "RV list: " + self.data.keys())
        # TODO also need an observed target RV
        # If no conditional clause is given, then we condition on all variables except cf, target_rv
        if conditional is None:
            conditional = [rv_name for rv_name in unit.keys() if rv_name != cf[0] and rv_name != target_rv]
        elif cf[0] in conditional or target_rv in conditional:
            raise Exception("The conditional clause cannot contain the counterfactual or target RV variables.")

        # TODO sort conditional clause -- ascending order of independence from cf?

        print('Given Observation:\n' + str(pd.DataFrame(unit, index=[0]).to_string(index=False)))
        df_worlds = self.find_closest_worlds(unit, cf, conditional)
        data_worlds = df_worlds.to_dict('list')  # Convert the closest_worlds dataframe to a data dictionary
        print('\nNew Dataset (n=' + str(len(df_worlds)) + '):\n' + str(df_worlds))

        return self.prob.PredictDist(target_rv, data_worlds)

    def find_closest_worlds(self, unit, cf, conditional):
        """Finds the "closest worlds" to the given unit in which the given counterfactual(s) hold. Uses the Jaccard
        coefficient for categorical RVs and Euclidean distance for continuous RVs.

        :param unit: an observation that we wish to perform counterfactual queries on; a dictionary
        :param cf: the counterfactual to apply; a tuple
        :param conditional: conditional clause -- the dataset features to be used in finding the closest worlds; a list
        :return: a dataframe with the k closest observations sorted by similarity
        """
        df = pd.DataFrame.from_dict(self.data)
        df_obs = df[df[cf[0]] == cf[1]]  # Conditionalize on the counterfactual
        df_obs = df_obs[conditional].copy()  # Filter columns on the given cond clauses
        unit_cond = {k: v for k, v in unit.items() if k in conditional}  # Filter unit obs on the given cond clauses

        # Restrict the dataset based on the conditional clause -- find the k closest worlds (1%)
        num_obs = len(self.data[list(self.data)[0]])
        k = min(int(round(.01 * num_obs)), len(df_obs))

        print('\nGiven Dataset:\n', df)
        print('\nWe filter on counterfactual', cf[0], '=', cf[1], 'then find the', k, 'closest worlds based on on RVs',
              conditional, ':')

        # Split the dataset into continuous and discrete vars
        cols_disc, cols_cont = [], []
        [cols_disc.append(col) if self.prob.isDiscrete(col) else cols_cont.append(col) for col in df_obs.columns]
        df_disc, df_cont = df_obs[cols_disc].copy(),     df_obs[cols_cont].copy()

        unit_disc = np.array([v for k, v in unit_cond.items() if k in cols_disc]).flatten()
        unit_cont = np.array([v for k, v in unit_cond.items() if k in cols_cont]).flatten()

        # Use the Jaccard coefficient for discrete variables
        def jac_coeff(row_disc):
            row = row_disc.to_numpy().reshape(-1, 1).flatten()
            jac = metrics.jaccard_score(row, unit_disc, zero_division=1.0, average='weighted')
            # print('JAC / row:', row, 'vs unit:', unit_disc, ':', jac)
            return jac

        # Use a Euclidean distance transformation for continuous variables
        def euc_dist(row_cont):
            row = row_cont.to_numpy().flatten()
            euc = np.linalg.norm(row - unit_cont)
            sim = 1 / mp.exp(euc)  # src: https://tinyurl.com/5dfmyh3e
            # print('EUC / row:', row, 'vs unit:', unit_cont, ':', round(sim, 5))
            return sim
            # cos = 1 - spatial.distance.cosine(row, unit_cont)
            # print('COS / row:', row, 'vs unit:', unit_cont, ':', cos)
            # return cos

        jac_coeffs = df_disc.apply(jac_coeff, axis=1)
        df_disc_jac = pd.DataFrame(data={"match": jac_coeffs, "idx": jac_coeffs.index})

        euc_dists = df_cont.apply(euc_dist, axis=1)
        df_cont_dist = pd.DataFrame(data={"dist": euc_dists, "idx": euc_dists.index})

        # We want the minimum non-zero distance to be ~= 0.01; 0.01 = min_dist**power
        df_cont_dist_gt_0 = df_cont_dist[abs(df_cont_dist['dist'] - 0.0) > 0.1 ** 100]
        power = mp.log(0.01) / mp.log(df_cont_dist_gt_0['dist'].min())
        df_cont_dist["dist"] = df_cont_dist["dist"] ** power
        print('\nAfter scaling of Euclidean distances using power',
              str(round(power, 5)) + ':\n' + str(df_cont_dist["dist"]))

        df_combined = pd.merge(df_cont_dist, df_disc_jac, on='idx')
        idx = list(df_combined['idx'])
        averaged_sim = df_combined.apply(lambda row: round(row['dist'] * (len(cols_cont) / len(df_obs.columns)) +
                                                           row['match'] * (len(cols_disc) / len(df_obs.columns)), 5),
                                         axis=1)
        df_combined = pd.DataFrame(data={"sim": averaged_sim, "idx": idx})
        df_combined.sort_values(by="sim", ascending=False, inplace=True)
        print('\nTop 100 similarity scores descending:', df_combined.head(100)['sim'].to_list())

        k_closest_obs = pd.DataFrame()
        for i in range(0, k):
            idx_original = df_combined.iloc[i]["idx"]
            k_closest_obs = k_closest_obs.append(df.loc[int(idx_original)])

        return k_closest_obs[unit.keys()]
