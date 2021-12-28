
import pyemma

import numpy as np
import pandas as pd
try:
    from sklearn.cross_validation import KFold
except:
    from sklearn.model_selection import KFold
from sklearn.pipeline import Pipeline
from pyemma.coordinates.clustering import kmeans
from pyemma.coordinates.transform import tica
import logging

def _dict_compose(dict1, dict2):
    """
    Example
    -------
    >>> dict1 = {'a': 0, 'b': 1, 'c': 2}
    >>> dict2 = {0: 'A', 1: 'B'}
    >>> _dict_compose(dict1, dict2)
    {'a': 'A', 'b': 'b'}
    """
    return {k: dict2.get(v) for k, v in dict1.items() if v in dict2}

class MaximumLikelihoodMSM(pyemma.msm.MaximumLikelihoodMSM):
    def __init__(self, lag: int, score_k=1) -> None:
        super().__init__(lag=lag, score_k=score_k)
        #self._estimate(dtrajs)
    
    def score_GMRQ(self, dtrajs: list) -> float:
        # eigenvectors from the model we're scoring, `self`
        V = self.eigenvectors_right()

        # Note: How do we deal with regularization parameters like prior_counts
        # here? I'm not sure. Should C and S be estimated using self's
        # regularization parameters?
        m2 = self.__class__(**self.get_params())
        m2.fit(dtrajs)

        if self.mapping_ != m2.mapping_:
            V = self._map_eigenvectors(V, m2.mapping_)
            # we need to map this model's eigenvectors
            # into the m2 space

        # How well do they diagonalize S and C, which are
        # computed from the new test data?
        S = np.diag(m2.pi)
        C = S.dot(m2.transition_matrix)

        try:
            gmrq = np.trace(V.T.dot(C.dot(V)).dot(np.linalg.inv(V.T.dot(S.dot(V)))))
        except np.linalg.LinAlgError:
            gmrq = np.nan

        return gmrq
    
    def _map_eigenvectors(self, V: np.array, other_mapping: dict) -> np.array:
        """Mapping correspond

        Args:
            V (np.array): right eigvectors
            other_mapping (dict): mapping_ of new model

        Returns:
            np.array: [description]
        """
        self_inverse_mapping = {v: k for k, v in self.mapping_.items()}
        transform_mapping = _dict_compose(self_inverse_mapping, other_mapping)
        source_indices, dest_indices = zip(*transform_mapping.items())

        #print(source_indices, dest_indices)
        mapped_V = np.zeros((len(other_mapping), V.shape[1]))
        mapped_V[dest_indices, :] = np.take(V, source_indices, axis=0)
        return mapped_V

    @property
    def mapping_(self) -> dict:
        """[summary]

        Returns:
            dict: [description]
        """
        mapping = {}
        for i, state in enumerate(self.active_set):
            mapping[state] = i
        return mapping

class TICA(tica.TICA):

    @property
    def mean(self):
        n_features = self.data[0].shape[0]
        n_observations_ = 0
        n_sequences_ = len(self.data)
        _sum_0_to_TminusTau = np.zeros(n_features)
        _sum_tau_to_T = np.zeros(n_features)
        for X in self.data:
            _sum_0_to_TminusTau += X[:-self.lag].sum(axis=0)
            n_observations_ += X.shape[0]
            _sum_tau_to_T += X[self.lag:].sum(axis=0)
        two_N = 2 * (n_observations_ - self.lag * n_sequences_)
        means = (_sum_0_to_TminusTau + _sum_tau_to_T) / float(two_N)
        return means

class MSMPipeline(Pipeline):
    def fit(self, trainData: list) -> None:
        output = trainData
        for name, obj in self.steps[:-1]:
            obj.fit(output)
            if name == 'cluster':
                output = obj.dtrajs
            else:
                output = obj.transform(output)
        _, obj = self.steps[-1]
        obj.fit(output)
    
    def score(self, testData: list) -> float:
        output = testData
        for name, obj in self.steps[:-1]:
            if name == 'cluster':
                output = [ np.concatenate(d) for d in obj.transform(output)]
            else:
                output = obj.transform(output)
        _, obj = self.steps[-1]
        return obj.score(output)


class GMRQValid(object):

    def __init__(self) -> None:
        super().__init__()

    def tICA_para_nICs(self, data, nICs=None, lag=1, kinetic_map=True, commute_map=False, n_folds=5):
        cv = KFold(n_splits=n_folds).split(data)
        dim = data[0].shape[1]
        if nICs is None:
            nICs = range(5, dim)
        results = []
        #tica = pyemma.coordinates.tica(lag=lag, kinetic_map=kinetic_map, commute_map=commute_map)
        model = MSMPipeline([
                         ('tica', tica.TICA(lag=lag, kinetic_map=kinetic_map, commute_map=commute_map)),
                         ('cluster', kmeans.KmeansClustering(fixed_seed=43, n_clusters=100)),
                         ('msm', MaximumLikelihoodMSM(lag=1))
                        ])
        for nIC in nICs:
            print(nIC)
            print(model)
            model.set_params(tica__dim=nIC)
            for foldIndex, (trainIndex, testIndex) in enumerate(cv):
                trainData = [data[i] for i in trainIndex]
                testData = [data[i] for i in testIndex]

                trainScore = 0
                testScore = 0
                #try:
                model.fit(trainData)
                trainScore = model.score(trainData)
                testScore = model.score(testData)
                #except:
                #    logging.warn('Failed train the data for fold %d at nIC %d'%(foldIndex, nIC))
                results.append({
                    'train_score': trainScore,
                    'test_score': testScore,
                    'nIC': nIC,
                    'fold': foldIndex})
        results = pd.DataFrame(results)
        avgs = (results
            .groupby('nIC')
            .aggregate(np.median)
            .drop('fold', axis=1))
        print(results)
        print(avgs)
        #bestN = avgs['test_score'].argmax()
        #bestScore = avgs.loc[bestN, 'test_score']
        #print(bestN, "gives the best score:", bestScore)
        return results