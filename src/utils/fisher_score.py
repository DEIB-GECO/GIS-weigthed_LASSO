""" Script to compute the Fisher's score.

Reference:
-----
Li, J., Cheng, K., Wang, S., Morstatter, F., Trevino, R. P., Tang, J., & Liu, H. (2017)
Feature Selection: A Data Perspective,
ACM Computer Surveys, 50(6),
https://doi.org/10.1145/3136625

Original code:
-----
https://github.com/jundongl/scikit-feature/blob/master/skfeature/function/similarity_based/fisher_score.py

"""


import numpy as np
from scipy.sparse import *
from utils. construct_W import construct_W


def fisher_score(X, y):
    """
    This function implements the fisher score feature selection, steps are as follows:
    1. Construct the affinity matrix W in fisher score way
    2. For the r-th feature, we define fr = X(:,r), D = diag(W*ones), ones = [1,...,1]', L = D - W
    3. Let fr_hat = fr - (fr'*D*ones)*ones/(ones'*D*ones)
    4. Fisher score for the r-th feature is score = (fr_hat'*D*fr_hat)/(fr_hat'*L*fr_hat)-1

    Input
    -----
    X: {numpy array}, shape (n_samples, n_features)
        input data
    y: {numpy array}, shape (n_samples,)
        input class labels

    Output
    ------
    score: {numpy array}, shape (n_features,)
        fisher score for each feature

    Reference
    ---------
    He, Xiaofei et al. "Laplacian Score for Feature Selection." NIPS 2005.
    Duda, Richard et al. "Pattern classification." John Wiley & Sons, 2012.
    """

    # Construct weight matrix W in a fisherScore way
    kwargs = {"neighbor_mode": "supervised", "fisher_score": True, 'y': y}
    W = construct_W(X, **kwargs)
    #print(W.shape)
    # build the diagonal D matrix from affinity matrix W
    D = np.array(W.sum(axis=1))
    #print(D.sum())
    L = W
    tmp = np.dot(np.transpose(D), X)
    D = diags(np.transpose(D), [0])
    Xt = np.transpose(X)
    t1 = np.transpose(np.dot(Xt, D.todense()))
    t2 = np.transpose(np.dot(Xt, L.todense()))
    #print(t1.shape)

    #print(t2.shape)
    # compute the numerator of Lr
    #print(np.sum(np.multiply(t1, X), 0).shape)
    #print(np.multiply(tmp, tmp)[0].shape)
    #print(np.multiply(tmp, tmp)[0]/D.sum())
    D_prime = np.sum(np.multiply(t1, X), 0) - np.multiply(tmp, tmp)[0]/D.sum()
    # compute the denominator of Lr
    L_prime = np.sum(np.multiply(t2, X), 0) - np.multiply(tmp, tmp)[0]/D.sum()
    #print(D_prime.shape)
    #print(L_prime.shape)
    # avoid the denominator of Lr to be 0
    D_prime[D_prime < 1e-12] = 10000
    lap_score = 1 - np.array(np.multiply(L_prime, 1/D_prime))

    # compute fisher score from laplacian score, where fisher_score = 1/lap_score - 1
    score = 1.0/lap_score - 1
    return np.transpose(score)


def feature_ranking(score):
    """
    Rank features in descending order according to fisher score, the larger the fisher score, the more important the
    feature is
    """
    idx = np.argsort(score, 0)
    return idx[::-1]