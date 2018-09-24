import numpy as np
from scipy.stats import spearmanr,pearsonr,kendalltau

def make_pairwise_correlation_matrix(counts_matrix,correlation='pearson',
                                     lower_triangle_only=True):
    """
    Computes a matrix of pairwise correlation coefficients.

    Parameters
    ----------
    counts_matrix : np.ndarray
        A 2 x 2 matrix where the columns are different libraries, and the rows
        are gene counts across libraries. 
    correlation : {'pearson', 'spearman','kendalltau'}
        Controls which correlation will be used.
    lower_triangle_only : bool
        If True, only the lower triangle of the correlation matrix will be
        filled in. 
    Returns
    -------
    np.ndarray
        The pairwise correlation matrix.
    """
    # resolve correlation
    corr_fn = None
    if correlation == 'pearson':
        corr_fn = pearsonr
    elif correlation == 'spearman':
        corr_fn = spearmanr
    elif correlation == 'kendalltau':
        corr_fn = kendalltau

    # transpose column-major counts_matrix to row major
    if counts_matrix.shape[0] > counts_matrix.shape[1]:
        counts_matrix = counts_matrix.T

    # compute matrix of pairwise correlation coefficients
    correlation_matrix = np.zeros((len(counts_matrix), len(counts_matrix)))
    for i in range(len(correlation_matrix)):
        for j in range(i + 1):
            if i == j:
                correlation_matrix[i, j] = 1.0
            else:
                corr_value = corr_fn(counts_matrix[i], counts_matrix[j])[0]
                correlation_matrix[i, j] = corr_value
                if not lower_triangle_only:
                    correlation_matrix[j, i] = corr_value

    return correlation_matrix
