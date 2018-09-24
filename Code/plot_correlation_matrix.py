"""
Module for plotting pairwise correlation matrices.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def plot_correlation_matrix(matrix, filename, lower_triangle_only=False, label_values=None, cluster=False,
                            cbar=False):
    """
    Plots a pairwise corrrelation matrix as a heatmap.

    Parameters
    ----------
    matrix : np.ndarray
        The pairwise correlation matrix to visualize.
    filename : str
        String reference to the file write the plot to.
    lower_triangle_only : bool
        If True, only the lower triangle of the correlation matrix will be
        plotted. 
    label_values : Optional[List[str]]
        A list of strings labeling the columns of the matrix. If not passed, no
        labels will be included.
    cluster : bool
        Pass True to perform heirarchical clustering on the rows and columns of
        the matrix.
    cbar : bool
        Pass True to include a colorbar.
    """
    plt.clf()
    if len(matrix) > 6:
        plt.figure(figsize=(len(matrix), len(matrix)))
    if cluster:
        cm = sns.clustermap(data=matrix, metric='cosine', square=True,
                            annot=True, xticklabels=label_values,
                            yticklabels=label_values)
        if label_values is not None:
            plt.setp(cm.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
            plt.setp(cm.ax_heatmap.xaxis.get_majorticklabels(), rotation=45)
            plt.setp(cm.ax_heatmap.xaxis.get_majorticklabels(), ha='right')
    else:
        if lower_triangle_only:
            mask = np.zeros_like(matrix)
            mask[np.triu_indices_from(mask,k=1)] = True
            with sns.axes_style('white'):
                sns.heatmap(data=matrix, mask=mask, square=True, annot=True,
                    xticklabels=label_values,yticklabels=label_values,
                    cbar=cbar)
        else:
            sns.heatmap(data=matrix, square=True, annot=True,
                    xticklabels=label_values, yticklabels=label_values,
                    cbar=cbar)
        if label_values is not None:
            plt.yticks(rotation=0)
            plt.xticks(rotation=45, ha='right')
    plt.savefig(filename, bbox_inches='tight')
