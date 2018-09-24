from __future__ import division
import numpy as np

def normalize_by_tpm(gene_order,gene_info,counts_matrix):
    """
    Normalizes libraries per kilobase of total transcript length per million
    reads.

    Parameters
    ----------
    gene_order : List[str]
        The list of genes that correspond to each row of the counts_matrix
    gene_info : Dict[str,Dict[str,any]]
        A superdict containing information about each gene. The outer keys are
        the UCSC gene names, and the values contain the following fields:
        
            {
            'name' : str
            'accession' : str
            'chr' : str
            'start' : int
            'end' : int
            'length : int
            }

    counts_matrix : np.ndarray
        The counts matrix of all sequencing libraries. The rows are the counts
        for each gene and the columns represent the different libraries.

    Returns
    -------
    np.ndarray
        A 2 x 2 matrix like counts_matrix, but normalized for transcript length
        and sequencing depth. 
    """

    normalized_matrix = np.zeros_like(counts_matrix)

    for i in range(counts_matrix.shape[0]):
        try:
            k = gene_info[gene_order[i]]['length'] / 1000
            normalized_matrix[i] = counts_matrix[i] / k
        except KeyError:
            print 'gene information is not available for %s'%gene_order[i]
            normalized_matrix[i,0] = 0
            normalized_matrix[i,1] = 0
            normalized_matrix[i,2] = 0
            normalized_matrix[i,3] = 0


    rpm = np.nansum(normalized_matrix,axis=0) / 1000000
    tpm_matrix = normalized_matrix / rpm
    

    return tpm_matrix 
