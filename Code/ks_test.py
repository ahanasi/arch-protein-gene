from scipy import stats
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as pyplot


def ks_test(file_string):


    normed = np.loadtxt(file_string, usecols=(1, 2, 3, 4), dtype = 'float', skiprows=1)
    gene_name = np.loadtxt(file_string, usecols=(0,), skiprows=1, dtype = 'string')


    pvalues_ks_test = np.zeros((np.shape(normed)[0], 2))
    ks_stats = np.zeros((np.shape(normed)[0],2))


    for i in range(0,(np.shape(normed)[0])):
        ks_stats[i,0] = np.mean(normed[i,2:4])
        ks_stats[i,1] = np.mean(normed[i,0:2])
        ks_output = stats.ks_2samp(normed[i,0:2],normed[i,2:4])
        print ks_output[0],ks_output[1]
        pvalues_ks_test[i,0] = -10*np.log2(ks_output[1])
        pvalues_ks_test[i,1] = np.log2(np.divide(((np.mean(normed[i,2:4]))+1),((np.mean(normed[i,0:2]))+1)))


    #Make MA plots to look at fold change with respect to abundance
    m = np.log2(np.divide((ks_stats[:,0]+1),((ks_stats[:,1]+1))))
    a = a = 0.5*(np.log2((ks_stats[:,0]+1) + ((ks_stats[:,1]+1))))

    pyplot.figure()
    pyplot.scatter(a,m)
    pyplot.scatter(a[np.where(pvalues_ks_test[:,0] >= -10*np.log2(0.05))],m[np.where(pvalues_ks_test[:,0] >= -10*np.log2(0.05))],color='red')
    pyplot.title('MA Plot: ES vs. NPC - pvalues from K-S Test Code')
    pyplot.show()
    


    # Make the volcano plots to look at pvalues with respect to fold change

    pyplot.figure()
    pyplot.scatter(pvalues_ks_test[:,1],pvalues_ks_test[:,0])
    pyplot.scatter(pvalues_ks_test[np.where(pvalues_ks_test[:,0]>=-10*np.log2(0.05)),1],pvalues_ks_test[np.where(pvalues_ks_test[:,0]>=-10*np.log2(0.05)),0],color='red')
    pyplot.title('Volcano plot: ES vs. NPC - pvalues from K-S Test Code')
    pyplot.show()

    return pvalues_ks_test
