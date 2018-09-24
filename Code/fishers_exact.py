from scipy import stats
import numpy as np
from scipy.stats import fisher_exact
import matplotlib.pyplot as pyplot
import statsmodels.api as sm
import statsmodels.sandbox


def fishers_exact(file_string):


    normed = np.loadtxt(file_string, usecols=(1, 2, 3, 4, 5, 6, 7, 9, 9, 10), dtype = 'float', skiprows=1)
    normed = normed.astype(int)
    gene_name = np.loadtxt(file_string, usecols=(0,), skiprows=1, dtype = 'string')
    column_sum_wt1=np.sum(np.sum(normed[:,0]))
    column_sum_wt2=np.sum(np.sum(normed[:,1]))
    column_sum_wt3=np.sum(np.sum(normed[:,2]))
    column_sum_ctcf4=np.sum(np.sum(normed[:,3]))
    column_sum_ctcf5=np.sum(np.sum(normed[:,4]))
    column_sum_ctcf6=np.sum(np.sum(normed[:,5]))
    column_sum_med7=np.sum(np.sum(normed[:,6]))
    column_sum_med8=np.sum(np.sum(normed[:,7]))
    column_sum_smc9=np.sum(np.sum(normed[:,8]))
    column_sum_smc10=np.sum(np.sum(normed[:,9]))


    print normed
    print np.shape(normed)[1],np.shape(normed)[0]

    
    p_values = np.zeros((np.shape(normed)[0], 1))
    pvalues_fishers_exact = np.zeros((np.shape(normed)[0], 2))

    fisher_stats = np.zeros((np.shape(normed)[0],4))


    for i in range(0,(np.shape(normed)[0])):
        
        #Wild Type
        fisher_stats[i,0] = np.mean(normed[i,0:3])
        
        #CTCF
        fisher_stats[i,1] = np.mean(normed[i,3:6])
        
        #MED
        fisher_stats[i,2] = np.mean(normed[i,6:8])
        
        #SMC
        fisher_stats[i,3] = np.mean(normed[i,8:10])
        
        #For WT vs MED
        
        oddsratio, pvalue = stats.fisher_exact([[np.sum(normed[i,0:3]),np.sum(normed[i,3:6])],[((column_sum_wt1+column_sum_wt2+column_sum_wt3)-np.sum(normed[i,0:3])),((column_sum_ctcf4+column_sum_ctcf5+column_sum_ctcf6)-np.sum(normed[i,3:6]))]])
        #print oddsratio, pvalue, "Fisher's Exact pvalue and Odd's ratio"
        p_values[i,0] = pvalue
        #pvalues_fishers_exact[i,0] = -10*np.log2(pvalue)
        #pvalues_fishers_exact[i,1] = np.log2(np.divide((fisher_stats[i,3]+1),((fisher_stats[i,0]+1))))
        
    
    
##    pvalues_corrected = statsmodels.sandbox.stats.multicomp.multipletests(p_values[:,0], alpha=0.05, method='fdr_bh')
##    p_corr = pvalues_corrected[1]
##        
##    #Find sig increse and decrease
##    
##    inc = 0
##    dec = 0
##    for i in range(0,(np.shape(normed)[0])):    
##        if (p_values[i,0] <= (0.05)):
##            if fisher_stats[i,0] < fisher_stats[i,3]:
##                inc += 1
##            else:
##                dec += 1
##    
##    print(inc,dec)
##
##    #Find sig increse and decrease after correction
##    print 'After Correction'
##    inc = 0
##    dec = 0
##    for i in range(0,(np.shape(normed)[0])):    
##        if (p_corr[i] <= (0.05)):
##            if fisher_stats[i,0] < fisher_stats[i,3]:
##                inc += 1
##            else:
##                dec += 1
##    
##    print(inc,dec)
        

    #Make Histogram of Fisher's Exact p-values
    
    

##    pyplot.figure()
##    pyplot.hist(pvalues_fishers_exact[:,0], density = True, log = True)
##    pyplot.xlim(0,1)
##    pyplot.xlabel('p-values')
##    pyplot.ylabel('Frequency')
##    pyplot.title('WT vs. SMC')
##    pyplot.show()

    pyplot.figure()
    pyplot.hist(p_values[:,0], density = True)
    pyplot.xlim(0,1)
    pyplot.xlabel('p-values')
    pyplot.ylabel('Frequency')
    pyplot.title('WT vs. CTCF')
    pyplot.show()

    pyplot.figure()
    pyplot.hist(p_values[:,0])
    pyplot.xlim(0,1)
    pyplot.xlabel('p-values')
    pyplot.ylabel('Frequency')
    pyplot.title('WT vs. CTCF')
    pyplot.show()

##    pyplot.figure()
##    pyplot.hist(pvalues_fishers_exact[:,0])
##    pyplot.xlim(0,1)
##    pyplot.xlabel('p-values')
##    pyplot.ylabel('Frequency')
##    pyplot.title('WT vs. SMC')
##    pyplot.show()

    
##    #Make MA plot to look at fold change with respect to abundance (Case vs WT)
##        
##    m = ((fisher_stats[:,0]) - (fisher_stats[:,3]))
##    a = 0.5 * (np.log2(fisher_stats[:,3]+1) + np.log2(fisher_stats[:,0]+1))
##    pyplot.figure()             
##    pyplot.scatter(a,m)
##
##    pyplot.scatter(a[np.where(pvalues_fishers_exact[:,0]>=-10*np.log2(0.05))],m[np.where(pvalues_fishers_exact[:,0]>=-10*np.log2(0.05))],color='red')
##    pyplot.title('MA Plot: WT vs. SMC - pvalues from Fishers Exact Test Code')
##    pyplot.show()
##
##
##  # Make the volcano plots to look at pvalues with respect to fold change
##
##    pyplot.figure()
##    pyplot.scatter(pvalues_fishers_exact[:,1],pvalues_fishers_exact[:,0])
##    pyplot.scatter(pvalues_fishers_exact[np.where(pvalues_fishers_exact[:,0]>=-10*np.log2(0.05)),1],pvalues_fishers_exact[np.where(pvalues_fishers_exact[:,0]>=-10*np.log2(0.05)),0],color='red')
##    pyplot.title('Volcano plot: WT vs. SMC - pvalues from Fishers Exact Test Code')
##    pyplot.show()
                                   

    return pvalues_fishers_exact
			
			
