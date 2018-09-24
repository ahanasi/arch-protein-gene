from load_RNAseq_counts import load_RNAseq_counts
from make_pairwise_correlation_matrix import make_pairwise_correlation_matrix
from plot_correlation_matrix import plot_correlation_matrix
from tabulate import tabulate
import numpy as np
import matplotlib.pyplot as pyplot
import seaborn as sns


def main():
    print 'loading counts'
    counts_matrix,gene_order,rep_list= load_RNAseq_counts('ARCH_RNAseq_concatenated_counts_2_19_2016.txt',delimiter=' ')

    
    #Print information about the dimensions of your counts array and the number of independent samples
    print np.shape(counts_matrix), "size of counts array"
    print np.shape(gene_order), "size of gene order array"
    print rep_list, "Here are the replicates in this file\n\n"

    #Here, for every file, compute descriptive statistics and plot histograms and boxplots

    soc = ['Sum of Counts']
    mean = ['Mean']
    pop_var = ['Population Variance']
    samp_var = ['Sample Variance']
    median = ['Median']
    sd = ['Standard Deviation']
    max_count = ['Maximum Count']
    min_count = ['Minimum Count']
    for j in range(0,np.shape(counts_matrix)[1]): 
        #print j, np.shape(counts_matrix[:,j])       
        soc.append(np.sum(counts_matrix[:,j]))#, "here is the sum of counts in sample ", rep_list[j]     
        mean.append(np.mean(counts_matrix[:,j]))# "here is the mean of counts in sample ", rep_list[j]
        pop_var.append(np.var(counts_matrix[:,j]))#, "here is the population variance of counts in sample ", rep_list[j]
        samp_var.append(np.var(counts_matrix[:,j], ddof = 1))#e, "here is the sample variance of counts in sample ", rep_list[j]
        median.append(np.median(counts_matrix[:,j]))#, "here is the median of counts in sample ", rep_list[j]
        sd.append(np.std(counts_matrix[:,j]))#, "here is the sample standard deviation of counts in sample ", rep_list[j]
        max_count.append(np.max(counts_matrix[:,j]))#, "here is the maximum count in sample ", rep_list[j]
        min_count.append(np.min(counts_matrix[:,j]))#, "here is the minimum count in sample ", rep_list[j]


##        #Create a counts histogram for all four samples
##        pyplot.figure()
##        pyplot.hist(counts_matrix[:,j],2000)
##        pyplot.xlim(0,5000)
##        pyplot.title('Frequency distribution for counts in '+rep_list[j])
##        pyplot.xlabel('Genes')
##        pyplot.ylabel('Frequency')
##        pyplot.show() 
##
##        #Create a probability distribution for all four samples
##        pyplot.figure()
##        pyplot.hist(counts_matrix[:,j],2000,normed=True)
##        pyplot.xlim(0,5000)
##        pyplot.title('Probability distribution for counts in '+rep_list[j])
##        pyplot.xlabel('Genes')
##        pyplot.ylabel('Probability')
##        pyplot.show()
##
##        #Create a boxplot for all four samples
##        pyplot.figure()
##        pyplot.boxplot(counts_matrix[:,j])
##        pyplot.ylim(0,20)
##        pyplot.title('Distribution for counts in '+rep_list[j])
##        pyplot.xlabel(rep_list[j])
##        pyplot.ylabel('Gene Counts')
##        pyplot.show()

    # Display descriptive statistics in the form of a table
    data1 = [soc,mean,pop_var,samp_var,median,sd,max_count,min_count]
    print tabulate(data1,headers = ['Statistic',rep_list[0],rep_list[1],rep_list[2],rep_list[3],rep_list[4],rep_list[5],rep_list[6],rep_list[7],rep_list[8],rep_list[9]],tablefmt='orgtbl')


##    #Here, create a box plot of counts across all four samples
##    data = [counts_matrix[:,0],counts_matrix[:,1],counts_matrix[:,2],counts_matrix[:,3]]
##    pyplot.figure()
##    pyplot.boxplot(data)
##    pyplot.ylim(0,20)
##    pyplot.title('Probability distribution for counts in all samples')
##    pyplot.ylabel('Gene Counts')
##    pyplot.show()
##
##    # Create scatter plots to explore associations between samples
##    for i in range (0,3):
##            pyplot.figure()
##            pyplot.scatter(counts_matrix[:,0],counts_matrix[:,i+1])
##            pyplot.title('Scatterplot for counts in '+rep_list[i+1]+' vs. '+rep_list[0])
##            pyplot.xlabel(rep_list[0])
##            pyplot.ylabel(rep_list[i+1])
##            pyplot.show()
##
##
##    for i in range (1,3):
##            pyplot.figure()
##            pyplot.scatter(counts_matrix[:,1],counts_matrix[:,i+1])
##            pyplot.title('Scatterplot for counts in '+rep_list[i+1]+' vs. '+rep_list[1])
##            pyplot.xlabel(rep_list[1])
##            pyplot.ylabel(rep_list[i+1])
##            pyplot.show()
##
##    for i in range (2,3):
##            pyplot.figure()
##            pyplot.scatter(counts_matrix[:,2],counts_matrix[:,i+1])
##            pyplot.title('Scatterplot for counts in '+rep_list[i+1]+' vs. '+rep_list[2])
##            pyplot.xlabel(rep_list[2])
##            pyplot.ylabel(rep_list[i+1])
##            pyplot.show()
##
##    
##    print 'calculating correlations'
##    pearson_matrix = make_pairwise_correlation_matrix(counts_matrix,correlation='pearson',lower_triangle_only=True)
##    spearman_matrix = make_pairwise_correlation_matrix(counts_matrix,correlation='spearman',lower_triangle_only=True)
##    kendall_tau_matrix = make_pairwise_correlation_matrix(counts_matrix,correlation='kendalltau',lower_triangle_only=True)
##
##    print 'plotting correlation heatmap'
##    plot_correlation_matrix(pearson_matrix,'pearson_correlation.png',lower_triangle_only=True,label_values=rep_list)
##    plot_correlation_matrix(spearman_matrix,'spearman_correlation.png',lower_triangle_only=True,label_values=rep_list)
##    plot_correlation_matrix(kendall_tau_matrix,'kendall_tau_correlation.png',lower_triangle_only=True,label_values=rep_list)

if __name__ == '__main__':
    main()
