import numpy as np
import matplotlib.pyplot as pyplot

"""
This function primarily works to create MA plot from input data.
"""

def ma_plot(matrix,gene_order,lib_order):


    m = (matrix[:,3]-matrix[:,0])
    a = 0.5 * (np.log2(matrix[:,3]) + np.log2(matrix[:,0]))

    pyplot.figure()
    pyplot.scatter(a,m)    
    pyplot.title('MA '+lib_order[0]+' vs. '+lib_order[3])
    pyplot.show()


    m = (matrix[:,4]-matrix[:,0])
    a = 0.5 * (np.log2(matrix[:,4]) + np.log2(matrix[:,0]))

    pyplot.figure()
    pyplot.scatter(a,m)                    
    pyplot.title('MA '+lib_order[0]+' vs. '+lib_order[4])
    pyplot.show()


    m = (matrix[:,5]-matrix[:,0])
    a = 0.5 * (np.log2(matrix[:,5]) + np.log2(matrix[:,0]))

    pyplot.figure()
    pyplot.scatter(a,m)                    
    pyplot.title('MA '+lib_order[0]+' vs. '+lib_order[5])
    pyplot.show()



    m = (matrix[:,6]-matrix[:,1])
    a = 0.5 * (np.log2(matrix[:,6]) + np.log2(matrix[:,1]))

    pyplot.figure()
    pyplot.scatter(a,m)                    
    pyplot.title('MA '+lib_order[1]+' vs. '+lib_order[6])
    pyplot.show()


    m = (matrix[:,7]-matrix[:,1])
    a = 0.5 * (np.log2(matrix[:,7]) + np.log2(matrix[:,1]))

    pyplot.figure()
    pyplot.scatter(a,m)                    
    pyplot.title('MA '+lib_order[1]+' vs. '+lib_order[7])
    pyplot.show()


    m = (matrix[:,8]-matrix[:,2])
    a = 0.5 * (np.log2(matrix[:,8]) + np.log2(matrix[:,2]))

    pyplot.figure()
    pyplot.scatter(a,m)                    
    pyplot.title('MA '+lib_order[2]+' vs. '+lib_order[8])
    pyplot.show()
    
    m = (matrix[:,9]-matrix[:,2])
    a = 0.5 * (np.log2(matrix[:,9]) + np.log2(matrix[:,2]))

    pyplot.figure()
    pyplot.scatter(a,m)                    
    pyplot.title('MA '+lib_order[2]+' vs. '+lib_order[9])
    pyplot.show()


