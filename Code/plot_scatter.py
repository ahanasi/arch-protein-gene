import numpy as np
import matplotlib.pyplot as pyplot

"""
This function primarily works to create a histogram(s) from an input data.
"""

def scatterplot(matrix,gene_order,lib_order):


    pyplot.figure()
    pyplot.scatter(matrix[:,0],matrix[:,3])    
    pyplot.title('Scatterplot comparing for counts in '+lib_order[0]+' vs. '+lib_order[3])
    pyplot.xlim(0,10000)
    pyplot.ylim(0,10000)
    pyplot.xlabel('Counts in ' + lib_order[0])
    pyplot.ylabel('Counts in ' + lib_order[3])
    pyplot.show()
    
    pyplot.figure()
    pyplot.scatter(matrix[:,0],matrix[:,4])    
    pyplot.title('Scatterplot comparing for counts in '+lib_order[0]+' vs. '+lib_order[4])
    pyplot.xlim(0,10000)
    pyplot.ylim(0,10000)
    pyplot.xlabel('Counts in ' + lib_order[0])
    pyplot.ylabel('Counts in ' + lib_order[4])
    pyplot.show()
    
    pyplot.figure()
    pyplot.scatter(matrix[:,0],matrix[:,5])    
    pyplot.title('Scatterplot comparing for counts in '+lib_order[0]+' vs. '+lib_order[5])
    pyplot.xlim(0,10000)
    pyplot.ylim(0,10000)
    pyplot.xlabel('Counts in ' + lib_order[0])
    pyplot.ylabel('Counts in ' + lib_order[5])
    pyplot.show()
    
    pyplot.figure()
    pyplot.scatter(matrix[:,1],matrix[:,6])
    pyplot.title('Scatterplot comparing for counts in '+lib_order[1]+' vs. '+lib_order[6])
    pyplot.xlim(0,10000)
    pyplot.ylim(0,10000)
    pyplot.xlabel('Counts in ' + lib_order[1])
    pyplot.ylabel('Counts in ' + lib_order[6])
    pyplot.show()
    
    pyplot.figure()
    pyplot.scatter(matrix[:,1],matrix[:,7])
    pyplot.title('Scatterplot comparing for counts in '+lib_order[1]+' vs. '+lib_order[7])
    pyplot.xlim(0,10000)
    pyplot.ylim(0,10000)
    pyplot.xlabel('Counts in ' + lib_order[1])
    pyplot.ylabel('Counts in ' + lib_order[7])
    pyplot.show()

    pyplot.figure()
    pyplot.scatter(matrix[:,2],matrix[:,8])
    pyplot.title('Scatterplot comparing for counts in '+lib_order[2]+' vs. '+lib_order[8])
    pyplot.xlim(0,10000)
    pyplot.ylim(0,10000)
    pyplot.xlabel('Counts in ' + lib_order[2])
    pyplot.ylabel('Counts in ' + lib_order[8])
    pyplot.show()
    
    pyplot.figure()
    pyplot.scatter(matrix[:,2],matrix[:,9])
    pyplot.title('Scatterplot comparing for counts in '+lib_order[2]+' vs. '+lib_order[9])
    pyplot.xlim(0,10000)
    pyplot.ylim(0,10000)
    pyplot.xlabel('Counts in ' + lib_order[2])
    pyplot.ylabel('Counts in ' + lib_order[9])
    pyplot.show()


   
