#! /soe/mollyzhang/anaconda/bin/python
"""
this program reads difference postions
and plots the distribution of distances between the differences

Molly Zhang(mollyzhang@ucsc.edu)
Mon Aug  4 20:14:32 PDT 2014
"""

import sys
from collections import Counter
import matplotlib.pyplot as plt

def main():
    """ this function reads a file difference locations and plots
    the distribution of distances between the differences
    """ 

    num1 = 0
    num2 = 0
    distances = []

    f = open("chr21_test","r")
    #calculate the distances between each location and store in a list    
    for line in f:
        items = line.split()        
        position = int(items[1])
        num1 = num2
        num2 = position
        distances.append(num2-num1)
        
    count = Counter(distances)
    x = count.keys()[0:20]
    y = count.values()
    percent_y = []    
    for i in y[0:20]:    
        percent_y.append(i/float(sum(y))*100)
    
    plt.plot(x,percent_y)
    plt.xlabel('distance between SNPs')
    plt.ylabel('percentange occurance')
    plt.title('distribution of distances between SNPs in chr21 after clean')
    plt.show()



if __name__ == "__main__":
    sys.exit(main())
