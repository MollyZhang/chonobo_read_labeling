#! /soe/mollyzhang/anaconda/bin/python
"""make_histogram.py

this program will plot a file's first two columns
Created by Molly Zhang on 08/21/14 at 17:32
"""

import sys
import argparse
import numpy
from matplotlib import pyplot as plt


def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("make_histogram.py")
    parser.add_argument("data_file", help = "the file with desired data")
    parser.add_argument("-min", "--min_x_axis", type = int, default = 0,
                        help = "the minimum x axis number to display")
    parser.add_argument("-max", "--max_x_axis", type = int, default = 10000, 
                        help = "the maximum x axis number to display") 
    parser.add_argument("-a", "--accumulative", action = "store_true", 
                        help = "display accumulative histogram or not") 
    parser.add_argument("-t", "--title", default = "coverage distribution", 
                        help = "title of the plot")
    parser.add_argument("-x", "--xaxis", default = "coverage", 
                        help = "x axis title")
    parser.add_argument("-y", "--yaxis", default = "percentage",
                        help = "y axis title")
    return parser.parse_args()

def main(args):
    """Main method"""
    args = parse_args()
    
    x = list()
    y = list()
    y_accumulative = list()
    y_percentile = list()
    y_accu_percentile = list()

    f = open(args.data_file, "r")
 
    for line in f:
        items = line.strip().split()
        
        x.append(int(items[1]))
        y.append(int(items[0]))
        y_accumulative.append(sum(y))
   
 
    # calculate the percentile list and accumulative percentile list

    sum_all = float(y_accumulative[-1])    
    for i in y:
        y_percentile.append( i / float(sum_all) * 100)
    for i in y_accumulative:
        y_accu_percentile.append( i / float(sum_all) * 100)
            

    if args.accumulative:
        plt.plot(x[args.min_x_axis:args.max_x_axis], 
                 y_accu_percentile[args.min_x_axis:args.max_x_axis])
    else:   
        plt.plot(x[args.min_x_axis:args.max_x_axis],
                 y_percentile[args.min_x_axis:args.max_x_axis])
    plt.xlabel(args.xaxis)
    plt.ylabel(args.yaxis)
    
    upper_limit = calculate_percentage_limit(0.975, x, y_accumulative)
    lower_limit = calculate_percentage_limit(0.025, x, y_accumulative)
    plt.title(args.title)
    plt.show()


def calculate_percentage_limit(percent, x, y_accumulative):
    """this method calculates the 97.5% percentile of the coverage"""
    limit = percent * float(y_accumulative[-1])
    for i in range(len(y_accumulative)):
        if y_accumulative[i] > limit:
            index = i
            break
    return x[index] 
 
if __name__ == "__main__" :
    sys.exit(main(sys.argv))

