#! /soe/mollyzhang/anaconda/bin/python

"""
venn_diagram.py

Molly Zhang
Thu Aug 14 16:17:20 PDT 2014

this program takes three sam files:
1) chonobo reads mapped to chimpanzee reference
2) chonobo reads mapped to human reference
3) chonobo reads mapped to bonobo reference

and then it saves a venn diagram of the intersection
"""

from matplotlib_venn import venn3
from matplotlib import pyplot as plt
import pysam
import collections
import argparse
import sys

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("calculates read mapping intersection")
    parser.add_argument("f_chimp", help = 
                        "samfile of chonobo reads mapped to chimp reference")
    parser.add_argument("f_bono", help = 
                        "samfile of chonobo reads mapped to bonobo referece")
    parser.add_argument("f_human", help = 
                        "samfile of chonobo reads mapped to human reference")
    parser.add_argument("plot_title", help = "title of the plot")
    return parser.parse_args()


def main():
    """ main function """
    args = parse_args()

    f_chimp = pysam.Samfile(args.f_chimp, "rb")
    f_bono = pysam.Samfile(args.f_bono, "rb")
    f_human = pysam.Samfile(args.f_human, "rb")
    
    # stores the ID of each mapped reads in each samfile in a list
    reads_chimp = list()
    reads_bono = list()
    reads_human = list()

    for read in f_chimp:
        reads_chimp.append(read.qname)

    for read in f_bono:
        reads_bono.append(read.qname)

    for read in f_human:
        reads_human.append(read.qname)

    # find out if there is duplicate in the read IDs
    




    # find the intersections between the three
    overlap_bono_chimp = set(reads_bono).intersection(set(reads_chimp))
    overlap_bono_human = set(reads_bono).intersection(set(reads_human))
    overlap_chimp_human = set(reads_chimp).intersection(set(reads_human))
    overlap_bono_chimp_human = (set(reads_bono).
                                intersection(set(reads_chimp)).
                                intersection(set(reads_human)))

    #venn3() takes a list of 7 numbers:
    #venn[6] -> number of reads mapped to human, chimp and bonobo  
    #venn[5] -> number of reads mapped to chimp and human but not bonobo
    #venn[4] -> number of reads mapped to bonobo and human but not chimp
    #venn[3] -> number of reads mapped to chimp and bonobo but not human
    #venn[2] -> number of reads mapped only to human
    #venn[1] -> number of reads mapped only to chimp
    #venn[0] -> number of reads mapped only to bonobo

    venn = [0]*7
    venn[6] = overlap_bono_chimp_human
    venn[5] = overlap_chimp_human - overlap_bono_chimp_human
    venn[4] = overlap_bono_human - overlap_bono_chimp_human
    venn[3] = overlap_bono_chimp - overlap_bono_chimp_human
    venn[2] = f_human.mapped - venn[4] - venn[5] - overlap_bono_chimp_human
    venn[1] = f_chimp.mapped - venn[3] - venn[5] - overlap_bono_chimp_human
    venn[0] = f_bono.mapped - venn[3] - venn[4] - overlap_bono_chimp_human

    #plot and save venn diagrams   
    venn3(subsets=venn, set_labels = ("Bonobo", "Chimpanzee", "Human"))
    plt.title(args.plot_title)
    plt.show()

if __name__ == "__main__" :
    sys.exit(main())

