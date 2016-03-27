#! /soe/mollyzhang/anaconda/bin/python

"""dif_cov.py
Created by Molly Zhang on 08/13/14 at 15:29

this program calculates the coverage of a bam file at each difference 
location between chimpanzee (panTro4) and bonobo (panpan1) reference

it return a histogram distribution of the coverage
and it also returns the upper 95% cutoff value for coverage   

"""

import sys
import argparse
import pysam


def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("calcaulate coverage at differences")
    parser.add_argument("f_dif", help = "the difference file")
    parser.add_argument("f_samfile", help = "the samfile")
    
    return parser.parse_args()

def main(args):
    """Main method"""
    args = parse_args()

    samfile = pysam.Samfile(args.f_samfile, "rb")
    f_dif = open(args.f_dif, "r")
    cov = dict()
    # loop through difference file and save the coverage at each difference
    # to the cov dictionary
    for difference in f_dif:
        items = difference.strip().split()
        chrm = items[0] # chrm is the chromosome number
        p_1b = int(items[1]) # p_1b stands for position 1-based

        # get number of reads mapped to each difference using pileup()
        for pileupcolumn in samfile.pileup(chrm, p_1b -1, p_1b,
                                          truncate = True):
            n_reads = pileupcolumn.n
            
            # add the number of reads to the cov dictionary
            if n_reads in cov.keys():
                cov[n_reads] += 1
            else: 
                cov[n_reads] = 1
    
    for keys in cov:
        print keys, cov[keys]

if __name__ == "__main__" :
    sys.exit(main(sys.argv))

