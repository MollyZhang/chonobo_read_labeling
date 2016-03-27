#! /soe/mollyzhang/anaconda/bin/python
"""faSize_calc.py.py
this script calculate the percentage lengths of references that are
labeled as chr1-chrX, not counting the ones that looks like: 
chr10_AACZ03166820_random

Created by Molly Zhang on 09/17/14 at 12:15
"""
import sys
import argparse

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("faSize_calc.py")
    parser.add_argument("file", help = "the faSize -detail -tab output")
    return parser.parse_args()

def main(args):
    """Main method"""
    args = parse_args()
    chroms = ['1', '2A', '2B', '3', '4', '5', '6', '7', '8', '9', '10', '11', 
    '12', '13', '14', '15', '16', '17', '18', '19', '20', '21','22', 'X']
    for i in range(len(chroms)):
        chroms[i] = "chr" + chroms[i]
    length_chrms = 0
    length_non_chrms = 0
    f = open(args.file, "r")    
    for line in f:
        items = line.strip().split()
        if items[0] not in chroms:
            length_non_chrms += int(items[1])
        else:
            length_chrms += int(items[1])
    percentage = length_chrms/float(length_chrms + length_non_chrms)
    print "total length of bases in intact chromosomes: %d" %length_chrms
    print "total length of bases in randome scanfolds:  %d" % length_non_chrms
    print "percentage of intact chromosomes: %f5.2" %percentage
        
        
if __name__ == "__main__" :
    sys.exit(main(sys.argv))

