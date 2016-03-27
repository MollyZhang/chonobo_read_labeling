#! /soe/mollyzhang/anaconda/bin/python
"""find_het_site_dif.py
after labeling heterozygous sites, this scripts find the difference between:
A. previously labeled difference sites from chimp and bono reference
B. the part of newly heterozygous sites that overlaps with A

Created by Molly Zhang on 12/17/14 at 14:42
"""
import sys
import argparse

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("find_het_site_dif.py")
    parser.add_argument("labeled_het", help = "the newly labeled het sites")
    parser.add_argument("original_dif", help = "the original difference file")
    return parser.parse_args()

def main(args):
    """ Main method (very descriptive isn't it)"""
    args = parse_args()
    f_het = open(args.labeled_het, "r")
    f_dif = open(args.original_dif, "r")

    het = dict()
    dif = dict()
    
    for line in f_het:
        items = line.strip().split("\t")
        key = items[0] + ":" + str(items[1])
        value = items[2] + items[3]
        het[key] = value
    f_het.close()

    for line in f_dif:
        items = line.strip().split("\t")  
        key = items[0] + ":" + str(items[1])
        value = items[2] + items[3]
        dif[key] = value    
    
    f_dif.close()

    print "items in het dictionary: %d" %len(het.viewitems())
    print "items in dif dictionary: %d" %len(dif.viewitems())
  
    labeled_dif = set() # dif sites that can be found in het sites
    correct_labeled_dif = set() # correctly labeled dif sites
    wrong_labeled_dif = set() # wrongly labeled dif sites
    un_labeled_dif = set() # difs that are not found in het dictionary
    all_het = set()
    dif_not_in_all_het = set()
    
    f = open("../all_heterozygous_sites/diana.all_heterozygous_sites", "r")
    for line in f:
        items = line.strip().split("\t")  
        key = items[0] + ":" + str(items[1])
        all_het.add(key)
    
    for key in dif:
        if key in het:
            labeled_dif.add(key)
            if dif[key] == het[key]:
                correct_labeled_dif.add(key)
            else:
                wrong_labeled_dif.add(key)
        else:
            un_labeled_dif.add(key)

        if key not in all_het:
            dif_not_in_all_het.add(key)

    print "dif sites that can be found in het sites: %d" %len(labeled_dif)
    print "correctly labeled dif sites in het: %d" %len(correct_labeled_dif)
    print "wrongly labeled dif sites in het: %d" %len(wrong_labeled_dif)
    print "dif sites not found in labeled het sites: %d" %len(un_labeled_dif)
    print "dif sites not found in all het sites: %d" %len(dif_not_in_all_het)
    difs = un_labeled_dif.symmetric_difference(dif_not_in_all_het)
    print "difference between the above two sets: %d" %len(difs)

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
