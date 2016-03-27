#! /soe/mollyzhang/anaconda/bin/python
"""read_table_stats.py
this program reads in a tab delimited read labeling file
outputs the following:
1. number of total differences labeled in reads
2. number and percentage of non-conflicting differences labeled in reads
3. number and percentage of read pairs labeled
4. number and percentage of read pairs labeled as chimp
5. number of percentage of read pairs labeled as bonobo
6. add label to unambiguous reads as fourth column

Created by Molly Zhang on 10/16/14 at 21:59
"""
import sys
import argparse
import re

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("read_table_stats.py")
    parser.add_argument("f_catalogue", help = "the catalog file to start with")
    parser.add_argument("f_labeled", help = "the labled catalog file to end up with")
    return parser.parse_args()

def main(args):
    """Main method"""
    args = parse_args()
    total_read_pairs = 192365898.0
    f = open(args.f_catalogue, "r")
    f_labeled = open(args.f_labeled,"w")
    nc = 0  # number of read pairs labeled exclusively as chimp
    nb = 0  # number of read pairs labeled exclusively as bonobo
    n_both = 0  # number of read pairs labeled as both
    difs = set()
    difs_w_conflict = set()

    for line in f:
        # split with either \t or \n as delimiator
        # will always generate an empty string as the last item
        item = re.split("\t|\n", line)
                
        # seperate chimp differece positions and save to difs if not already there
        cp = item[1].split(",")
        for position in cp:
            difs.add(position)
                
        # seperate bonobo differece positions and save to 
        bp = item[2].split(",")
        for position in bp:
            difs.add(position)
 
        if cp != [''] and bp == ['']:
            item[3] = "c"
            nc += 1
            
        elif cp == [''] and bp != ['']:
            item[3] = "b"
            nb += 1 
            
        elif cp != [''] and bp != ['']:
            n_both += 1
            item[3] = "conflict"

            for position in cp:
                difs_w_conflict.add(position)
            
            for position in bp:
                difs_w_conflict.add(position)

        f_labeled.write(item[0] + "\t" + item[1] + "\t" + item[2] + "\t" 
                        + item[3] + "\n")    
        
    print ("number of read pairs labeled exlusively as chimp : %d" %nc +
           ", %.2f" %(nc/total_read_pairs*100) + "%")
    print ("number of read pairs labeled exlusively as bonobo :%d" %nb +
           ", %.2f" %(nb/total_read_pairs*100) + "%")
    print ("number of read pairs labeled as both chimpand bonobo : %d" %n_both +
           ", %.2f" %(n_both/total_read_pairs*100) + "%")
    print "number of difference positions recorded in the table: %d" %len(difs)
    print ("number of difference positions with conflict: %d" %len(difs_w_conflict) +
           ", %.2f" %(len(difs_w_conflict)/float(len(difs))*100) + "%")

if __name__ == "__main__" :
    sys.exit(main(sys.argv))

