#! /soe/mollyzhang/anaconda/bin/python
"""label_stats.py
this program will read a bam file and return the following statistics:
--what's the distribution of number of differences mapped to each read.
--if one read is labeled, the chance its read pair is labeled as well
--percentage of reads that none of the pairs are labeled
--out of labeled pairs, what's the percentage of contradciting label
--Created by Molly Zhang on 09/09/14 at 15:46
"""
import sys
import argparse
import pysam

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("label_stats.py")
    parser.add_argument("bam_file", help = "the bam file to analyze")
    return parser.parse_args()

def main(args):
    """Main method"""
    #args = parse_args()
    f_bam = pysam.Samfile("../chonobo_data/chr21.qual_filter.labeled.r1.bam", "rb")
    sys.stdout.write("distribution of differences mapped to each labeled read: \n")
    dif_distribution = dict()
    for read in f_bam:        
        labels = read.opt('cp') + read.opt('bp')
        num_difs = len(labels.split(',')) - 1        
        try: dif_distribution[num_difs] += 1
        except: dif_distribution[num_difs] = 1
    sys.stdout.write(str(dif_distribution) + "\n")

if __name__ == "__main__" :
    sys.exit(main(sys.argv))

