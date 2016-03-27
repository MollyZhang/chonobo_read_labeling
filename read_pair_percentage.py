#! /soe/mollyzhang/anaconda/bin/python
"""read_pair_percentage.py
Created by Molly Zhang on 09/15/14 at 21:29
"""
import sys
import argparse

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("read_pair_percentage.py")
    parser.add_argument("folder", help = "folder where all the bam files are located")
    return parser.parse_args()

def main(args):
    """Main method"""
    args = parse_args()
    chroms = ['1', '2A', '2B', '3', '4', '5', '6', '7', '8', '9', '10', '11',
              '12', '13', '14', '15', '16', '17', '18', '19', '20', '21',
              '22', 'X']
    
    num_reads = 24 * ['', 0, 0]
    # 2D array num_reads saves all results
    # num_reads[i][0]: chromosome name
    # num_reads[i][1]: number of all reads
    # num_reads[i][2]: number of reads that are not paired
    i = 0
    for chrm in chroms:
        f = pysam.Samfile(args.folder/args.folder + "/chr" + chrm + ".labeled.r1.bam", "rb")
        num_reads[i][0] = "chr" + chrm
        num_reads[i][1] = f.mapped + f.unmapped
        for read in f:
            if read.is_paired: 
                num_reads[i][2] += 1
        f.close()
        i += 1





if __name__ == "__main__" :
    sys.exit(main(sys.argv))

