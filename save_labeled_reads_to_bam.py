#! /soe/mollyzhang/anaconda/bin/python
"""save_labeled_reads_to_bam.py
take the read catalogue and original diana bam file and save the labeled reads 
(without conflict) to a bam file. Each read is affixed with a tag "cb", in which a "c" label means this read is labeled as chimp, and a "b" label means that this read is labeled as bonobo
Created by Molly Zhang on 11/02/14 at 15:29
"""
import sys
import argparse
import pysam


def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("save_labeled_reads_to_bam.py")
    parser.add_argument("f_bam", help = "the starting diana bam file")
    parser.add_argument("f_read", help = "the read catalog")
    parser.add_argument("f_labeled", help = "the reads that are unanimously labeled as either chimp or bonobo")
    return parser.parse_args()

def main(args):
    """Main method"""
    args = parse_args()
    
    # save all labeled reads into two different sets, chimp reads or bonobo reads
    saved_reads = reads_to_set(args.f_read)
    chimp_reads = saved_reads[0]
    bono_reads = saved_reads[1]

    # open the bam files to read and to write
    f_bam = pysam.Samfile(args.f_bam, "rb")
    f_labeled = pysam.Samfile(args.f_labeled, "wb", template=f_bam)
    
    # loop through each read in the original bam file and save read to corresponding sub bam-file
    for read in f_bam:
        read_name = read.qname
        if read_name in chimp_reads:
            read.setTag("cb", "c")
            f_labeled.write(read)
        elif read_name in bono_reads:
            read.setTag("cb", "b")
            f_labeled.write(read)

def reads_to_set(filename):
    """ save reads to either a set of chimp reads or a set of bonobo reads """

    f = open(filename, "r")
    chimp_reads = set()
    bono_reads = set()

    for line in f:
        items = line.strip().split("\t")
        if items[-1] == "c":
            chimp_reads.add(items[0])
        elif items[-1] == "b":
            bono_reads.add(items[0])
    f.close()
    return chimp_reads, bono_reads


if __name__ == "__main__" :
    sys.exit(main(sys.argv))

