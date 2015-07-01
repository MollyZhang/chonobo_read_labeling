#! /soe/mollyzhang/anaconda/bin/python
"""sub_sampling_dif.py
pick designated differences between chimp/bonobo reference by 
picking one row out of every n rows
save both the subsample and the leftover after subsampling in new files 
Created by Molly Zhang on 12/17/14 at 15:37
"""
import sys
import argparse

def parse_args():
    parser = argparse.ArgumentParser("sub_sampling_dif.py")
    parser.add_argument("n", type=int, 
        help="pick one row out of every n-th row")
    parser.add_argument("subsample_filename")
    parser.add_argument("subsample_leftover_filename")
    return parser.parse_args()


def main(args):
    args = parse_args()
    f = open("../clean_dif_catalog/pp1_pT4_difs_cleaned.baseq-c.tab", "r")
    f_new = open(args.subsample_filename, "w")
    f_leftover = open(args.subsample_leftover_filename, "w")
    nline = 1 
    for line in f:
        if nline%args.n == 1:
            f_new.write(line)
        else:
            f_leftover.write(line)
        nline += 1


if __name__ == "__main__" :
    sys.exit(main(sys.argv))

