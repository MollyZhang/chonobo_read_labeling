#! /usr/bin/python
#
# Based on a perl script written by Ed, we take a fasta file from input and
# outputs the sequences, split into smaller sequences of the specified 
# length and overlap. Also adding a feature that does not output sequences
# that are all N's
#
# Sam Vohr (svohr@soe.ucsc.edu)
# Fri Feb 14 15:18:17 PST 2014

import sys
import argparse
import re

import fasta

n_pat = re.compile('^NN*$')

def process_sequences( fasta_in, fasta_out, win_size, overlap, \
                       root_seq_id=None, keep_ns=False ):
   
    root_id = root_seq_id
    for seq_id, comment, seq in fasta.read_seq( fasta_in ):
        if root_seq_id is None:
            root_id = seq_id
        
        pos = 1
        win = 1

        seq_len = len(seq)
        while pos < seq_len:
            win_id = '%s.%d' % (root_id, win)
            end = pos + win_size - 1
            if seq_len < (pos + win_size):
                end = seq_len
            desc = '%d to %d' % (pos,end)

            win_seq = seq[ pos-1 : end ]
            
            if not keep_ns:
                all_ns = re.match( n_pat, win_seq )
                if all_ns is None:
                    fasta.write_seq( fasta_out, win_id, desc, win_seq )
            else:
                fasta.write_seq( fasta_out, win_id, desc, win_seq )
            
            pos += (win_size - overlap)
            win += 1
    


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser("Split and tile sequences by size and overlap.")
    parser.add_argument('-f', '--file', dest='in_file', type=str, \
                        metavar='in.fa', help="Fasta input file.")
    parser.add_argument('-i', '--id', dest='root_id', type=str, \
                        metavar='root_id', help="Fasta input file.")
    parser.add_argument('-w', '--window', dest='win_size', default=10000, \
                        type=int, metavar='win_size', help="Window size.")
    parser.add_argument('-o', '--overlap', dest='overlap', default=0, \
                        type=int, metavar='overlap', help="Overlap size.")
    parser.add_argument('-n', dest='keep_ns', default=False, \
                        action="store_true", \
                        help="Output subsequences of all Ns.")
    args = parser.parse_args()

    fasta_in = sys.stdin
    if args.in_file is not None:
        fasta_in = open(args.in_file,'r')
    
    fasta_out = sys.stdout

    process_sequences( fasta_in, fasta_out, args.win_size, args.overlap, \
                       root_seq_id=args.root_id, keep_ns=args.keep_ns )

    fasta_in.close()

    return 0

if __name__ == "__main__":
    sys.exit(main())
