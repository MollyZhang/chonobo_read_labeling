#! /usr/bin/python
""" 
find_difference.py
This program reads input from samtools mpileup and writes only the
positions that differs between chimp and bonobo.
We only look at sites with coverage 1.

Molly Zhang (mollyzhang@ucsc.edu)
Thu Jul 17 11:46:19 PDT 2014
"""
import sys

def main():
    """ this function reads input from samtools mpileup and writes only the
    positions that differs between chimp and bonobo"""
    for line in sys.stdin:
        items = line.split()
        chrom_num = items[0]
        position = items[1]
        ref_base = items[2].upper()
        coverage = int(items[3])
        if coverage == 1 and ref_base != 'N':
            base = items[4].upper()
            if is_snp(base):
                sys.stdout.write("%s\t%s\t%s\t%s\n" % (chrom_num, position,
                                 ref_base, extract(base)))
    return 0

# to calculate the difference of the snps with >30 mapping quality:
# samtools mpileup -q 30 -f /projects/redser/ftp/UCSC/bwa-indeces/panTro4.fa \ 
# panpan1_256_tile.pt4.s.bam | ./find_difference.py | less
        
def is_snp(string):
    """this is a function to tell if a string is a snp or not"""
    if string[0].upper() in 'AGTC':
        return True
    elif string[0] == '^':
        if string[2].upper() in 'AGTC':
            return True
        else:
            return False
    else:
        return False

def extract(string):
    """ this is a function to get rid of ^F or $ in the samtools mpileup output     and only to return the desired base"""
    if len(string) == 1:
        return string
    elif len(string) == 2 and string[1] == '$': # end of reads looks like G$
        return string[0] 
    elif len(string) == 3 and string[0] == '^': # begining of reads is like ^FG
        return string[2]
    else: 
        return string

if __name__ == "__main__":
    sys.exit(main())

