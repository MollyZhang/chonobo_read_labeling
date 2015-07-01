#! /soe/mollyzhang/anaconda/bin/python
"""all_heter_sites_pysam.py

use pysam to select for all the diana heterozygous sites using 
the criteria: coverage >= 4 and coverage <=23
mapping quality of the read >= 30, base quality >= 20
at least one forward base and reverse base match one base
at least one forward base and reverse base match anothter base

Created by Molly Zhang on 09/17/14 at 15:59
"""
import sys
import argparse
import pysam

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("all_heter_sites_pysam.py")
    parser.add_argument("f_samfile", help = "the input samfile")
    return parser.parse_args()

def main(args):
    """Main method"""
    args = parse_args()
    samfile = pysam.Samfile(args.f_samfile, "rb")    
    for pileupcolumn in samfile.pileup(truncate = True):
        chonobo_bases = list() # store the aligned bases
        for pileupread in pileupcolumn.pileups:
            base_q = pileupread.alignment.qual[pileupread.qpos]
            # check for if aligned read's mapq and base_q
            if (pileupread.alignment.mapq >= 30) and (ord(base_q) - 33 >= 20):
                base = pileupread.alignment.seq[pileupread.qpos]
                if pileupread.alignment.is_reverse:
                    base = base.lower() # reverse reads --> lower case
                else:
                    base = base.upper() # forward reads --> upper case
                chonobo_bases.append(base)
        if len(chonobo_bases) < 4 or len(chonobo_bases) > 23:
            continue
        elif not has_multiple_bases(''.join(chonobo_bases)):
            continue
        else:
            het_result = is_heter(''.join(chonobo_bases))
            if het_result[0]:
                print(samfile.getrname(pileupcolumn.tid) + "\t" + 
                      str(pileupcolumn.pos+1) + "\t" +
                      ''.join(het_result[1]) + "\t" + 
                      str(len(chonobo_bases)) + "\t" + 
                      ''.join(chonobo_bases))



def is_heter(bases):
    """ determine if a site is heterozygous, return a list
    first item in list: True or False (Tri-allelic is false)
    second and third item in list: two heterozygous bases
    """
    table = dict() # save the occurnace of all bases in a table

    #initialize table
    for i in "agtcAGTC":
        table[i] = 0

    #save number of occurances in the base string
    for x in bases:
        if x in "agtcAGTC":
            table[x] += 1

    s = [] #only save the bases that occured at least once 
    het_bases = []
    for base in table:
        if table[base] > 0:
            s.append(base)

    for letter in "AGTC":
        if (letter.lower() in s) and (letter.upper() in s):
            het_bases.append(letter)

    if len(het_bases) == 2:
        return True, het_bases
    else:
        return False, het_bases



def has_multiple_bases(bases):
    """ this method checks whether two more more different bases exist
    in te base string
    """
    bases = bases.upper()
    if ("A" in bases) and ("G" in bases):
        return True
    elif ("A" in bases) and ("T" in bases):
        return True
    elif ("A" in bases) and ("C" in bases):
        return True
    elif ("G" in bases) and ("T" in bases):
        return True
    elif ("G" in bases) and ("C" in bases):
        return True
    elif ("T" in bases) and ("C" in bases):
        return True
    else:
        return False













if __name__ == "__main__" :
    sys.exit(main(sys.argv))

