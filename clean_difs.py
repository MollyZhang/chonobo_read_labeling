#! /soe/mollyzhang/anaconda/bin/python 

"""this program pileup the chonobo reads mapped to each difference specified
in the chimp-bono difference tab file based on the base of chonobo reads , then
this program evaluates whether the difference meets the clean up criteria, and
it only prints out the differences that meet the criteria.

Molly zhang (mollyzhang@ucsc.edu)
Thu Jul 24 15:13:50 PDT 2014
"""
import pysam
import sys
import argparse

def main():

    #argparse: to add command line options
    parser = argparse.ArgumentParser(
    "clean_disf.py")
    parser.add_argument("f_dif", help="the difference file")
    parser.add_argument("f_samfile", help="the samfile")
    parser.add_argument("-c", "--chonobo", default=False,
                        action = "store_true",
                        help = "print out chonobo coverage and bases")
    args = parser.parse_args()
    
    samfile = pysam.Samfile(args.f_samfile, "rb")
    f_dif = open(args.f_dif, "r")
    for difference in f_dif:
        """"loop through each difference and return the relevant information
        for each difference """
        items = difference.strip().split()
        chrm = items[0] # chrm is the chromosome number
        p_1b = int(items[1]) # p_1b stands for position 1-based
        chimp_base = items[2]
        bono_base = items[3]
        chonobo_bases = list() # stores the aligned bases in chonobo reads

        # get reads mapped to each difference using samfile.pielup(),
        # converting to 0-based coordinate
        # truncate = True: outputs the column in the exact region specified
        for pileupcolumn in samfile.pileup(chrm, p_1b -1, p_1b,
                                          truncate = True):
            
            for pileupread in pileupcolumn.pileups:
                # check for if aligned read's mapq
                if pileupread.alignment.mapq >= 30:
                    base = pileupread.alignment.seq[pileupread.qpos]
                    if pileupread.alignment.is_reverse:
                        base = base.lower() # reverse reads --> lower case
                    else:
                        base = base.upper() # forward reads --> upper case
                    base_q = pileupread.alignment.qual[pileupread.qpos]
                    # only save the base if base quality is higher than 20
                    if ord(base_q) - 33 >= 20:
                        chonobo_bases.append(base)

        #save the difference in an ChimpBonoDiff objec
        this_dif = ChimpBonoDiff(chrm, p_1b, chimp_base,
                                 bono_base, chonobo_bases)

        #if provided -c or --chonobo argument, prints out chonobo bases
        #otherwise, only output the dif location, ref chimp and ref bono bases
        if this_dif.meets_criteria():
            if args.chonobo:
                print this_dif
            else:
                print this_dif.str_no_chonobo_bases()

    f_dif.close()


class ChimpBonoDiff():
    """ this class stores base information of one difference"""

    def __init__(self, chrm, p_1b, chimp_base, bono_base, chonobo_bases):
        self.chrm = chrm # chromosome number
        self.p_1b = p_1b  # position with 1-based coordinate
        self.chimp_base = chimp_base
        self.bono_base = bono_base
        self.chonobo_bases = chonobo_bases # chonobo reads mapped to dif
        self.coverage = len(self.chonobo_bases)

    def __str__(self):
        """ return the contents of a ChimpBonoDiff object in a nice format"""
        return "%s\t%d\t%s\t%s\t%d\t%s\t" %(self.chrm, self.p_1b,
                                            self.chimp_base, self.bono_base,
                                            len(self.chonobo_bases),
                                            ''.join(self.chonobo_bases))

    def str_no_chonobo_bases(self):
        """ return the contents of a ChimpBonoDiff object in a nice formate
        without printing out the chonobo coverages nd chonobo bases"""
        return "%s\t%d\t%s\t%s\t" %(self.chrm, self.p_1b,
                                    self.chimp_base, self.bono_base,)

    def meets_criteria(self):
        """ this function checks if a dif object meets the criteria of
        selection, the criteria is:
        coverage >= 4 and coverage <= 23
        number of fowrad reads >= 2
        number of reverse reads >= 2
        forward reads contain at least one chimp base and one bonobo base
        reverse reads contain at least one chimp base and one bonobo base"""

        #seperate chonobo_bases into forward bases and reverse bases
        f_bases = list()
        r_bases = list()
        for x in self.chonobo_bases:
            if x.isupper():
                f_bases.append(x)
            else:
                r_bases.append(x)
        
        # remove the higher and lower 2.5% coverage
        if self.coverage < 4 or self.coverage > 23:
            return False
        elif len(f_bases) <2 or len(r_bases) < 2:
            return False
        elif self.chimp_base.upper() not in f_bases:
            return False
        elif self.bono_base.upper() not in f_bases:
            return False
        elif self.chimp_base.lower() not in r_bases:
            return False
        elif self.bono_base.lower() not in r_bases:
            return False
        else:
            return True

if __name__ == "__main__":
    sys.exit(main())
