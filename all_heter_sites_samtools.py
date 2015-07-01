#! /soe/mollyzhang/anaconda/bin/python
"""all_heter_sites.py

take samtools mpileup result and select for the heterozygous sites using 
the criteria: coverage >= 4 and coverage <=23
at least one forward base and reverse base match one base
at least one forward base and reverse base match anothter base
if only one read has indel, ignore this read
if more than read has indel, ignore the position

Created by Molly Zhang on 09/17/14 at 15:59
"""
import sys
import argparse

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("all_heter_sites.py")
    parser.add_argument("mpileup", 
                        help = "the input file of samtools mpileup")
    parser.add_argument("out", 
                        help = "file name of the output")
    return parser.parse_args()

def main(args):
    """Main method"""
    args = parse_args()
    f_het = open(args.out, "w")
    f_pileup = open(args.mpileup, "r")
    unprocessed_cases = 0
    for line in f_pileup:    
        items = line.strip().split()
        chrm = items[0]
        position = items[1]
        coverage = int(items[2])
        if coverage < 4 or coverage > 23:
            continue
        else:
            bases = items[3]
            if not has_multiple_bases(bases):
                continue
            elif bases.count('*') > 1:
                continue
            elif has_other_symbol(extract(bases)):
                unprocessed_cases += 1
                continue
            else:
                bases = extract(bases)
                if is_heter(bases)[0]:
                    f_het.write(chrm + "\t" + position + "\t" + 
                    ''.join(is_heter(bases)[1]) + "\t"
                    + str(coverage) + "\t" + items[3] + "\n")
    print "number of unprocessed cases: %d" %unprocessed_cases

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


def has_other_symbol(string):
    """ check wether a base string has any other uncommon symbols"""
    symbols = list("agtcAGTC^$+-*0123456789")
    for x in string:
        if x not in symbols:
            return True
    return False


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
    
def extract(string):
    """ handle start of a line: ^* , end of a line: $ and indels"""
    bases = list(string)
    # handle situations with ^:    
    # enumerate all bases and save the indices i+1 (immediately after "^")  
    indices = [i+1 for i, j in enumerate(bases) if j == '^']
    if len(indices) > 0:
        # save all items in bases excet for the ones in indices in a new list 
        bases = [i for j, i in enumerate(bases) if j not in indices]
    
    # handle situations with + and -:
    while "+" in bases:
        bases = handle_indel("+", bases)
    while "-" in bases:        
        bases = handle_indel("-", bases)     
        
    return "".join(bases)

def handle_indel(indel, bases):
    """ delete a number of bases after + or - sign"""
    n = [] # the number of bases to delete
    i = 1 # the number of digits that has
    index = bases.index(indel)
    while bases[index+i].isdigit():
        n.append(bases[index+i])
        i += 1            
    n = int("".join(n))
    del bases[index: (index+i+n)]
    return bases

if __name__ == "__main__" :
    sys.exit(main(sys.argv))

