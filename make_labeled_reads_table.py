#! /soe/mollyzhang/anaconda/bin/python
"""make_labeled_reads_table.py
this script will read a bam file and a difference file
and it will write the read pairs into a tab delimied text file 
in which the positions of chimp or bonobo bases mapped a read pair are listed

Created by Molly Zhang on 10/7/14 at 14:44
"""
import sys
import argparse
import pysam

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("make_labeled_reads_table.py")
    parser.add_argument("f_bam", help = "this is the bam file")
    parser.add_argument("f_dif", help = "this is the difference file")
    return parser.parse_args()

def main(args):
    """Main method"""
    args = parse_args()
    f_bam = pysam.Samfile(args.f_bam, "rb")
    # save all the differences in one table/dictionary
    difs_table = store_dif_data(args.f_dif)  
    last_read_pair = ["qname",["cp"],["bp"]]    
    
    for read in f_bam:
        # if this read name is different from last read name
        # write out the previous read and reset the last_read_pair to current read
        if last_read_pair[0] != read.qname:
            sys.stdout.write(last_read_pair[0] + "\t" + ",".join(last_read_pair[1]) + 
                            "\t" + ",".join(last_read_pair[2]) + "\n")
            last_read_pair[0] = read.qname
            last_read_pair[1] = []
            last_read_pair[2] = []
        
        # pass this read if its mapping quality is less than 30
        if read.mapq < 30:
            continue 

        # check if any read.positions is in the difference table
        # convert to 1-based coorndinate           
        for pos in read.positions:     
            position = f_bam.getrname(read.tid) + ":" + str(pos+1)             
            if position in difs_table:
                base = get_base(pos+1, read)[0]
                base_qual = get_base(pos+1, read)[1]
                # if base quality is less than 20, then skip this base
                if ord(base_qual) - 33 < 20 :
                    continue
                # if the base is the chimp base
                elif base == difs_table[position][0]:
                    last_read_pair[1].append(position)
    
                # if the base is bonobo base
                elif base ==  difs_table[position][1]:
                    last_read_pair[2].append(position)

    # write out the last read (pair)
    sys.stdout.write(last_read_pair[0] + "\t" + ",".join(last_read_pair[1]) +
                            "\t" + ",".join(last_read_pair[2]) + "\n")
                    
def get_base(p_ref_1b, read):
    """ return the index of base in a read aligned to reference position"""
    seq = list(read.seq)
    qual = list(read.qual)
    tracking = 0
    for section in read.cigar:
        if section[0] == 1 or section[0] == 4:
            # Insertion or soft-clip
            del seq[tracking: tracking + section[1]]
            del qual[tracking: tracking + section[1]]        
        elif section[0] == 2:
            # Deletion
            seq.insert(tracking, section[1] * 'x')
            qual.insert(tracking, section[1] * 'x')
            tracking += 1
        elif section[0] == 0:
            # match
            tracking += section[1]
        elif section[0] == 5:
            # hard-clip            
            pass
    seq = ''.join(seq)
    qual = ''.join(qual)
    return [seq[p_ref_1b -1 - read.pos], qual[p_ref_1b-1-read.pos]]



def store_dif_data(filename):
    """save save all the differences a dictionary
    key: chrm + positions (str), values: chimp base + bonobo base (concacnated)
    """
    difs = dict()
    f_dif = open(filename, "r")
    for line in f_dif:
        items = line.strip().split()
        difs[items[0] + ":" + items[1]] = items[2] + items[3]
    return difs        

if __name__ == "__main__" :
    sys.exit(main(sys.argv))
