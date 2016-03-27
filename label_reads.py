#! /soe/mollyzhang/anaconda/bin/python
"""add_tag.py
this script will add tags to reads in a samfile
specifying the number of detected chimp bases and bonobo bases

Created by Molly Zhang on 08/22/14 at 14:44
"""
import sys
import argparse
import pysam

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("tag_reads.py")
    parser.add_argument("f_bam", help = "this is the bam file")
    parser.add_argument("f_dif", help = "this is the difference file")
    parser.add_argument("bam_out", help = "the output bam file")
    parser.add_argument("stat_out", help = "the output stat file")
    return parser.parse_args()

def main(args):
    """Main method"""
    args = parse_args()
    f_bam = pysam.Samfile(args.f_bam, "rb")
    f_labeled = pysam.Samfile(args.bam_out, "wb", template=f_bam)
    f_stat = open(args.stat_out, "w")    
    # save all the differences in one table/dictionary
    difs_table = store_dif_data(args.f_dif)
    
    num_read_labeled = 0
    num_read_labeled_chimp = 0
    num_read_labeled_bonobo = 0
    num_read_labeled_both = 0
    for read in f_bam:
        read_labeled = False
        read_labeled_chimp = False
        read_labeled_bonobo = False
        # add empty tag of cp (chimp positions) and bp (bonobo positions)
        # only if there isn't any cp or bp tag already     
        try: read.opt("cp")        
        except: read.setTag("cp", "")
        
        try: read.opt("bp")
        except: read.setTag("bp", "")
        # pass this read if its mapping quality is less than 30
        if read.mapq < 30:
            continue
        # read.positions: the list of reference positions read aligns to
        for pos in read.positions:
            # check if any read.positions is in the difference table
            # convert to 1-based coorndinate        
            position = f_bam.getrname(read.tid) + str(pos+1)             
            if position in difs_table:
                base = get_base(pos+1, read)[0]
                base_qual = get_base(pos+1, read)[1]
                # if quality_filter is set true and base quality is less
                # than 20, then skip this base
                if ord(base_qual) - 33 < 20 :
                    continue
                # if the base is the chimp base
                elif base == difs_table[position][0]:
                    read_labeled = True
                    read_labeled_chimp = True
                    if read.opt("cp") != "":
                        read.setTag("cp", read.opt("cp") + "," + str(pos+1) )
                    else:
                        read.setTag("cp", str(pos+1))
                # if the base is bonobo base
                elif base ==  difs_table[position][1]:
                    read_labeled = True
                    read_labeled_bonobo = True
                    if read.opt("bp") != "":
                        read.setTag("bp", read.opt("bp") + "," + str(pos+1) )
                    else:
                        read.setTag("bp", str(pos+1))
        f_labeled.write(read)
            
        if read_labeled:
            num_read_labeled += 1
        if read_labeled_chimp:
            num_read_labeled_chimp += 1
        if read_labeled_bonobo:
            num_read_labeled_bonobo += 1
        if read_labeled_chimp and read_labeled_bonobo:
            num_read_labeled_both += 1
    
    f_stat.write("total number of mapped reads in bam file: %d\n"
                %f_bam.mapped )
    f_stat.write("total number of difference positions: %d\n" %len(difs_table))
    f_stat.write("total number of labled reads: %d, percentage: %.3f"
                %(num_read_labeled, num_read_labeled/float(f_bam.mapped)*100)
                + " %\n"
                )
                
    f_stat.write("total number of reads labeled as chimpanzee: %d\n"
                %num_read_labeled_chimp)
    f_stat.write("total number of reads labeled as bonobo: %d\n"
                %num_read_labeled_bonobo)
    f_stat.write("reads labeled as both chimp & bonobo: %d, percentage: %.3f"
                % (num_read_labeled_both,
                   num_read_labeled_both/float(num_read_labeled)*100)
                + " %\n"
                )

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
        difs[items[0] + items[1]] = items[2] + items[3]
    return difs        

if __name_ == "__main__" :
    sys.exit(main(sys.argv))
