#! /soe/mollyzhang/anaconda/bin/python
"""label_het_sites.py
this script label the hetertozygous sites in diana:

Created by Molly Zhang on 11/02/14 at 16:44
"""
import sys
import argparse
import pysam
import collections

def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("description")

    parser.add_argument("f_bam", help = "bam file of labeled reads")
    parser.add_argument("f_het", help = "the complete heterozygous sites")
    parser.add_argument("het_labeled", help = "output labeled sites")
    parser.add_argument("het_unlabeled", help = "output unlabeled sites")
    
    parser.add_argument("reads_wrong_label", 
                        help = "the wrongly labeled reads")

    return parser.parse_args()

def main(args):
    """Main method"""
    args = parse_args()

    f_bam = pysam.Samfile(args.f_bam, "rb")
    f_het = open(args.f_het, "r")
    het_labeled = open(args.het_labeled, "w")
    het_unlabeled = open(args.het_unlabeled, "w")
    f_reads_wrong_label = open(args.reads_wrong_label, "w")

    for line in f_het:
        items = line.strip().split("\t")
        chrm = items[0]
        p_1b = int(items[1])     # position is 1-based
        het_bases = items[2]
        
        pileup_bases = []
        pileup_labels = []
        reads_dumped = set()       # dump reads giving wrong info

        # print out unlabeled sites if there is no pileupcolumn 
        # (no reads mapped to the site)
        no_pileupcolumn = True
        for pileupcolumn in f_bam.pileup(chrm, p_1b-1, p_1b, truncate = True):
            no_pileupcolumn = False
            reads_not_dumped = []   # use reads giving right info
            for pileupread in pileupcolumn.pileups:
                
                read_label = pileupread.alignment.opt("cb")
                mapQ = pileupread.alignment.mapq
                cigar_parsed = get_base(p_1b, pileupread.alignment)
                base = cigar_parsed[0]
                baseQ = ord(cigar_parsed[1]) - 33 
                
                if baseQ < 20 or baseQ > 70 or mapQ < 30 or base not in het_bases:
                    reads_dumped.add(pileupread)
                else:
                    reads_not_dumped.append(pileupread)
                    pileup_bases.append(base.upper())
                    pileup_labels.append(read_label)

            # check for the status of the site: label, unlabel, conflict
            site_label = check_site_label(pileup_bases, pileup_labels, het_bases)
    
            if site_label[0]:
                new_line = (chrm + "\t" + str(p_1b) + "\t" +
                           site_label[1] + "\t" + site_label[2] + "\t" + 
                           items[3] + "\t" + items[4] + "\n")
                het_labeled.write(new_line)
            else:
                het_unlabeled.write(line)

        # prints out sites without any read mapped to it (no pileupcolumn)
        if no_pileupcolumn:
            het_unlabeled.write(line)

        # write out all dumped reads in one file
        for reads in reads_dumped:
            f_reads_wrong_label.write(reads.alignment.qname + "\n")


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


def check_site_label(bases, labels, het_bases):
    """ check whether a site should be labeled, unlabled, or conflict
    returns either of the three following lists:
    ["unlabeled",'','', label_table]
    ["labeled", chimp_base, bonobo_base, label_table]
    """
    is_labeled = True
    chimp_base = ''
    bono_base = ''

    # label table looks like {"c":"AAAAA", "b":"TTTTT"}
    label_table = dict()
    label_table1 = dict()

    for i in range(len(labels)):
        if labels[i] not in label_table:
            label_table[labels[i]] = bases[i]
        else:
            label_table[labels[i]] = label_table[labels[i]] + bases[i]

    if len(label_table) == 0:
        is_labeled = False

    else:
        # get the most common charater in each label
        for key in label_table:
            label_table1[key] = collections.Counter(label_table[key]).most_common(1)[0][0]


        if len(label_table) == 1 and is_labeled:
            if "c" in label_table:
                chimp_base = label_table1["c"]
                bono_base = het_bases.replace(chimp_base, '')      
            elif "b" in label_table:
                bono_base = label_table1["b"]
                chimp_base = het_bases.replace(bono_base, '')

        elif len(label_table) == 2 and is_labeled:
            chimp_base = label_table1["c"]
            bono_base = label_table1["b"]

    return [is_labeled, chimp_base, bono_base] 


if __name__ == "__main__" :
    sys.exit(main(sys.argv))
