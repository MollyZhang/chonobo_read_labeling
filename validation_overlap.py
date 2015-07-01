#! /soe/mollyzhang/anaconda/bin/python
"""validation_overlap.py
edited by Molly Zhang on 06/22/15
this script will calculate how well the newly labeled het sites agree with the old ones

Here is a more specific example to describe what it does: 

file1 content
id  label
1   none
2   none 
3   none   
4   none  
5   none
6   none
7   none     
8   none
9   none
10  none

file 2 content:
id  label
2   x
3   x
5   x
8   x
10  x

file 2a --> get a subset of rows from file 2:
id  label
2   x
5   x
10  x

file 2b --> the rest of rows from file 2 not included in file3a:
id  label
3   x
8   x

file 3 --> use file 3a to generate a set of label for file 1:
id  label
1   x
2   x1
3   x
4   none
5   x1
6   none
7   x
8   x
9   x
10  x


what this script does are:
1) compare file3 with file2a to make sure the label agree with each other
2) compare file3 with file2b to find out the labeling accuracy
"""

import sys
import argparse


def parse_args():
    """Set up and parse command-line arguments."""
    parser = argparse.ArgumentParser("validation_overlap.py")
    parser.add_argument("file2a", help="read doc string of script for info")
    parser.add_argument("file2b", help="read doc string of script for info")
    parser.add_argument("file3", help="read doc string of script for info")
    parser.add_argument("-v", "--verbose", help="print out all mislabeling",
                        action="store_true")
    return parser.parse_args()

def main(args):
    args = parse_args()
    file2a = save_het_to_dict(args.file2a)
    file2b = save_het_to_dict(args.file2b)
    file3 = save_het_to_dict(args.file3)
    compare_labels_validate(file2a, file3)
    compare_labels_prediction(file2b, file3)

def compare_labels_validate(origin_file, new_file):
    n_p = 0     # number of overlaps between positions
    n_pl = 0    # number of overlaps between positiosn AND labels
    f_output1 = open("validation_original_label_nonexist", "w")
    f_output1.write("original label not included in new label:\n")
    f_output2 = open("validation_new_label_conflict", "w")
    f_output2.write("new label in conflict with old label:\n")
    f_output2.write("position\told_label\tnew_label\n")
    for key in origin_file:
        if key in new_file:
            n_p += 1
            if origin_file[key] == new_file[key]:
                n_pl += 1
            else:
                f_output2.write(key + "\t"+ origin_file[key] + "\t" + 
                                new_file[key] + "\n")
        else:
            f_output1.write(key + "\n")

    print "number of common sites/original sites: %d/%d = %f" % \
          (n_p, len(origin_file), n_p/float(len(origin_file)))
    print "number of correct labeling/all labeling: %d/%d = %f" % \
          (n_pl, n_p, n_pl/float(n_p))



def compare_labels_prediction(origin_file, new_file):
    n_p = 0     # number of overlaps between positions
    n_pl = 0    # number of overlaps between positiosn AND labels
    f_output = open("prediction_incorrect_label", "w")
    f_output.write("position\told_label\tnew_label\n")
    for key in origin_file:
        if key in new_file:
            n_p += 1
            if origin_file[key] == new_file[key]:
                n_pl += 1
            else:
                f_output.write(key + "\t"+ origin_file[key] + "\t" + 
                                new_file[key] + "\n")
    
    print "number of common sites/original sites: %d/%d = %f" % \
          (n_p, len(origin_file), n_p/float(len(origin_file)))
    print "number of correct labeling/all labeling: %d/%d = %f" % \
          (n_pl, n_p, n_pl/float(n_p))

def save_het_to_dict(filename):
    """ save a labeled het sites file into dictionary
    key: chrm:position
    value: chimpbase_bonobobase
    example:
    het["chr1:969032"] = "C_T"
    """
    het = dict()
    f = open(filename, "r")
    for line in f:
        items = line.strip().split("\t")
        key = items[0] + ":" + str(items[1])
        value = items[2] + "_" + items[3]
        het[key] = value
    f.close()
    return het

if __name__ == "__main__" :
    sys.exit(main(sys.argv))

