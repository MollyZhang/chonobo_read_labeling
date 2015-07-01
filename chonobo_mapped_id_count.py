#! /usr/bin/python

import collections

#read files content to list
with open("chonobo_100k.bono.s.mapped_ids.txt", 'r') as f:
    bono = [line.rstrip('\n') for line in f]
with open("chonobo_100k.pt4.s.mapped_ids.txt", 'r') as f:
    chimp = [line.rstrip('\n') for line in f]
with open("chonobo_100k.hg19.s.mapped_ids.txt", 'r') as f:
    human = [line.rstrip('\n') for line in f]
     

#check for possible repeats in each list
def check_for_repeat(list_name, string): 
    	repeat_list = [x for x, y in collections.Counter(list_name).items() if y > 1]
	print "repeats in file %s: %d" %(string, len(repeat_list))

check_for_repeat(bono, "bonobo")
check_for_repeat(chimp, "chimpanzee")
check_for_repeat(human, "human")


#count numbers of mapped reads in a list
print "the number of 100k chonobo mapped to bonobo: %d" %len(bono)
print "the number of 100k chonobo mapped to chimpanzee: %d" %len(chimp)
print "the number of 100k chonobo mapped to human: %d" %len(human)


#print the number of intersections between lists:
overlap_bono_chimp = set(bono).intersection(set(chimp))
overlap_bono_human = set(bono).intersection(set(human))
overlap_chimp_human = set(chimp).intersection(set(human))
overlap_bono_chimp_human = set(bono).intersection(set(chimp)).intersection(set(human))

print "the number of overlap between bonobo and chimp: ", len(overlap_bono_chimp)
print "the number of overlap between bonobo and human: ", len(overlap_bono_human)
print "the number of overlap between chimp and human: ", len(overlap_chimp_human)
print "the number of overlap between bonobo, chimp and human", len(overlap_bono_chimp_human)


