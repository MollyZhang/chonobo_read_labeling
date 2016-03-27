### chonobo_read_labeling -- A pipeline to phase the genome of a chimpanzee-bonobo hybrid


####Step1. label the reads with SNPs they map to based on heterozygous sites

`make_labeled_reads_table.py bamfile_input heterozygous_site_input > step1_output`

####Step2. add a forth column in the labeling table as either "b", "c" or "conflict"

`read_labeling_and_stats.py step1_output step2_output`

####Step3. save all labeled reads into a bam file

`save_labeled_reads_to_bam.py bamfile_input step2_output step3_output.bam`

####Step4. use samtools to sort and index the labeled reads bam file

`samtools sort step3_output.bam step3_output_r1.s`

`samtools index output_labeled.bam`

####step5. label all heterozygous sites with labeled reads

`label_het_sites.py output_labeled.bam all_het_sites_input output1_labeled_sites output2_unlabeled_sites output3_wrongly_labeled_reads`
