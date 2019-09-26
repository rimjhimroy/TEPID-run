Reformats the genotyping output from genotype.py in TEPID as a matrix with some classification of the loci (which I did not use)

To run the reformat_tepav_matrix_del.py first flip the genotyped deletion calls, making called samples as absence (0) and no calls as presence (1).
To do this just exchange the last two columns of genotyped_del_with_second_pass.bed
This is necessary to merge the insertion and deletion calls of TEPID into one dataset.