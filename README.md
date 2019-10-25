# TEPID-run
The Scripts were run using IBM LSF HPC cluster. Needs to be adapted according to the type of cluster.

### 01. Scripts to run TEPID
**Create Yaha and Bowtie index**  
yaha -g refgenome.fasta -L 11 -H 2000  
bowtie -f refgenome.fasta --threads 5 refgenome

**Step1:** Run TEPID_map using submit.tepid-map.lsf  
Output: Bam files and split and unmapped reads files

**Step2:** Run TEPID_discover, separately for insertion and deletion using submit.tepid-discover_del.lsf and submit.tepid-discover_ins.lsf  
Output: bed files for insertion and deletion

Run the above for all samples in a population, keep separate folders for each sample

**Step3:** Run merge_insertions.py and merge_deletions.py from the main folder with folders of each sample
merge_insertions.py -f insertions
merge_deletions.py -f deletions

**Step4:** Run TEPID_refine, separately for insertion and deletion (don't flip deletions) using	submit.tepid-refine_ins.lsf and submit.tepid-refine_noflipdel.lsf

**Step5:** Run genotype.py, separately for insertion and deletion with the arguments:
-d run on deletions
-i run on insertionss
-a ambiguous variants filename
-r name of the reference accession


flip the genotyped deletion file and concatenate with genotyped insertion file

### 02. Reformat polyTEs
Reformats the genotyping output from genotype.py in TEPID as a matrix with some classification of the loci (which I did not use)

To run the reformat_tepav_matrix_del.py first flip the genotyped deletion calls, making called samples as absence (0) and no calls as presence (1).
To do this just exchange the last two columns of genotyped_del_with_second_pass.bed
This is necessary to merge the insertion and deletion calls of TEPID into one dataset.
