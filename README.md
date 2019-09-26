# TEPID-run
The Scripts were run using IBM LSF HPC cluster. Needs to be adapted according to the type of cluster.

### Scripts to run TEPID
**Create Yaha and Bowtie index**  
yaha -g refgenome.fasta -L 11 -H 2000
bowtie -f refgenome.fasta --threads 5 refgenome

**Step1:** Run TEPID_map.lsf
Output: Bam files and split and unmapped reads files

**Step2:** Run TEPID_discover.lsf, separately for insertion and deletion
Output: bed files for insertion and deletion

Run the above for all samples in a population, keep separate folders for each sample

**Step3:** Run merge_insertions.py and merge_deletions.py from the main folder with folders of each sample
merge_insertions.py -f insertions
merge_deletions.py -f deletions

**Step4:** Run TEPID_refine.lsf, separately for insertion and deletion (don't flip deletions)

**Step5:** Run genotype.py, separately for insertion and deletion with the arguments:
-d run on deletions
-i run on insertionss
-a ambiguous variants filename
-r name of the reference accession


flip the genotyped deletion file and concatenate with genotyped insertion file
