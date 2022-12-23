################################################################################################################
###    Script to extract basecalled fast5, fastq and alignments from full length reads from two samples to   ###
### create a sample with known stoichiometry and a random population of modified reads across stoichiometries ###
################################################################################################################


##Read inputs from the command line:
modified=$1
unmodified=$2

reads_mod=$3
reads_unmod=$4

output_folder=$5
num_readsModified=$6
num_readsUnmodified=$7

#Generate fast5 folder:
echo 'Fast5 processing'
mkdir $output_folder
cd ./$output_folder
mkdir fast5_files
mkdir fastq_files
mkdir alignment
mkdir QC_files
cd fast5_files

echo 'Generating readID list - Modified'
sort -R ./../../$6 | head -n $num_readsModified > list_reads_mod.txt

echo 'Subsetting fast5'
fast5_subset -i ./../../$modified/fast5_files -s $PWD -l list_reads_mod.txt -f batch_modified -n 4000 --recursive 
mv filename_mapping.txt filename_mapping_modified.txt

echo 'Generating readID list - Unmodified'
sort -R ./../../$7 | head -n $num_readsUnmodified > list_reads_unmod.txt

echo 'Subsetting fast5'
fast5_subset -i ./../../$unmodified/fast5_files -s $PWD -l list_reads_unmod.txt -f batch_unmodified -n 4000 --recursive
mv filename_mapping.txt filename_mapping_unmodified.txt


#Generate final_summary.stats file:
echo 'Generating QC_files'
cat filename_mapping_modified.txt filename_mapping_unmodified.txt > final_summary.stats
rm filename_mapping_*
sed -i '1 i\read_id\tfilename' final_summary.stats
mv final_summary.stats ./../QC_files

#Generate fastq files: 
cd ./../fastq_files
echo 'Extract modified and unmodified reads from fastq'
seqkit grep -f ./../fast5_files/list_reads_mod.txt ./../../$modified/fastq_files/*.fq.gz > Modified.fq
seqkit grep -f ./../fast5_files/list_reads_unmod.txt ./../../$unmodified/fastq_files/*.fq.gz > Unmodified.fq
cat Modified.fq Unmodified.fq > $output_folder.fq
gzip $output_folder.fq
rm *.fq

#Generate bam files:
cd ./../alignment
echo 'Generating bam files'
samtools view -H ./../../$unmodified/alignment/*.bam > header.sam
samtools view ./../../$modified/alignment/*.bam | fgrep -w -f ./../fast5_files/list_reads_mod.txt > modified.sam 
samtools view ./../../$unmodified/alignment/*.bam | fgrep -w -f ./../fast5_files/list_reads_unmod.txt > unmodified.sam 
cat header.sam modified.sam unmodified.sam > concat.sam
samtools view -Sb concat.sam > concat.bam
samtools sort -o $output_folder.sorted.bam concat.bam
samtools index $output_folder.sorted.bam
rm concat.bam
rm *.sam
