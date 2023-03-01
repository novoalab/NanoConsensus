############################################################################################################
### Script to extract basecalled fast5, fastq and alignments from full length reads from a single sample ###
############################################################################################################


##Read inputs from the command line:
modified=$1
output_folder=$2
list_reads=$3

#Generate fast5 folder:
echo 'Fast5 processing'
mkdir $output_folder
cd ./$output_folder
mkdir fast5_files
mkdir fastq_files
mkdir alignment
mkdir QC_files
cd fast5_files

echo 'Subsetting fast5'
fast5_subset -i ./../../$modified/fast5_files -s $PWD -l $list_reads -n 4000 --recursive 

#Generate final_summary.stats file:
echo 'Generating QC_files'
mv filename_mapping.txt final_summary.stats
sed -i '1 i\read_id\tfilename' final_summary.stats
mv final_summary.stats ./../QC_files

#Generate fastq files: 
cd ./../fastq_files
echo 'Extract modified and unmodified reads from fastq'
seqkit grep -f $list_reads ./../../$modified/fastq_files/*.fq.gz > $output_folder.fq
gzip $output_folder.fq

#Generate bam files:
cd ./../alignment
echo 'Generating bam files'
samtools view -H ./../../$modified/alignment/*.bam > header.sam
samtools view ./../../$modified/alignment/*.bam | fgrep -w -f $list_reads > modified.sam 
cat header.sam modified.sam > concat.sam
samtools view -Sb concat.sam > concat.bam
samtools sort -o $output_folder.sorted.bam concat.bam
samtools index $output_folder.sorted.bam
rm concat.bam
rm *.sam
