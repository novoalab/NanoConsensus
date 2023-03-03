## Analysis of the DRS datasets
### Basecalling and demultiplexing:
The module [mop_preprocess](https://biocorecrg.github.io/master_of_pores/nanopreprocess.html) module from [Master of Pores2](https://github.com/biocorecrg/MOP2) pipeline was used to perform basecalling, with guppy (version 3.1.5, model rna_r9.4.1_70bps_hac), and demultiplexing - when needed - with Deeplexicon (version 1.2.0). Then, it also mapped the basecalled reads to the suitable rRNA transcriptome (for E.coli: `ref/Escherichia_coli.rRNA.fa` and for S.cerevisiae: `ref/Saccharomyces_cerevisiae.rRNA.fa`) with graphmap (version 0.5.2) with unspliced parameters (-v1 -K fastq). To execute this module, please add these parameters into both `params.config` and `tools_opt.tsv` files and run the command below: 

```
nextflow run mop_preprocess.nf -with-singularity -bg > log.txt
```

### Extracting full-length reads:
For each sample and transcript, only full-length reads were included in the different analysis done in this study:
* Generate bed file from bam file to know the alignment start and end for each read:
```
bedtools bamtobed -i Ecoli.sorted.bam > Ecoli.sorted.bam.bed
```
* Extracting reads' ID from full length reads:
```
#Usage: 
bash ./scripts/Full_Length/Extract_FullLength_IDs.sh /path/to/bed chr start_position end_position output_name

#Examples:
bash ./scripts/Full_Length/Extract_FullLength_IDs.sh Ecoli.sorted.bam.bed 16S 50 1525 Ecoli
bash ./scripts/Full_Length/Extract_FullLength_IDs.sh Ecoli.sorted.bam.bed 23S 50 2894 Ecoli
```
* Extracting basecalled fast5, fastq and alignments from full-length reads:

  a) **Downsampling with a fixed population of modified reads for benchmarking**: Script to merge reads from modified and unmodified samples to build datasets with known stoichiometry and the same population of modified reads across stoichiometries:
  ```
  #Usage: 
  bash ./scripts/Downsampling/Fixed_Modified_Population/FixedDownsampling_Basecalled_Fastq_Bam.sh /path/modified/sample /path/unmmodified/sample number_modified_reads number_unmodified reads output_folder full-length_reads_modified.txt full-length_reads_unmodified.txt

  #Examples:
  bash ./scripts/Downsampling/Fixed_Modified_Population/FixedDownsampling_Basecalled_Fastq_Bam.sh Parent rsmA_KO 100 100 16S_Parent-50_rsmAKO-50_Rep1 Parent.16S.read_IDs.txt rsmA-KO.16S.read_IDs.txt
  ```
  
  b) **Random downsampling for benchmarking**: Script to merge reads from modified and unmodified samples to build datasets with known stoichiometry and a random population of modified reads across stoichiometries:
   ```
  #Usage: 
  bash ./scripts/Downsampling/Random/RandomDownsampling_Basecalled_Fastq_Bam.sh /path/modified/sample /path/unmmodified/sample number_modified_reads number_unmodified reads output_folder full-length_reads_modified.txt full-length_reads_unmodified.txt

  #Examples:
  bash ./scripts/Downsampling/Random/RandomDownsampling_Basecalled_Fastq_Bam.sh Parent rsmA_KO 100 100 16S_Parent-50_rsmAKO-50_Rep1 Parent.16S.read_IDs.txt rsmA-KO.16S.read_IDs.txt
  ```
  
  c) **For antibiotic exposure analysis**: Script to extract basecalled fast5, fastq and alignments from full-length reads from single samples: 
  ```
  #Usage: 
  bash ./scripts/Full_Length/FullLength_Basecalled_Fastq_Bam.sh /path/input/sample output_folder full-length_readsIDs.txt

  #Examples:
  bash ./scripts/Full_Length/FullLength_Basecalled_Fastq_Bam.sh Ecoli_Untreated_1h Ecoli_Untreated_1h_FullLength Ecoli_Untreated_1h.16S.read_IDs.txt 
  ```

### Running RNA modification detection softwares: 

EpiNano, Nanopolish, Tombo and Nanocompore were run by the module [mop_mod](https://biocorecrg.github.io/master_of_pores/nanomod.html) module from [Master of Pores2](https://github.com/biocorecrg/MOP2) with the default options:

```
nextflow run mop_mod.nf -with-singularity -bg > log.txt
```

### Running NanoConsensus:

NanoConsensus was run by the module mop_consensus from [Master of Pores2](https://github.com/biocorecrg/MOP2) with the options `--MZS 3.75` and `--NC 4`.

```
nextflow run mop_consensus.nf -with-singularity -bg > log.txt
```
