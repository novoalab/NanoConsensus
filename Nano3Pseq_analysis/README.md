## Analysis of Nano3P-seq datasets
### Basecalling and demultiplexing:
The module [mop_preprocess](https://biocorecrg.github.io/master_of_pores/nanopreprocess.html) module from [Master of Pores2](https://github.com/biocorecrg/MOP2) pipeline was used to perform basecalling and demultiplexing with guppy (version 4.0, model dna_r9.4.1_70bps_hac). To execute this module, please add these parameters into both `params.config` and `tools_opt.tsv` files and run the command below: 

```
nextflow run mop_preprocess.nf -with-singularity -bg > log.txt
```

Before proceding to mapping, all reads whose barcode could not be identified by guppy (ie: classified as unknown), went through a second round of demultiplexing using readucks (version 0.0.3) using the code below:

```
readucks -i /path/to/unknown_bin -o /path/to/output --limit_barcodes_to 1 2 3 4 5 6 --adaptter_threshold 73 --threshold 50
```

Finally, fastq files from both runs of demultiplexing for every barcode were merged and compressed:
```
#Merging fastq:
zcat Guppy_BC1.fq.gz Readucks_BC1.fq.gz > BC1.fq

#Compressing the merged data:
zcat BC1.gz
```

### Mapping:
**To the bacterial genome:**

```
#Mapping spliced:
minimap2 -ax splice -k14 Ecoli_BW25113_NCBI.fa WT_1h_Rep1.fq.gz | samtools view -F 4 -Sb > WT_1h_Rep1.genome.bam
```

**To the rRNA reference:**

```
#Mapping unspliced:
minimap2 -ax map-ont Ecoli_BW25113_NCBI.fa WT_1h_Rep1.fq.gz | samtools view -F 4 -Sb > WT_1h_Rep1.rRNA.bam
```

**To the transcriptome:**

Before mapping, the transcriptome reference has to be generated using the genome (.fa) and annotation (.gtf) files and the code below:
```
#Filter the gtf file to only include gene features:
MISSING: To fill code to generate only genes gtf. 

#Extract the reference transcriptome:
bedtools getfasta -fi Ecoli_BW25113_NCBI.fa -bed Ecoli_BW25113_NCBI.OnlyGenes.gtf -s -name > Ecoli_BW25113_NCBI.Transcriptome.fa

#Mapping unspliced:
minimap2 -ax map-ont Ecoli_BW25113_NCBI.Transcriptome.fa WT_1h_Rep1.fq.gz | samtools view -F 4 -Sb > WT_1h_Rep1.transcriptome.bam
```

### Counting:
In this study, three different set of counts were generated:

**1. Counts based on genome alignments with htseq-count (version 0.13.5) with default parameters for cDNA:**

**2. Counts based on genome alignments with htseq-count (version 0.13.5) with additional parameters to take into account the bacterial operons:**

**3. Counts based on transcriptome alignments with salmon (version 1.9.0) with default parameters for ONT datasets:**

```
salmon quant --ont -t ./ref/Ecoli_BW25113_NCBI.Transcriptome.fa -l SR -a WT_1h_Rep1.transcriptome.bam -o WT_1h_Rep1
```

### Differential expression analysis and visualization:
