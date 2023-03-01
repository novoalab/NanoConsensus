## Analysis of Nano3P-seq datasets:

### Basecalling and trimming:

### Mapping to rRNA reference: 

### Mapping to the transcriptome: 
* Generate the transcriptome reference from genome (.fa) and annotation (.gtf) files:
```
MISSING: To fill code to generate only genes gtf. 
bedtools getfasta -fi Ecoli_BW25113_NCBI.fa -bed Ecoli_BW25113_NCBI.OnlyGenes.gtf -s -name > Ecoli_BW25113_NCBI.Transcriptome.fa
```

* Mapping and counting with salmon:
```
minimap2 -ax map-ont Ecoli_BW25113_NCBI.Transcriptome.fa WT_1h_Rep1.fq.gz | samtools view -F 4 -Sb > WT_1h_Rep1.transcriptome.bam
salmon quant --ont -t /no_backup_isis/enovoa/reference_fasta/bacterial_references/Ecoli_BW25113_NCBI.Transcriptome.fa -l SR -a WT_1h_Rep1.transcriptome.bam -o WT_1h_Rep1
```
