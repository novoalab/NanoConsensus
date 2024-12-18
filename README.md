

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5805806.svg)](https://doi.org/10.5281/zenodo.5805806)

# NanoConsensus: consensus prediction of RNA modifications from direct RNA nanopore sequencing data

*NanoConsensus* is an software to robustly identify putative RNA modified sites from direct RNA sequencing datasets. 

## Table of Contents  
- [General Description](#General-Description)
- [Installation](#Installation)
- [Running the code](#Running-the-code)
- [Expected output](#Expected-output)
- [Required dependencies](#Required-dependencies)
- [Citation](#Citation) 
- [Contact](#Contact) 

## General Description

NanoConsensus performs pairwise comparisons between two conditions (e.g. WT vs KO) using four different RNA modification detection softwares ([Epinano](https://github.com/enovoa/EpiNano), [Nanopolish](https://github.com/jts/nanopolish), [Tombo](https://github.com/nanoporetech/tombo) and [Nanocompore](https://github.com/tleonardi/nanocompore)) at per-transcript level. Then, it combines all results generated to produce a consensus prediction, which is more robust than those from individual softwares.  

![NanoConsensus_scheme](/img/NanoConsensus_scheme.png)


*NanoConsensus* performs the following steps: 
* **1. Running RNA modification detection softwares.** The first step of *NanoConsensus*  runs the four different RNA modification detection algorithms (EpiNano, Tombo, Nanopolish and Nanocompore) in a pairwise manner (i.e. comparing a WT/control/condition1 against an IVT/Knockout/Condition2 sample. This step is implemented in the form of Nextflow pipeline, and has been embedded [Master of Pores](https://github.com/biocorecrg/master_of_pores), in the form of an updated [NanoMod](https://biocorecrg.github.io/master_of_pores/nanomod.html) module. These algorithms will produce scores (p-values, differential error or differential current intensity) at per-transcript level for each individual position. 

* **2. Per-transcript Z-Score normalization.** *NanoConsensus* then performs Z-Score normalization for each dataset and method (EpiNano, Tombo, Nanopolish and Nanocompare), in a per-transcript manner.  Z-scores are then used to select candidate RNA modified positions for each independent software, which must be higher than the provided user-defined threshold (default value = 5). 

* **3. Flexible overlapping.** In the following step, candidate positions for each individual software identified in step 2 are extended into 5-mers. Flexible overlapping is then performed to identify overlapping k-mers across softares. The regions supported by two or more softwares are then saved as putative modified sites. 

* **4. Re-scaling of Z-scores.** Z-score values are affected by coverage of the transcript. Thus, to have final merged results that are comparable across transcripts, *NanoConsensus* will rescaled Z-score values between 0 and 1. Then, it calculates the *NanoConsensus score*, which is equal to the median of the rescaled Z-scores across all softwares. This is performed at all positions across the determined transript. 

* **5. Filtering of putative modified sites.** *NanoConsensus scores* from all previously identified putative modified sites are retrieved and compared to a specific threshold. This threshold is determined by the median of the *Nanoconsensus scores* across the entire transcript multiplied by an integer, with 5 as the default value. Those putative modified sites whose *NanoConsensus scores* are above the threshold are then reported whilst all unverified results are discarded.

## Installation 
If needed, install *Master of Pores* from its github repository:
```
git clone https://github.com/biocorecrg/master_of_pores_dev.git
```

Then, install *NanoConsensus* from its github repository:
```
git clone https://github.com/ADelgadoT/NanoConsensus.git
```

## Running the code 
Firstly, user should basecall and map direct RNA sequencing datasets using the [NanoPreprocess](https://biocorecrg.github.io/master_of_pores/nanopreprocess.html) module from [Master of Pores](https://github.com/biocorecrg/master_of_pores) pipeline. To run this module, please fill in `params.config` file with the required information and then, use the command below. If more information is needed, please click [here](https://biocorecrg.github.io/master_of_pores/nanopreprocess.html). 
```
nextflow run nanopreprocess.nf -with-singularity -bg > log.txt
```

The following step is to run the four different RNA modification detection algorithms (EpiNano(v1.1), Tombo(v1.5), Nanopolish(v0.11.1) and Nanocompore(v1.0.0)) in a pairwise manner with [NanoMod](https://biocorecrg.github.io/master_of_pores/nanomod.html) module from [Master of Pores](https://github.com/biocorecrg/master_of_pores) pipeline. Fill in both `params.config` and `comparison.tsv` file and then use the command below. If more information is needed, please click [here](https://biocorecrg.github.io/master_of_pores/nanomod.html).

```
nextflow run nanomod.nf -with-singularity -bg > log.txt
```

Finally, to perform NanoConsensus analysis, together with the files generated by *NanoMod*, the user must provide the **following inputs**:
* Reference file (*.fa) ([Check full input](example_input/Reference.fa))
* Transcript name - should match a fasta ID in the reference file - i.e: *18S*
* Start position - must be an integer - i.e: *50*
* End position - must be an integer - i.e: *1492*

These inputs - generated by *NanoMod* - are also required for the analysis:
* **Epinano output**: both for WT and IVT samples ([Check full input](example_input/Epinano/epinano_WT.tsv.per.site.var.csv))
```
#Ref,pos,base,cov,q_mean,q_median,q_std,mis,ins,del
23S,2515,C,45497.0,5.36995,4.00000,3.97797,0.0822032221904741,0.18715519704595907,0.2058377475437941
23S,2516,A,45504.0,5.38207,4.00000,4.71619,0.17128164556962025,0.20497099156118143,0.07733386075949367
23S,2517,C,45529.0,6.92130,5.00000,5.04250,0.06165301236574491,0.1505633771881658,0.13540820136616222
23S,2518,A,45545.0,6.49821,5.00000,5.47485,0.10802503018992206,0.10855198155670216,0.2082775277198375
23S,2519,T,45557.0,6.51247,5.00000,4.81853,0.09386043857145993,0.14792457800118533,0.2033057488421099
23S,2520,C,45609.0,10.72107,9.00000,7.03936,0.03990440483237957,0.1465938740160933,0.020105680896314326
```

* **Nanopolish output**: both for WT and IVT samples ([Check full input](example_input/Nanopolish/Nanopolish_WT.tsv.per.site.var.csv))
```
contig	position	reference_kmer	read_name	median	coverage
23S	2515	ACATC	1	79.03	271377
23S	2516	CATCC	1	87.0	118486
23S	2517	ATCCT	1	74.47	124448
23S	2518	TCCTG	1	73.21	65149
23S	2519	CCTGG	1	88.2	84488
23S	2520	CTGGA	1	107.84	125206
```

* **Tombo output** ([Check full input](example_input/Tombo/WT_IVT_plus_Tombo_Output.tsv))
```
"Ref_Position"	"Chr"	"Position"	"Tombo_SiteScore"	"Coverage_Sample"	"Coverage_IVT"	"Tombo_KmerScore"
"23S_2515"	"23S"	"2515"	"0.0000"	"37042"	"37884"	0
"23S_2516"	"23S"	"2516"	"0.0000"	"37039"	"37880"	5e-04
"23S_2517"	"23S"	"2517"	"0.0000"	"37053"	"37895"	0.0012
"23S_2518"	"23S"	"2518"	"0.0005"	"37050"	"37892"	0.0018
"23S_2519"	"23S"	"2519"	"0.0007"	"37064"	"37907"	0.0018
"23S_2520"	"23S"	"2520"	"0.0006"	"37141"	"37988"	0.0018
```

* **Nanocompore output** ([Check full input](example_input/Nanocompore/WT_IVT_nanocompore_results.tsv))
```
pos	chr	genomicPos	ref_id	strand	ref_kmer	GMM_logit_pvalue	GMM_logit_pvalue_context_4	KS_dwell_pvalue	KS_dwell_pvalue_context_4	KS_intensity_pvalue	KS_intensity_pvalue_context_4	MW_dwell_pvalue	MW_dwell_pvalue_context_4	MW_intensity_pvalue	MW_intensity_pvalue_context_4	TT_dwell_pvalue	TT_dwell_pvalue_context_4	TT_intensity_pvalue	TT_intensity_pvalue_context_4	GMM_cov_type	GMM_n_clust	cluster_counts	Logit_LOR
2515	NA	NA	23S	NA	ACATC	0.9999145725030434	0.9999997713740346	1.0	1.0	1.0	1.0	0.9998988895345491	0.9999961317288804	0.9997948160734869	0.9999943103243797	0.9998967130410295	0.999999547434852	0.9999171708631587	0.9999935368706631	full	2	Parent10_JW3718-90_1:2137/5121__JW3718_1:2284/5121	-0.06649521174331559
2516	NA	NA	23S	NA	CATCC	0.9999145725030434	0.9999997713740346	1.0	1.0	1.0	1.0	0.9998988895345491	0.9999961317288804	0.9997948160734869	0.9999943103243797	0.9998967130410295	0.999999547434852	0.9999171708631587	0.9999935368706631	full	2	Parent10_JW3718-90_1:2946/4290__JW3718_1:3033/4347	-0.015898013921205018
2517	NA	NA	23S	NA	ATCCT	0.9999145725030434	0.9999997713740346	1.0	1.0	1.0	1.0	0.9998988895345491	0.9999961317288804	0.9997948160734869	0.9999943103243797	0.9998967130410295	0.999999547434852	0.9999171708631587	0.9999935368706631	full	2	Parent10_JW3718-90_1:5699/1344__JW3718_1:5818/1374	0.001397467230704439
2518	NA	NA	23S	NA	TCCTG	0.9999145725030434	0.9999997713740346	1.0	1.0	1.0	1.0	0.9998988895345491	0.9999961317288804	0.9997948160734869	0.9999943103243797	0.9998967130410295	0.999999547434852	0.9999171708631587	0.9999935368706631	full	2	Parent10_JW3718-90_1:159/6944__JW3718_1:150/7091	0.07883939042455489
2519	NA	NA	23S	NA	CCTGG	0.9999145725030434	0.9999997713740346	1.0	1.0	1.0	1.0	0.9998988895345491	0.9999961317288804	0.9997948160734869	0.9999943103243797	0.9998967130410295	0.999999547434852	0.9999171708631587	0.9999935368706631	full	2	Parent10_JW3718-90_1:4179/2943__JW3718_1:4325/2949	-0.032296111509103595
2520	NA	NA	23S	NA	CTGGA	0.9999145725030434	0.9999997713740346	1.0	1.0	1.0	1.0	0.9998988895345491	0.9999961317288804	0.9997948160734869	0.9999943103243797	0.9998967130410295	0.999999547434852	0.9999171708631587	0.9999935368706631	full	2	Parent10_JW3718-90_1:3925/3177__JW3718_1:3981/3235	0.003922821293985011
```

### Usage

* Default command:
```
Rscript NanoConsensus.R -Epi_Sample ./example_input/Epinano/epinano_WT.tsv.per.site.var.csv -Epi_IVT ./example_input/Epinano/epinano_IVT.tsv.per.site.var.csv -NP_Sample ./example_input/Nanopolish/Nanopolish_WT_processed_perpos_median.tsv -NP_IVT ./example_input/Nanopolish/Nanopolish_IVT_processed_perpos_median.tsv -Tombo ./example_input/Tombo/WT_IVT_plus_Tombo_Output.tsv -Nanocomp ./example_input/Nanocompore/WT_IVT_nanocompore_results.tsv -ini_pos 50 -fin_pos 1492 -output output_name -fasta ./example_input/Reference.fa -chr 18S 
```

* Use the following command to change the Z-Score threshold (default = *5*) to identify putative modified sites from individual software data:
```
Rscript NanoConsensus.R -Epi_Sample ./example_input/Epinano/epinano_WT.tsv.per.site.var.csv -Epi_IVT ./example_input/Epinano/epinano_IVT.tsv.per.site.var.csv -NP_Sample ./example_input/Nanopolish/Nanopolish_WT_processed_perpos_median.tsv -NP_IVT ./example_input/Nanopolish/Nanopolish_IVT_processed_perpos_median.tsv -Tombo ./example_input/Tombo/WT_IVT_plus_Tombo_Output.tsv -Nanocomp ./example_input/Nanocompore/WT_IVT_nanocompore_results.tsv -ini_pos 50 -fin_pos 1492 -output output_name -fasta ./example_input/Reference.fa -chr 18S -plot --MZS_thr 4
```

* Use the following command to change the NanoConsensus score threshold (default = median(NanoConsensus_Score) * *5*) to filter the final putative modified sites:
```
Rscript NanoConsensus.R -Epi_Sample ./example_input/Epinano/epinano_WT.tsv.per.site.var.csv -Epi_IVT ./example_input/Epinano/epinano_IVT.tsv.per.site.var.csv -NP_Sample ./example_input/Nanopolish/Nanopolish_WT_processed_perpos_median.tsv -NP_IVT ./example_input/Nanopolish/Nanopolish_IVT_processed_perpos_median.tsv -Tombo ./example_input/Tombo/WT_IVT_plus_Tombo_Output.tsv -Nanocomp ./example_input/Nanocompore/WT_IVT_nanocompore_results.tsv -ini_pos 50 -fin_pos 1492 -output output_name -fasta ./example_input/Reference.fa -chr 18S  -plot --NC_thr 3
```

## Expected output 

By default, *NanoConsensus* generates all files listed below:
* **Raw_kmers** file: it contains results for all positions across the transcript - ZScores for all softwares, NanoConsensus score, kmer sequence and if the kmer contains the RRACH motif ([Check full output](example_output/m7G-100_16S_MZS-5_rep1_Raw_kmers.txt))
```
"Chr"	"Start"	"End"	"Epinano_RawScore"	"Nanopolish_RawScore"	"Tombo_RawScore"	"Nanocompore_RawScore"	"Epinano_Score"	"Nanopolish_Score"	"Tombo_Score"	"Nanocompore_Score"	"Epinano_Status"	"Nanopolish_Status"	"Tombo_Status"	"Nanocompore_Status"	"Merged_Score"	"Kmer"	"RRACH_motif"
"16S"	50	54	0.0369794364631788	1.82818269809314	1.61603605003496	4.7381541131775	-0.153261931378237	1.23445293486926	0.0486326005850889	0.0773097901992757	"NO"	"NO"	"NO"	"NO"	0.0110980252190282	"CTAACACAT"	FALSE
"16S"	51	55	1.19239687383871	0.593153945318812	1.98174190039624	4.77181029462782	0.0306194048103809	-0.606426948302039	0.0775030320334447	0.0780058429162199	"NO"	"NO"	"NO"	"NO"	0.0127467047756284	"TAACACATG"	FALSE
"16S"	52	56	0.40074893192265	2.34822804462824	2.09237821459094	5.37449733439723	-0.0953690705593066	2.00960980034528	0.086237150236644	0.0904701788397033	"NO"	"NO"	"NO"	"NO"	0.0137390036070691	"AACACATGC"	FALSE
"16S"	53	57	1.15732557192689	0.515553739466873	1.66024395928832	4.40949919767311	0.0250379087650485	-0.722094423704443	0.0521225677603841	0.0705127877760644	"NO"	"NO"	"NO"	"NO"	0.0110201820433646	"ACACATGCA"	FALSE
"16S"	54	58	0.529541201180387	0.526011294229817	1.51169295315049	3.84379605690622	-0.0748721542104668	-0.70650684961192	0.0403952966897964	0.0588133260086669	"NO"	"NO"	"NO"	"NO"	0.00989066788460121	"CACATGCAA"	FALSE
"16S"	55	59	0.17142513069125	0.518767445453894	1.02431823479139	1.21726635433112	-0.131865288831783	-0.717304214012041	0.00191978862191106	0.0044933450473563	"NO"	"NO"	"NO"	"NO"	0.0055430731860283	"ACATGCAAG"	FALSE
"16S"	56	60	0.133031821411215	0.573909937206806	0.646103643850517	0.582421724730132	-0.137975472720951	-0.635111225753392	-0.027938137932284	-0.0086360508088972	"NO"	"NO"	"NO"	"NO"	0.00333820080630401	"CATGCAAGT"	FALSE
"16S"	57	61	1.11383029056232	1.26490169831062	0.440136741511926	0.379594344462311	0.0181157608066584	0.394850894234233	-0.0441980728737736	-0.0128307794745501	"NO"	"NO"	"NO"	"NO"	0.00507170736103087	"ATGCAAGTC"	FALSE
```

* **Supported_kmers** file: it contains only those regions supported by two or more RNA detection algorithms and whose NanoConsensus score is higher than the specified threshold. All scores reported here are the maximum value found within the reported region. ([Check full output](example_output/m7G-100_16S_MZS-5_rep1_Supported_kmers.txt))
```
"Chr"	"Start"	"End"	"Epinano_RawScore"	"Nanopolish_RawScore"	"Tombo_RawScore"	"Nanocompore_RawScore"	"Epinano_Score"	"Nanopolish_Score"	"Tombo_Score"	"Nanocompore_Score"	"Epinano_Status"	"Nanopolish_Status"	"Tombo_Status"	"Nanocompore_Status"	"Kmer"	"RRACH_motif"
"16S"	511	536	177.76447909029	2.63628455239892	112.828762333929	600.28039307292	28.1315544965619	2.43897424160257	8.82825531428257	12.3938821290665	"YES"	"NO"	"YES"	"YES"	"AACTCCGTGCCAGCAGCCGCGGTAATACGG"	FALSE
```

* **NanoConsensus_Scores** pdf file: it contains tracks with both Z-Scores for all RNA detection algorithms and NanoConsensus scores across the transcript. Coloured, positions supported by each one of the softwares. 
![Zscores_tracks](/img/NanoConsensus_Scores.png)

Additionally, it also creates two directories:
* **Bedgraph_tracks**: one bedgraph track per software (Epinano, Nanopolish, Tombo, Nanocompore and NanoConsensus) can be found. Each contains the score per-position from one specific software. These tracks help visualizing scores across the transcript. ([Check output](example_output/Bedgraph_tracks))
![Bedgraph_tracks](/img/Bedgraph_tracks.png)

* **Kmer_tracks**: one bedgraph track per software (Epinano, Nanopolish, Tombo and Nanocompore) can be found. These tracks show all kmers supported by individual algorithms. They help visualizing which softwares are the main contributors to the final results. ([Check output](example_output/Kmer_tracks))
![Kmer_tracks](/img/Kmers_tracks.png)

## Required dependencies
* R (version 4.1.1)
* The following R packages: GenomicRanges, plyr, dplyr, VennDiagram, ggplot2, argparse, stringr, scales, ggnewscale, ggrepel, gtable
* bedtools (version 2.29.2)

## Citation

If you find this work useful please cite: 
Delgado-Tejedor A, Medina R, Begik O, Cozzuto L, Ponomarenko J and Novoa EM. Native RNA nanopore reveals antibiotic-induced loss of rRNA modifications located in the A and P sites. bioRxiv 2023. 

## Contact 

Please open an issue in the GitHub repo if you have any questions/doubts/suggestions about how to use this software. Thanks!

