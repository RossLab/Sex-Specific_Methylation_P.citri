# Sex Specific DNA Methylation in *Planococcus citri*.
This is the repository for all of the code used in the below paper:

Bain, S. A., Marshall, H. and Ross. L. (2020). Sex-specific expression and DNA methylation in a species with extreme sexual dimorphism and paternal genome elimination.

Pre-print [here](https://doi.org/10.1101/2020.06.25.171488).

Each folder contains numbered scripts representing the sequence in which the scripts were executed. Scripts are then also grouped by letters if appropriate, e.g. A01.script and A02.script belong to the same pipeline.

---

## Project Structure and Scripts

**Alternative_splicing_etc**

Sex-specific alternative splicing events and identification of sex-biased, extremely sex-biased and sex limited genes. These scripts are based on Stevie's RSEM counts.

**De novo trancriptome**

Creation of a de novo transcriptome from the RNA-seq data in order to better annotate alternatvie splicing events. New analysis based on reviewer comments.

**Differential_methylation**

Full pipeline using whole genome bisulfite sequencing data from fastqc to differential methylation identification between sexes. It also includes annotation of differentially methylated sites. 

**GO_analysis**

GO enrichment analysis of gene sets.

**Genome_formating_and_extra_annotations**

Creation of putative promotors, intron annotation, exon 1-3 annotation and intergenic annotation. Also a script for calculating the number of CpGs per given annotation.

**Methylation_paired_with_expression**

Relationship of methylation with general expression, differential expression and alternative splicing. 

**Methylation_proxies**

Annotation of CpG islands. (Not used in the paper).

**Motif_enrichment_HOMER**

Identification of enriched motifs for methylation and unmethylation annotations. (Not used in the paper).

**Transposable_elements**

De novo annotation of transposable elements in *P. citri*.

---

If you have any questions/comments etc. please contact:
- email: hollie_marshall@hotmail.co.uk
- twitter: @mooholl
