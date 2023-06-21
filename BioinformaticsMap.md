# Hacky Hour Bioinformatics Links
Last updated 21.6.23

Guides to navigating the bioinformatics workflows 'help desk' space.  Always a work in progress.

# Workflow - workflow management tools

Current list (will no doubt grow)

[Snakemake](https://snakemake.github.io)

[cwltool](https://pypi.org/project/cwltool/) - Common Workflow Language, python implementation

[Toil](http://toil.ucsc-cgl.org) - also python

[Nextflow](https://www.nextflow.io)

Containers:

[Docker](https://www.docker.com)

[Singularity](https://singularity.hpcng.org) - HPC container system

# Workflow - common bioinformatics pipeline phases

## General steps that occur in most workflows.

(1) Seq & QA checks (may be done offset, by the big labs, but you supply data)

(2) Align 

(3) Index (reference index for counting - may be optional)

(4) Count 

(5) Analyse 

Various combinations of tools can be used, but since some tools are specialised, there are situations in which the same tool can be used at an earlier or later part of the process for the same workflows.  

As between the abbreviated tables with seq/align/count/analyse, it may be possible to skip from step 2 in one table to step 3 in the next.

Steps 2 to 4 are most data intensive but knowing the basis for 1 means knowing the scientific goals well.

Also, step 2 (alignment) is often performed with the assistance of a feature-based index, so this is something that must be obtained or created before performing the alignment.

In step 3, some packages are in python (e.g. https://htseq.readthedocs.io/en/master/counting.html), with differential expression (not variation calling or gene prediction) being the end goal. 

Noteworthy points about bioinformatics workflows:
1. You can mix and match a bit with software at each stage, but the more common approach is for vertical setups to be developed.  e.g. in R ecosystem, or Bioconductor etc.
2. Within workflows, there may sometimes be choices as to computing environment but alignment will generally push requirements to HPC.
3. Size of genomes and place in biological space (kingdom, species) often affect both assumptions of software development and HPC requirements e.g. different requirements for
(a) Mouse, Human, common medical context.
(b) Small plant v large plant genomes

The permutations involved in bionformatics workflows are considerable.  Even amongst these tools, and the named stages, researchers may be:
1. Attempting different analysis goals, including variant calling or differential gene expression (variant calling may qualify for stage 4 'analysis', omitting stage 3).  
2. Choosing to perform the work on a PC, rather than HPC
3. To work on a large rather than a small genome
4. To work on human or non-human animal genomes, or plant genomes, or bacteria.
5. To align short or long reads (some tools are for 'short reads' i.e.). 
6. DNA (with introns, gaps) or RNA (no introns, spliced) TopHat2, 
7. Single end or paired end reads
8. Whether a software package outputs files, structures, or some custom object.

(nb many of these choices find their way into paramaters of functions like featureCounts in RSubread)

It's necessary to understand the relationship between the type of genome (plant, animal, model species (human, mouse)), whether it is a DNA or RNA genomic study, and the consequences this has immediately for the alignment tools and constraints.  In the following workflows, the 4 stages of analysis may be the same, but these genomic considerations suggest that the use of appropriate alignment tools will be suggested.  If there are vertical workflows that pipe information from these aligners to specific counting programs, then these should also shop up in the workflow.

## RNA-seq studies

As ([Soneson et al](https://f1000research.com/articles/4-1521)) noted:

"most RNA-seq studies can be classified as either: 

1) differential gene expression (DGE) studies, where the overall transcriptional output of each gene is compared between conditions; 

2) differential transcript/exon usage (DTU/DEU) studies, where the composition of a gene’s isoform abundance spectrum is compared between conditions, or 

3) differential transcript expression (DTE) studies, where the interest lies in whether individual transcripts show differential expression between conditions."

The original counting methods used target genes to provide the static, non-overlapping areas that would be counted.  This required gene alignment, and then counts were expressed in terms of the counts for the target area.

An original measure, namely reads/fragments per kilobase per million reads, has been queried recently.

# Differential Expression workflows 

## Pseudo-alignment

Pseudo-alignment is a broad term used to describe the transcript counting methods of tools such as kallisto and Salmon.

These workflows illustrate possible RNA-seq transcriptome quantification and analysis pipelines

| | Stage  | WF1 | WF2 | WF3 |
|:-----|:-----|:-----|:-----| :-----|
| 1 | Seq & QA |  | 
| 2 | Align | | [Salmon](https://github.com/COMBINE-lab/salmon) | [Salmon](https://github.com/COMBINE-lab/salmon) |
| 3 | Count |  [Kallisto](http://pachterlab.github.io/kallisto/) | | [Alevin-Fry](https://github.com/COMBINE-lab/alevin-fry)  |
| 4 | Convert | | [wasabi](https://github.com/COMBINE-lab/wasabi) | [tximport](https://github.com/mikelove/tximport) |
| 5 | Analyse |  Kallisto/[Sleuth(R)](http://pachterlab.github.io/sleuth/) | [Sleuth(R)](http://pachterlab.github.io/sleuth/) |[Seurat 4.0](https://satijalab.org/seurat/) |

"On benchmarks with standard RNA-Seq data, kallisto can quantify 30 million human reads in less than 3 minutes on a Mac desktop computer using only the read sequences and a transcriptome index that itself takes less than 10 minutes to build" [RefLink](http://pachterlab.github.io/kallisto/)

Kallisto/Sleuth (Pachterlabs):

[GencoreNotes](https://learn.gencore.bio.nyu.edu/rna-seq-analysis/salmon-kallisto-rapid-transcript-quantification-for-rna-seq-data/)

[Geneviève Boucher's notes (2018)](https://bioinfo.iric.ca/understanding-how-kallisto-works/)

Kallisto will count the transcriptome (transcripts/sequences) information, but does so without alignment to, or referencing the genome?

Sleuth requires the kallisto output to work: [SciLifeLab_Sleuth](https://scilifelab.github.io/courses/rnaseq/labs/kallisto), but wasabi can convert Salmon output.

Salmon produces gene-level abundance estimates, and Alevin-Fry can do the transcript-evel.  If gene-level is suitable, it may be possible to take the Salmon data.

tximport ([Soneson et al](https://f1000research.com/articles/4-1521)) is an R package 'to help users integrate transcript-level abundance estimates from common quantification pipelines into count-based statistical inference engines.'.  It is useful for importing that output into *any* statistical analysis package that is R-based, like Seurat.

tximport can also be used to help translate the outputs of these workflows to common statistical inference programs like DESeq2, edgeR, limma and Seurat.


nb: 

1. Alevin-Fry was released 1 July 2021: [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.06.29.450377v1)
"Alevin-fry has been designed as the successor to alevin. "
"we encourage users of alevin to migrate to alevin-fry when feasible."
"The alevin-fry workflow quantifies single-cell data based on a reference index created by salmon."

2. What use does Alevin-fry make of Salmon?

From the 1 July 2021 bioRxiv paper:
```
"It makes use of salmon (4) for basic barcode and UMI parsing and the mapping of the reads to the constructed reference index. The output of salmon, when configured to produce output for alevin-fry, is a RAD (Reduced Alignment Data) format file, which is a chunk-based binary file optimized for machine parsing, that encodes the relevant information necessary for subsequent (post-mapping) processing of the data (Section S1)."
```

3. You can pipe the Alevin information into Seurat: [Tute](https://combine-lab.github.io/alevin-tutorial/2018/alevin-seurat/) (tutorial by )
The author of that tutorial is [Avi Srivastava](https://twitter.com/k3yavi), from the Satija Lab.

# RNA-sequencing and full alignment to genomes

## Some general RNA seq workflows that utilise Nextflow

The NextFlow-Core community is working on capturing some general bioinformatics workflows in Nextflow scripts (pipelines), from initial sequencing through to analysis.  Some of the pipelines have as many as 15 separate steps.  

Follow this general link:
[RNASeq_Nextflow](https://nf-co.re/rnaseq)

The above-mentioned pipeline includes use of STAR or HISAT2 at Step 7 (it is noted in link above that HISAT2 does not perform quantification at the alignment stage, and that other downstream software is needed to do this).  

For some possible alternatives RNA seq pipelines see below.

## Other

| | Stage  | WF4 | 
|:-----|:-----|:-----|
| 1 | Seq & QA | CellRanger/10x | 
| 2 | Align | CellRanger |  
| 3 | Count |  | 
| 4 | Convert | | 
| 5 | Analyse | Seurat | 

CellRanger is commercial, supplied with 10X data.

| | Stage  | WF5 | WF6 |  WF7 |
|:-----|:-----|:-----|:-----|:-----|
| 1 | Seq & QA | Other | | |
| 1.5  | Convert | Dropseq/Picard | Dropseq/Picard | Dropseq/Picard |  
| 2 | Align | [STAR](https://github.com/alexdobin/STAR) | STAR | STAR |  
| 3 | Index |  |  |
| 4 | Count | STARsolo | | ?? |   
| 5 | Analyse | | Seurat |  |  


WF4-WF6:

STAR is a long-read RNA-Seq aligner.  Another is [GMAP](http://research-pub.gene.com/gmap/).

STAR (and some other aligners) usually require unmapped BAM or SAM files - this means if you have FASTQ files as your sequencing output data format you have to insert another step to convert them to SAM or BAM, and this requires use of some of the FASTQ tools (e.g. Picard and DropUtils, which are relaetd to the DropSeq workflow.  There may be other tools.).  In any case, SAM and BAM are formats designed for the output of the alignment process, not just for the step preceding alignment.   Alignment may not immediately produce analysis ready BAM/SAM data: there may still be a need for removal of duplicates or bias.

An annotation file (GTF/GFF/GFF3) is needed for most alignment steps in Stage 2 (WF1-WF3), and also for some of the non-alignment counters like Kallisto and Salmon (WF4)

In step 3, the alignment free methods require tximport if DESEQ2 is next step.   The python program will not immediately interface with the R programs, but can be incorporated by reading it back into the DESeq2 functions.

| | Stage  | WF8 |
|:-----|:-----|:-----|
| 1 | Seq & QA |  |
| 2 | Align | [HISAT2](https://daehwankimlab.github.io/hisat2/) | 
| 3 | Count | HTSeq-count (py) |  
| 4 | Analyse |  DESeq2  |   

nb: HISAT2 is able to deal with spliced alignment (at cost of RAM needed).

HISAT2 superseded tophat, which was itself based on bowtie.  However, there is an important contextual limitation for HISAT2 - it's graph-based mapping algos for reads (DNA or RNA) are designed with reference to  the human genome.  So it's not useful as a general purpose aligner.  

HTSeq-count will count the transcriptome information, but does so with extra alignment step so that the transcriptome counts can be referenced back to the genome (exons overlapping with transcripts?).

Quick links:

Useful distinctions (as mentioned in beginner guide to DESeq2, by M. I. Love, W. Huber, S. Anders: function (https://workshop.veupathdb.org/bop/pdfs/beginner_DeSeq2.pdf), package, output, input function to next package).  Useful workshop: https://workshop.veupathdb.org/bop/ (2018, 2 lectures on iCloud).  Lecture 2 on counts has some good observations on counting, normalisation and why RPKM needs correction (also suggests Kallisto, Salmon offer some advantages)

| | Stage  | WF9 | WF10 | WF11 | 
|:-----|:-----|:-----|:-----|:-----|
| 1 | Seq & QA |   | | | 
| 2 | Align |  Subread (R) |Subread (R) | Subread (R) |
| 3 | Count |   featureCount (func. in RSubread) | scPipe (R) - uses subRead | |
| 4 | Analyse |   [Limma (R)](http://bioconductor.org/packages/release/bioc/html/limma.html) | singleCellExperiment (R) (R/Bioconductor) | |

featureCount will count the transcriptome information, but does so with extra alignment step so that the transcriptome counts can be referenced back to the genome (exons overlapping with transcripts, also multi-mapping reads of transcriptomes ignored for gene counts here, but contribute to transcript counts?).

Some HPC-based advice - Galaxy: https://help.galaxyproject.org/t/how-can-i-improve-very-low-assigned-rate-in-featurecounts/758
This [Galaxy tutorial](https://www.melbournebioinformatics.org.au/tutorials/tutorials/galaxy-workflows/galaxy-workflows/) by Simon Gladman of Melbourne University is also very helpful.

LIMMA: Linear Models for Microarray and RNA-seq Data

A Biconductor software package for the analysis of gene expression studies, especially the use of linear models for analysing designed experiments and the assessment of differential expression. The analysis methods apply to many different technologies, including microarrays, RNA-seq, quantitative PCR and many protein technologies.

LIMMA is available as part of Bioconductor project. To install limma from the R command line, use the BiocManager package from CRAN:

> library("BiocManager")
> install("limma")
> install("statmod")

http://bioinf.wehi.edu.au/limma/

LIMMA is a command driven package but menu driven interfaces are also available. See limmaGUI for two-colour arays or affylmGUI for Affymetrix arrays.

| | Stage  | WF12 | WF13 | WF14 |
|:-----|:-----|:-----|:-----| :-----| 
| 1 | Seq & QA |   | |  |
| 2 | Align |  Bowtie2 | Bowtie2 | Bowtie2 | 
| 3 | Count | MOSAIK, SHRiMP2, and Novoalign | packages producing SummarizedExperiment (Bioc) | [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) | 
| 4 | Analyse | | edgeR | [CummeRbund](http://compbio.mit.edu/cummeRbund/) | 

WF12-WF14 are examples of the use of a short-read aligner (Bowtie2).  You can mix this with Bioconductor packages for counting and analysis

WF14 - Cufflinks/CummeRbund (Colin Trapnell) software was last updated about 2014.

For Bowtie2 paper:
"Aligning sequencing reads to a reference genome is the first step in many comparative genomics pipelines, including pipelines for variant calling1, isoform quantitation2 and differential gene expression3. "
(Bowtie2 differs from Bowtie in that it can align even if the reads run across gaps.)

Quick links:

| | Stage  | WF15 |
|:-----|:-----|:-----|
| 1 | Seq & QA |   | 
| 2 | Align |  BWA | 
| 3 | Count |  RSEM  | 
| 4 | Analyse |  Blast2Go (DE) | 


# Gene prediction/annotation workflows

| | Stage  | WF17 | WF18| WF 19 |
|:-----|:-----|:-----|:-----|:-----|
| 1 | Seq & QA |  | | |
| 2 | Align | Cufflinks | | |
| 3 | Count | [Cufflinks](http://cole-trapnell-lab.github.io/cufflinks/) | | 
| 4 | Train |[BRAKER1&2](http://topaz.gatech.edu/Genemark/braker.html) | |[BIND](https://github.com/eswlab/orphan-prediction/tree/master/evidence_based_pipeline) |
|5 | Combine | [TSEBRA](https://github.com/Gaius-Augustus/TSEBRA) | |
| 6 | Predict | [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus)| MAKER | |
| 7 | Compare | | |[Mikado](https://github.com/lucventurini/mikado) | 


WF17 to 19 - MAKER and BRAKER are *ab initio* tools that infer gene structures by machine learning.  This means that they attempt to use the target genome as input and some are easily trained for a new genome.

GeneMark-ET is an self-training algorithm, developed at Georgia Institute of Tech, that integrates information on splice-aligned RNA-seq reads with genomic information.

WF18 - [AUGUSTUS](https://doi.org/10.1093/bioinformatics/btn013) is another algorithm for gene prediction.  For model parameter estimation, AUGUSTUS requires an expert curated training set of genes.  It incorporates several different evidence sources, including gene and transcript annotations from related species.  Some essential tools (including for BRAKER) include the C++ tool bam2hints.  

The output is a .gtf file.  

EXAMPLE HACKY HOUR PROBLEM:
For NCBI genome information you may wish to output from AUGUSTUS as a gff3 file.  There are some instructions [here](https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/).  

*However, as far as I know, Augustus currently (15 March 2022) has some issues with producing a valid gff3 that will pass the NCBI genometools validator.   Some suggested fixes will be included here in due course.  One tool that may help to prepare a suitable gff file is table2asn.  See [NCBI link](https://www.ncbi.nlm.nih.gov/genbank/table2asn/).  As advised by P.Bayer 15.3.22*

WF18 & 19 - *Automated genome annotation*

BRAKER1 was released in 2015, as a tool to train AUGUSTUS.  BRAKER is a collection of Perl and Python scripts and a Perl module.  It is "a combination of GeneMark-ET R2 and AUGUSTUS R3, R4, that uses genomic and RNA-Seq data to automatically generate full gene structure annotations in novel genome." :[RefLink](https://github.com/Gaius-Augustus/BRAKER#what-is-braker).  "If Biopython is installed, BRAKER can generate FASTA-files with coding sequences and protein sequences predicted by AUGUSTUS..."

BRAKER2 is a fully automated training tool, based on GeneMark-EX and AUGUSTUS.

For the latest on the BIND pipeline (direct inference) see Dept of Genetics, Iowa State Uni team: [BIND/MIND bioRxiv](https://www.biorxiv.org/content/10.1101/2019.12.17.880294v3).  The Direct Inference pipeline, with SnakeMake ([python script](https://github.com/eswlab/orphan-prediction/blob/master/evidence_based_pipeline/snakemake_config/slurm/slurm-submit.py)), is available from: https://github.com/eswlab/orphan-prediction/tree/master/evidence_based_pipeline  (there are instructions there for installing from the conda/bioconda channel).  There are some basic single node and HPC pipeline instructions too.

TSEBRA is the Transcript Selector for BRAKER. 

The goal is to try and use RNA-seq data to annotate protein-encoding genes for future reference, using machine learning techniques.  The latest methods experiment with hybrid approaches to predicting genes.   This method involves a step to align the transcript using Cufflinks.

Mikado is a python based prediction comparision tool.  It requires, as a minimu, SQLite access.

# Variant calling workflows (TO DO)

Variant calling software (to identify variants within reads in the experiment that are in the same aligned regions) includes GATK-Haplotype caller, playtpus, SAMTools, bcftools.  [Past this, the goals of variant calling is annotation, which will feed back into subsequent indexing for alignment, and differential gene expression work].  The meaning of variant itself varies - it could be a single nucleotide (SNPs) or a deletion or insertion (indels), or structural variations (SVs).

(Data Carpentry-based course on bioinformatics](https://hoytpr.github.io/bioinformatics-semester/about/acknowledgments/), by Peter Hoyt, including [variant calling workflow] (https://hoytpr.github.io/bioinformatics-semester/materials/genomics-variant-calling-workflow/)

[Haplotype Caller](https://hpc.nih.gov/training/gatk_tutorial/haplotype-caller.html#optimized-script-2)
https://hpc.nih.gov/training/gatk_tutorial/workflow-overview.html

BWA, GATK and ENSEMBL

[GATK](https://gatk.broadinstitute.org/hc/en-us) is a diverse, Java-based tool.

For now, see these Melbourne Bioinformatics resources:

1. [Variant Calling Pipeline](https://github.com/melbournebioinformatics/variant_calling_pipeline)
2. [Variant Calling Using Galaxy](https://www.melbournebioinformatics.org.au/tutorials/tutorials/variant_calling_galaxy_1/variant_calling_galaxy_1/)

Trimming may be required.  See [TrimmingIntro](https://dnatech.genomecenter.ucdavis.edu/faqs/when-should-i-trim-my-illumina-reads-and-how-should-i-do-it/)

[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

# PhyloFaN (False negative tweaking)

General Visualisation tool developed in the context of human genetic diversity studies.  

Developed by [Brenna Henn Lab](http://hennlab.ucdavis.edu/news.html).  Population genetics.

[Henn Lab Software](http://hennlab.ucdavis.edu/software.html)

[2016 - PhlyoFaN](https://www.biorxiv.org/content/10.1101/066043v2.full.pdf)

TItle: False Negatives Are a Significant Feature of Next Generation Sequencing Callsets
Single Nucleotide Variants (SNV's)  - analysis and tuning of false negative (FN) rate in haploid human NGS experiments.
[PhyloFaN Software](https://ecoevo.stonybrook.edu/hennlab/data-software/)

## Apricot

[Preprint_Apricot (Lin et al)](https://www.biorxiv.org/content/10.1101/2021.03.09.434643v1)

Title: Admixed Populations Improve Power for Variant Discovery and Portability in Genome-wide Association Studies

[Apricot Simulation Software](https://github.com/menglin44/APRICOT)

[Henn Lab Software](http://hennlab.ucdavis.edu/software.html)

# Metagenomics

Metagenome assembly and analysis 

## Some useful tools

See [Wei Shen Github](https://github.com/shenwei356), a bioinformation from China who works on microbial genomics and metagenomics.  These tools include [SeqKit](https://bioinf.shenwei.me/seqkit/) for FASTA/Q file manipulation.  SeqKit is written in Go Language, and in addition to Shen's code, makes use of BAM and watch functions from [Thist](https://github.com/botond-sipos/thist) go-lang package by [Botond Sipos](https://bsipos.github.io), from Oxford Nanopore.

Citation:
W Shen, S Le, Y Li*, F Hu*. SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PLOS ONE. doi:10.1371/journal.pone.0163962.


## One workflow example.

Na-Kyung Kim, Sang-Hoon Lee, Hyeokjun Yoon, Garam Jeong, You-Jung Jung, Moonsuk Hur, Byoung-Hee Lee, Hee-Deung Park,
Microbiome degrading linear alkylbenzene sulfonate in activated sludge,
Journal of Hazardous Materials,
Volume 418,2021, 126365, ISSN 0304-3894,
https://doi.org/10.1016/j.jhazmat.2021.126365.

 - The quality-filtered metagenomic sequences from sample datasets included together in the de novo assembly by MEGAHIT.
 - MAGs constructed (binning) using MetaBAT
 - Complete MAGs with conserved Phylogenetic marker genes, used to identify taxonomic affiliation of MAGs - PhyloPhlAn 3 with NCBI reference genomes.
 - PCA performed using R 'vegan' package. [CRAN vegan](https://rdrr.io/cran/vegan/man/vegan-package.html) [github vegan](https://github.com/vegandevs/vegan/)
 - Clustering and dendrogram preparation
 - Open reading frams on MAGs identified using prokka pipeline. [prokka github](https://github.com/tseemann/prokka)  [prokka paper 2014](https://research.monash.edu/en/publications/prokka-rapid-prokaryotic-genome-annotation)
 - The functional orthologs of annotated genes defined using GhostKOALA algorithm based on the KEGG database.

## Some Kraken2 resources

Kraken2 produces at least 2 formats.  The taxonomic trees in output files may need further work for statistical analysis.  For downstream processing (.report2 etc)

Good source for Kraken2Tools (there is a 2022 paper in Nature too):

[Jennifer Lu Kraken2Tools Github](https://github.com/jenniferlu717/KrakenTools)

Paper:
 ```
 [Lu J, Rincon N, Wood D E, Breitwieser F P, Pockrandt C, Langmead B, Salzberg S L, Steinegger M. 
 - Metagenome analysis using the Kraken software suite. Nature Protocols, doi: 10.1038/s41596-022-00738-y (2022)]
 ```
  [2022 Nature Paper, Jennifer Lu et al](https://www.nature.com/articles/s41596-022-00738-y)

Other install locations for [Kraken2Tools](https://bioconda.github.io/recipes/krakentools/README.html)

# MAGpy - workflow (Mick Watson Lab)

MAGpy = "a scalable and reproducible pipeline that takes multiple genome
assemblies as FASTA and compares them to several public databases, checks quality, suggests a
taxonomy and draws a phylogenetic tree."

[Stewart et al MAGpy 2018 paper](http://dx.doi.org/10.1093/bioinformatics/bty905)

[MAGpy_Github](https://github.com/WatsonLab/MAGpy)

## MAGpy internal workflow

As per 2018 paper, uses open source MAG analysis workflow, which can be defined by [Snakemake](https://snakemake.github.io):
 
1. [CheckM](http://ecogenomics.github.io/CheckM/) (Parks et al., 2015) is run to assess the completeness and contamination
of MAGs. [CheckM Paper](https://pubmed.ncbi.nlm.nih.gov/25977477/)
2. CheckM conservative assignment of taxonomic level to MAGs
to the MAGs, though in our experience this is often a conservative
3. predict protein coding sequences of MAGs using [Prodigal](https://github.com/hyattpd/Prodigal) (Hyatt et al., 2010). 
4. [DIAMOND aligner](https://github.com/bbuchfink/diamond) (Buchfink et al., 2015) BLASTP is used to compare the proteins to [UniProt](https://www.uniprot.org)
(UniProt Consortium, 2018). [DIAMOND paper](https://www.nature.com/articles/nmeth.3176/) [DIAMOND Docker](https://hub.docker.com/r/buchfink/diamond/)
5. Reports of the DIAMOND results at the level of the MAG and for each contig within each MAG. 
6. Proteins are compared to protein families in Pfam (Finn et al., 2014) using PfamScan [Pfamscan online tool](https://www.ebi.ac.uk/Tools/pfa/pfamscan/)
7. Create a tree using [PhyloPhlan](https://huttenhower.sph.harvard.edu/phylophlan) (Segata et al., 2013)
8. Visualise tree using [GraPhlAn](https://pypi.org/project/graphlan/) (Asnicar et al., 2015). 
9. The MAG genome sequences are also compared to over 100 000 public genomes using MinHash signatures as implemented in [Sourmash](https://sourmash.readthedocs.io/en/latest/index.html) (Brown and Irber, 2016). 

[JOSS entry Sourmash](https://www.theoj.org/joss-papers/joss.00027/10.21105.joss.00027.pdf)

[Sourmash github](https://github.com/dib-lab/sourmash/)

Per authors, [UniProt](https://www.uniprot.org) Hits at step 4:

(a) provide annotation of the putative proteins and may predict function; 

(b) protein hits can help define the closest sequenced genome;

(c) length of the predicted protein and that protein’s hits can be
used to detect truncated genes and proteins in the MAG annotation.

Both PhyloPhlan and GraPhlan are able to be installed via conda (python).

## PhyloPhlan 3 (from Huttenhower lab)

[PhyloPhlan3](https://huttenhower.sph.harvard.edu/phylophlan):

```
PhyloPhlAn 3.0 is released in Bioconda and can be installed using the following command:

conda install -c bioconda phylophlan
```

Bioconda is a channel for the conda package manager, with a focus on bioinformatics software.  It handles all the dependencies for installation (as an alternative to pip). 

[Bioconda Paper](https://www.nature.com/articles/s41592-018-0046-7)

Grüning, B., Dale, R., Sjödin, A. et al. Bioconda: sustainable and comprehensive software distribution for the life sciences. Nat Methods 15, 475–476 (2018). https://doi.org/10.1038/s41592-018-0046-7

# Visualisation of phylogenetic trees/dendrograms and fan-like visuals

## Dendrograms and phylogenetic fans.  

Dendrograms and phylogenetic fans/trees: 
[STHDA_dendrograms](http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning)
(*STHDA = Statistical tools for high-throughput data analysis.*)

Dendrogram options in R:
1. For hierarchical: ggraph
2. For clustering: hclust and plot.clust() R functions
3. Dendextend
4. ggdendro (extends ggplot2)
5. pvclust and dendextend

[Phylogenetic fan tip colour](https://stackoverflow.com/questions/46416143/apply-ifelse-function-to-colour-tips-of-phylogenetic-fan)

## Phytools

[PhyTools Blog](http://blog.phytools.org/2017/06/custom-tip-spacing-of-fan-tree-plot-in-r.html)

## Phylosignal

[CRAN phylosignal](https://cran.r-project.org/package=phylosignal)

[RefManual](https://cran.r-project.org/web/packages/phylosignal/phylosignal.pdf)

## GGTree (Bioconductor package)

[GGtree Home](https://guangchuangyu.github.io/software/ggtree/)

[GGtree vignette](https://bioconductor.org/packages/devel/bioc/vignettes/ggtree/inst/doc/ggtree.html)

[GGTree article](https://4va.github.io/biodatasci/r-ggtree.html)

# Things to know about Bioconductor:

1. It uses the SingleCellExperiment package (which creates an object of the same name) to store experimental data.
2. It doesn't care how you get the count matrix but refers to CellRanger, Rsubread, alevin (Salmon) and scPipe as possible programs to use to do the alignment and/or counting.

https://bioconductor.org/books/release/OSCA/data-infrastructure.html#background

In this respect, Bioconductor's SingleCellExperiment is another R package in the same space as Seurat package, but adopts a different set of R objects and functions to get the job done.

John Hopkins University has been operating with its own universe of R packages called 'Bioconductor'; a very regulated update and maintenance systems for the packages.

https://bioconductor.org/books/release/OSCA/index.html

Bioconductor has its own installation manager (BiocManager) to help with installation.
