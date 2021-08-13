# From Bam to consensus sequences and variant classification

---
## Initial Setup

*This document assumes [preproc htstream](./alignment) has been completed.*

---
## Building consensus sequences in R

### First lets setup our environment.

```bash
cd /share/workshop/alliance_covid/$USER/covid_swift
module load R
R
```

### In R Install needed packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("Biostrings", quietly = TRUE))
  BiocManager::install("Biostrings")

if (!requireNamespace("Rsamtools", quietly = TRUE))
  BiocManager::install("Rsamtools")

if (!requireNamespace("GenomicAlignments", quietly = TRUE))
  BiocManager::install("GenomicAlignments")

if (!requireNamespace("DECIPHER", quietly = TRUE))
  BiocManager::install("DECIPHER")

if (!requireNamespace("phangorn", quietly = TRUE))
  BiocManager::install("phangorn")
```

### Lets download the IUPAC Code map

```r
download.file("https://github.com/ucdavis-bioinformatics-training/2021-Alliance-Makerere_Covid/raw/main/data/IUPAC_CODE_MAP_extended.rds", "IUPAC_CODE_MAP_extended.rds")
IUPAC_CODE_MAP_extended <- readRDS("IUPAC_CODE_MAP_extended.rds")
```

### Build our consensus sequences from the aligned bam files

```r
#load the libraries
library("Biostrings")
library("Rsamtools")
library("GenomicAlignments")

# running samples in swift_samples.B5
samples <- readLines("swift_samples.B5.txt")

# load our reference sequence
sarscov2_seq <- readDNAStringSet(file.path("resources", "NC_045512.2.fasta"))

## Load the se and pe bam files, first specifying desired params
param <- ScanBamParam(what=c("qname","flag", "seq", "qual"),
                      flag=scanBamFlag(isUnmappedQuery = FALSE,
                                       hasUnmappedMate = FALSE,
                                       isSecondaryAlignment = FALSE,
                                       isSupplementaryAlignment = FALSE))

ga_align_se <- lapply(samples, function(x) readGAlignments(file.path("02-Bowtie2", paste0(x,"-se.bam")), param=param))
names(ga_align_se) <- basename(samples)

ga_align_pairs <- lapply(samples, function(x) readGAlignmentPairs(file.path("02-Bowtie2", paste0(x,"-pe.bam")), use.names = TRUE, param=param))
names(ga_align_pairs) <- basename(samples)

# First process is to combine the pairs and singles into a single genomic alignments object
ga_combined <- list()
for(s in basename(samples)){
  ga_combined[[s]] <- c(ga_align_se[[s]], first(ga_align_pairs[[s]]), second(ga_align_pairs[[s]]))
}

# Finally compute the consensus sequence using the basic process.
#  sequenceLayer -> consensusMatrix -> consensusString
sample_af <- DNAStringSet()
for(s in basename(samples)){
  gaa <- ga_combined[[s]]

  qseq_on_ref <- sequenceLayer(mcols(gaa)$seq, cigar(gaa))
  pos_by_ref <- start(gaa)
  cm <- consensusMatrix(qseq_on_ref, shift=pos_by_ref)
  cm <- cbind(cm,matrix(0,nrow=dim(cm)[1], ncol=width(sarscov2_seq)-dim(cm)[2]))

  col_sums <- colSums(cm)
  col_sums[col_sums == 0] <- 1  # to avoid division by 0
  cm_prop <- cm / rep(col_sums, each=nrow(cm))

  idx <- colSums(cm_prop) == 0
  cm_prop["N", idx] <- 1

  # How does changing the threshold change the results? here because their are 5 single character bases (A,C,G,T,N) in IUPAC_CODE_MAP_extneded, the threshold has to be in the range [0-0.2]
  dss <- DNAStringSet(consensusString(cm_prop, threshold=0.15,  ambiguityMap=IUPAC_CODE_MAP_extended) )
  names(dss) <- s
  sample_af <- c(sample_af,dss)
}
sample_af

# Finally write out the consensus object
dir.create("03-ConsensusSeqs", showWarnings = FALSE)
writeXStringSet(sample_af, file.path("03-ConsensusSeqs", "swift_samples.B5.fasta"))
q()
```

**quit R**

## Variant classification with Pangolin

We are going to use [Pangolin](https://cov-lineages.org/index.html) to type each variants

```bash
module load anaconda3
cd /share/workshop/alliance_covid/$USER/
git clone https://github.com/cov-lineages/pangolin.git
cd pangolin
conda env create --prefix /share/workshop/alliance_covid/$USER/conda_pangolin -f environment.yml

conda activate /share/workshop/alliance_covid/$USER/conda_pangolin
pip install .

# This should be run each time before you use Pangolin
pangolin --update

cd /share/workshop/alliance_covid/$USER/covid_swift
pangolin 03-ConsensusSeqs/swift_samples.B5.fasta -t 16 --outfile 04-Pangolin/swift_samples.B5.pangolin
```

Transfer the file

04-Pangolin/swift_samples.B5.pangolin

to your computer and view the result on your computer. What variants do we have????


Now do the same for the other dataset. What do you see about the Variants??  B5 was back in April and B19 is recent.

<!--
## Phylogenetics
```r
library(DECIPHER)
library(phangorn)

consensus_seqs <- readDNAStringSet(file.path("03-ConsensusSeqs", "swift_samples.B5.fasta"))
alignment = AlignSeqs(consensus_seqs, anchor=NA, processors=16)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm)

fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))

``` -->
