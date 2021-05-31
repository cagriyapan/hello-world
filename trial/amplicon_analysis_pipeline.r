####16 S rRNA Amplicon sequencing Analysis on Ollie####

# Script by: Batuhan Cagri Yapan.
# contact: byapan@mpi-bremen.de
# Last update: 31/05/2021
# Based on David Benito's and Marwa Baloza's scripts & DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
# Script works in R environment (R 4.0.6) on Linux Distro CentOS.
# Cutadapt is working on terminal; an R function does it automatically.
# R 4.0 and Cutadapt 3.2 are installed on Ollie, they must be called on BASH by module function e.g.:
# module load bio/R/4.0.0
# module load bio/cutadapt/3.2
# Working interactively is not recommended on Ollie's login nodes (ollie0 or ollie1)
# To work interactively ask for a proper node e.g.:
# salloc -N 1-1 -c 36 -p fat --mem 425G
# Then you can start R and work, but beware that it is hard to see graphs etc.
# If you want to put your script on queue please check script  "amplicon_analysis.sl"
# Then send it to queue via:
# sbatch amplicon_analysis.sl


##########################################################################
### Four main parts in this script:
### 1)Cutadapt - for trimming primers and adapters of paired-end Illumina reads.
### 2)DADA2 - 16 S analysis from paired-end Illumina reads
### 3)Phyloseq - Analyses of microbial community structure !!not completed yet!!
##########################################################################

### Part 1 - Cutadapt
### for trimming primers and adapters of paired-end Illumina rads.
#Install and call required libraries. dada2, ggplot are already installed on Ollie.
#Install BiocManager for maintaining Bioconductor Packages' compatibility

print(Sys.time()) #to check runtime of the script

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#call dada2
if (!is.loaded("dada2")) {
    library(dada2)
    }
#call ShortRead
if (!is.loaded("ShortRead")) {
    library(ShortRead)
    }
#call Biostrings
if (!is.loaded("Biostrings")) {
    library(Biostrings)
    }
#Install phyloseq; it will be called later.
BiocManager::install("phyloseq")

dir <- "/work/ollie/byapan/playground/bact_1088/Renamed" # path to the directory containing the raw sequences in fastq.gz format
#Ollie works much faster if working directory is on "/work/ollie" directory.
setwd(dir)
getwd()

list.files(dir) #all raw readings (forward and reverse) should be listed.

## Create vectors for forward and reverse readings of the samples.
# If  forward and reverse fastq filenames have a name format:
# SAMPLENAME_R1_001.fastq.gz and SAMPLENAME_R2_001.fastq.gz
# If they do not; they could be arranged by file_naming_210401.sh
# Then continue as below:

# forward readings:
R1_list <- sort(list.files(dir, pattern="_R1.fastq.gz", full.names = TRUE))
# reverse readings:
R2_list <- sort(list.files(dir, pattern="_R2.fastq.gz", full.names = TRUE))


## Define forward and reverse primers.
# In order to learn primer sets used by DSET group you could check the txt file:
# primers_list_16S_tag_seq_HABITAT.txt

# In Deep-sea Mining Project primers and target regions are listed below.
# DO NOT FORGET comment out the redundant primers.

FWD <- "CCTACGGGNGGCWGCAG"  ## 341F Bacteria V3-V4 forward primer sequence.
#FWD <- "GYGCASCAGKCGMGAAW"  ## 349F Archaea V3-V5 forward primer sequence.

REV <- "GACTACHVGGGTATCTAATCC"  ## 785R Bacteria V3-V4 reverse primer sequence.
#REV <- "GTGCTCCCCCGCCAATTCCT"  ## 915R Archaea V3-V5 reverse primer sequence.

## Arrange orientation of primer sequences.
allOrients <- function(primer) {
  #Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)

REV.orients <- allOrients(REV)

FWD.orients

REV.orients

#Printing checkpoints may make it easier to see failure point if you run script on ollie by queing.
print("checkpoint0")
print(Sys.time())

### Cutadapt part
# After arranging primers it is time to make arrangements for using Cutadapt.
cutadapt <- "/global/AWIsoft/bio/cutadapt/3.2/bin/cutadapt"  # This is the path for Cutadapt on AWI servers
system2(cutadapt, args = "--version") # Run shell commands from R

dir.cut <- file.path(dir, "cutadapt")
if(!dir.exists(dir.cut)) dir.create(dir.cut) #Create a subfolder for trimmed sequences.

R1_list.cut <- file.path(dir.cut, basename(R1_list))
R2_list.cut <- file.path(dir.cut, basename(R2_list))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

#Trim FWD and the reverse-complement of REV off of R1 (forward reads)
# Flags define options and variables for cutadapt function.
R1.flags <- paste("-g", FWD, "-a", REV.RC)

#Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
# Flags define options and variables for cutadapt function.
R2.flags <- paste("-G", REV, "-A", FWD.RC)

#Run Cutadapt
for(i in seq_along(R1_list)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", R1_list.cut[i], "-p", R2_list.cut[i], # output files
                             R1_list[i], R2_list[i])) # input files
}

#rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
#      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
#      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
#      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))


print("checkpoint 1")
print(Sys.time())

R1_list.cut <- sort(list.files(dir.cut, pattern = "_R1.fastq.gz", full.names = TRUE))
R2_list.cut <- sort(list.files(dir.cut, pattern = "_R2.fastq.gz", full.names = TRUE))

save.image(file="reseq_amplicon.RData")

#plotQualityProfile(R1_list.cut[1:2])

pdf(file="quality_graph_r2.pdf")

for (i in length(R2_list.cut))
{
  plotQualityProfile(R1_list.cut[i])
}
dev.off()


#Forward and reverse fastq filenames have the format respectively:
#SAMPLE_NAME_1_R1_001.fastq.gz
#SAMPLE_NAME_1_R2_001.fastq.gz
#Extract sample names
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]

sample.names <- unname(sapply(R1_list.cut, get.sample.name))


R1_filtered <- file.path(dir.cut, "filtered", basename(R1_list.cut))
R2_filtered <- file.path(dir.cut, "filtered", basename(R2_list.cut))


#FILTER and TRIM. In this step sequences are trimmed from the end side.
#Parameters for trimming are set manually according to quality of reads.
#The quality should be checked by checking fastqc results.
#FIGARO package could help on parameter decision: https://github.com/Zymo-Research/figaro#figaro
#In order to speed up downstream processes maxEE (max number of expected errors)
#could be tightened; if the number of passing reads is too low (could be checked by parameter "out" after running)
#maxEE parameter could be increased
#For more information on error rates: https://academic.oup.com/bioinformatics/article/31/21/3476/194979
out <- filterAndTrim(R1_list.cut, R1_filtered, R2_list.cut, R2_filtered, truncLen=c(260,200),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)


# Error rate calculation is one of the most time consuming part of the run
# If you work on a HPC, multithread=TRUE makes job distributed among cores and makes it faster.
# You should be sure you have configured high number of cores and high memory
# in slurm script amplicon_analysis.sl
R1_error <- learnErrors(R1_filtered, multithread=TRUE)
R2_error <- learnErrors(R2_filtered, multithread=TRUE)

#Error rates could be plotted.
#error_plot_R1 <- plotErrors(R1_error, nominalQ=TRUE); error_plot_R1
#error_plot_R2 <- plotErrors(R2_error, nominalQ=TRUE); error_plot_R2

save.image(file="reseq_amplicon.RData")
print("checkpoint 2")
print(Sys.time())

# DEREPLICATION. To decrease memory and CPU demand for large datasets and reduce time.
# In dereplication, identical sequences counted as one unique sequence and
# according to their number of those sequence a "abundance" value is given them.
# Further steps are done by using unique sequences and abundance value,
# By that way computational resources and time are saved.
# For further information:
# see https://rdrr.io/bioc/dada2/man/derepFastq.html
# see https://benjjneb.github.io/dada2/tutorial_1_8.html
R1_dereplicated <- derepFastq(R1_filtered, verbose=TRUE)
R2_dereplicated <- derepFastq(R2_filtered, verbose=TRUE)

# Use the same sample names in the dereplicated data:
names(R1_dereplicated) <- sample.names
names(R2_dereplicated) <- sample.names


# DADA; starring of the show. Sample inference algorithm:
# dada infers sequences and resolves differences as fine as 1 nucleotide difference.
# For further information (beware paywall): https://www.nature.com/articles/nmeth.3869#methods
# pool is a very important parameter for determination of ASVs and time/resources management.
# pool=FALSE - by default samples are processed individually
# pool=TRUE - pooling to increase sensitivity (https://benjjneb.github.io/dada2/pool.html#pooling-for-sample-inference)
# pool="pseudo" - pseudo-pooled where samples are still processed individually (https://benjjneb.github.io/dada2/pseudo.html#Pseudo-pooling)

# Pooling might be important in cases such as samples from different environments are studied together.
# A sequence which is rare in an environment could be abundant in others; however, if pooling is not applied
# that rare taxa might be missed.
# if pooling is allowed by pool=TRUE it might take so long time to process high number of reads.
# pseudo pooling option could be the best option; time could be saved while keeping inference relatively sensitive

R1_dada <- dada(R1_dereplicated, err=R1_error, multithread=TRUE, pool=TRUE)
R2_dada <- dada(R2_dereplicated, err=R2_error, multithread=TRUE, pool=TRUE)

R1_dada[[1]] # To access sample "1"
R2_dada[[1]]

save.image(file="reseq_amplicon.RData")
print("checkpoint 3")
print(Sys.time())

# MERGING. Merging of matching forward and reverse reads to obtain sequence of
# the region of interest.
# In case of low number of merging reads; the parameters
# maxMismatch could be increased, minOverlap could be decreased (default is 20)
# justConcatenate could be witched to TRUE.
R1R2_merged <- mergePairs(R1_dada, R1_dereplicated, R2_dada, R2_dereplicated, justConcatenate=FALSE, verbose=TRUE)

# Inspect the merged data.frame from the first sample
head(R1R2_merged[[1]])

# Construct the amplicon sequence variant table (ASV)
# This is analogous to an OTU table but with resolution up to single-nucleotide level.
ASV_table <- makeSequenceTable(R1R2_merged)
dim(ASV_table)
# Inspect distribution of sequence lengths
table(nchar(getSequences(ASV_table)))
#hist(nchar(getSequences(ASV_table)))

# Remove chimaeras
ASV_nochim <- removeBimeraDenovo(ASV_table, method="consensus", multithread=TRUE, verbose=TRUE)
dim(ASV_nochim)
# Proportion of non-chimaeric sequences:
sum(ASV_nochim)/sum(ASV_table)

#Read counts throughout the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(R1_dada, getN), sapply(R2_dada, getN), sapply(R1R2_merged, getN), rowSums(ASV_nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
#track
write.csv(track,"track_reads.csv")
saveRDS(ASV_nochim, "/work/ollie/byapan/playground/asv_bact_1088.rds")
save.image(file="reseq_amplicon.RData")
print("checkpoint 4")
print(Sys.time())

# For the studies which comprised samples from different sequencing runs
# it is better to recall ASV_tables saved as rds and continue the steps after here
# by merging them.

# TAXONOMIC ANNOTATION
# The database in the correct format can be found in the dada2 website.
# You should download the database to a directory you chhose and call it from there.
unite.ref <- "/work/ollie/byapan/playground/silva/silva_nr99_v138.1_train_set.fa.gz"
ASV_taxonomy <- assignTaxonomy(ASV_nochim,unite.ref, multithread=TRUE, verbose=TRUE)
# Add species:
unite.ref_sp <-"/work/ollie/byapan/playground/silva/silva_species_assignment_v138.1.fa.gz"
ASV_sp <- addSpecies(ASV_taxonomy, unite.ref_sp, verbose=TRUE)
saveRDS(ASV_sp, "/work/ollie/byapan/playground/asv_sp_bact_1088.rds")
#ASV_taxonomy_check <- ASV_taxonomy
#rownames(ASV_taxonomy_check) <- NULL


save.image(file="bacteria_amplicon.RData")
print("checkpoint 5")
print(Sys.time())
q()
##############################################################################################

# Analyse and plot results with phyloseq.

# Import results to phyloseq:
#library(phyloseq); packageVersion("phyloseq")
#library(ggplot2); packageVersion("ggplot2")
#theme_set(theme_bw())

# Extract the sample and ASV names:
#samples.out <- rownames(ASV_nochim)
#ASVs <- colnames(ASV_nochim)
# ASVs ID table:
#ASVs_ID <- cbind(ASVs, paste("asv", c(1:ncol(ASV_nochim)), sep=""))

# rename the ASV to asv#:
#colnames(ASV_nochim) <- paste("asv", c(1:ncol(ASV_nochim)), sep="")
#rownames(ASV_taxonomy) <- paste("asv", c(1:nrow(ASV_taxonomy)), sep="")
#ASV_taxonomy[is.na(ASV_taxonomy[,1])] <- "Unclassified" # Replace empty taxons (domain/kingdom level) with "Unclassified".

# Add sample names:
#head (samples.out)
#samples.out3<- #cbind(samples.out,c("nod111","nod112","nod113","nod117","nod118","nod119","nod173","nod75","no76","nod77","nod81","nod83","dn1","dn2","w49","w151","w159","120","121","127","128","134","135","2","22","23","43","44","50","51","57","58","78","79","80","85",,"86","9","92","93","dis05","dis06","dis07","dis08","dis14","dis15"),c("BGR_ref","BGR_ref","BGR_ref","BGR_ref","BGR_ref","BGR_ref","BGR_ref","GSR_ref","GSR_ref","GSR_ref","GSR_ref","GSR_ref","DISCOL","DISCOL","GSR_tri","GSR_ref","BGR_ref","BGR_r#ef","BGR_ref","BGR_ref","BGR_ref","BGR_ref","BGR_ref","GSR_tri","GSR_tri","BGR_tri","GSR_tri","GSR_tri","BGR_tri","GSR_tri","GSR_tri","DIS","DIS","DIS","DIS","DIS","DIS"))
#colnames(samples.out3) <- c("ID", "Sample","Location")
#rownames(samples.out3) <- samples.out3[,1] # Row names are samples IDs.
#samples.out3 <- as.data.frame(samples.out3)

#OTU_phyloseq3 <- otu_table(ASV_nochim, taxa_are_rows = FALSE)
#OTU_phyloseq3<- filter_taxa(OTU_phyloseq3, function (x) {sum(x > 0) > 1}, prune=TRUE) #Remove singletons
#SAMPLE_phyloseq3 <- sample_data(samples.out3)
#TAX_phyloseq3 <- tax_table(ASV_taxonomy)
#TAX_phyloseq3<- subset_taxa(TAX_phyloseq3, Kingdom=="Archaea") #Remove archaeal reads
#dada2_phyloseq4<- phyloseq(OTU_phyloseq3, TAX_phyloseq3, SAMPLE_phyloseq3) # Create a phyloseq object.

#phyloseq-class experiment-level object w/ singletons
#otu_table()   OTU Table:         [ 34916 taxa and 30 samples ]
#sample_data() Sample Data:       [ 30 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 34916 taxa by 6 taxonomic ranks ]

#phyloseq-class experiment-level object wtihout singletons
#otu_table()   OTU Table:         [ 15682 taxa and 30 samples ]
#sample_data() Sample Data:       [ 30 samples by 3 sample variables ]
#tax_table()   Taxonomy Table:    [ 15682 taxa by 6 taxonomic ranks ]


#phy = dada2_phyloseq4

#NameTax <- function(x, ind){
#  if(is.na(x[ind])){
#    x[ind] <- x[ind]
#  } else {
#    if(ind==1){x[ind] <- paste("d", x[ind], sep="__")} else{                                  # Domain
#      if(ind==2){x[ind] <- paste("p", x[ind], sep="__")} else{                                # Phylum
#        if(ind==3){x[ind] <- paste("c", x[ind], sep="__")} else{                              # Class
#          if(ind==4){x[ind] <- paste("o", x[ind], sep="__")} else{                            # Order
#            if(ind==5){x[ind] <- paste("f", x[ind], sep="__")} else{                          # Family
#              if(ind==6){x[ind] <- paste("g", paste(x[ind-1], x[ind], sep="_"), sep="__")}    # Genus
#              }
#            }
#          }
#        }
#      }
#    }
#  }
#}

#tax.tab <- data.frame(tax_table(phy))

#for (i in 1:7) {
#  tax_table(phy)[,i] <- apply(tax.tab, 1, NameTax, ind=i)
#}

#ModifyTax <- function(x,ind){
#  #   xth row in the dataframe
#  #   ind taxonomy level to change
#  if(is.na(x[ind])){
#    nonNa <- which(!is.na(x[-ind])) # which taxa are not NA excepting the one we're interested in.
#    maxNonNa <- max(nonNa)
#    x[ind] <- x[maxNonNa]
#  }else{x[ind] <- x[ind]}
#}

#for (i in 1:7) {
#  tax_table(phy)[,i] <- apply(tax.tab,1,ModifyTax,ind=i)
#}

#phy_rare <- phy

#wh0 = genefilter_sample(phy_samples, filterfun_sample(function(x) x > 5), A=0.5*nsamples(phy_samples))
#GP1 = prune_taxa(wh0, phy_samples)
#GP2 = transform_sample_counts(GP1, function(x) 100 * x/sum(x))
#GP20 <- transform_sample_counts(phy_rare, function(x) sqrt(x / sum(x)))
#GP30 <- transform_sample_counts(phy_rare, function(x) if (x>=1) {x=1} else {x=0} )

#GP2 = transform_sample_counts(GP1, decostand(otu_table(GP1),"hellinger"))

#GP.ord_rare <- ordinate(phy_rare, "NMDS", "Jaccard")
#Run 17 stress 0.06218149
#... Procrustes: rmse 8.369671e-05  max resid 0.0003976089
#... Similar to previous best
#Run 18 stress 0.06218146
#... Procrustes: rmse 4.971322e-05  max resid 0.0002354726
#... Similar to previous best
#Run 19 stress 0.230557
#Run 20 stress 0.0645415
#*** Solution reached

#p11 = plot_ordination(GP2, GP.ord2, type="taxa", color="Phylum", title="taxa")
#p20 = plot_ordination(GP20, GP.ord_rare, type="samples", color="location")

#vegan_otu <- function(GP2) {
#    OTU <- otu_table(GP2)
#    if (taxa_are_rows(OTU)) {
#        OTU <- t(OTU)
#    }
#    return(as(OTU, "matrix"))
#}

#metadata<- as(sample_data(phy_rare), "data.frame")
#sim_GP_rare <- anosim(otu_table(phy_rare), metadata$location, permutations = 999, distance = "jaccard", strata = NULL)

#p42 <- plot_richness(phy_rare, x="location", measures=c("Observed", "Shannon"),color="location")+
#scale_x_discrete(limits=c("DIS","BGR_ref","BGR_tri","GSR_ref","GSR_tri"))

#Call:
#anosim(x = otu_table(GP20), grouping = metadata$location, permutations = 999,      distance = "bray", strata = NULL)
#Dissimilarity: bray

#ANOSIM statistic R: 0.8224
#      Significance: 0.001

#Permutation: free
#Number of permutations: 999

#Upper quantiles of permutations (null model):
#   90%    95%  97.5%    99%
#0.0764 0.1067 0.1296 0.1481

#Dissimilarity ranks between and within classes:
#        0%    25%   50%    75% 100%   N
#Between  4 159.75 254.5 345.25  435 360
#BGR_ref 47  72.50  92.0 101.50  118  15
#BGR_tri 75  93.50 139.0 171.50  256  15
#DIS      1  18.00  29.0  36.00   58  15
#GSR_ref  2   8.50  32.0  48.50   67  15
#GSR_tri 16  32.50  56.0  87.50  154  15


#TopNOTUs <- names(sort(taxa_sums(GP2), TRUE)[1:10])
#ent10   <- prune_taxa(TopNOTUs, GP2)
#ent10 = transform_sample_counts(ent10, function(x) 100 * x/sum(x))
#sample_names(ent10) <- c("neg0","neg1","neg2","1","120","121","127","128","134","135","2","22","23","43","44","50","51","57","58","78","79","8","85","86","9","92","93","d05","d06","d07","d08","d14","d15")
#p35=plot_bar(ent10, fill="Order", x="sample_Sample")
#print("checkpoint6")
#save.image(file="archaea_0.RData")
