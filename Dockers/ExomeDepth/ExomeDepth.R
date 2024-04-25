# Script    : ExomeDepth.R
# Objective : To call CNVs using ExomeDepth
# Written by: egustavsson

# docker run --rm -it -v /mnt:/mnt -v /cmg1scratch:/cmg1scratch alesmaver/exome_depth R

library(optparse)
library(ExomeDepth)
library(GenomicRanges)
library(tidyverse)
library(rtracklayer)

# Function Definitions -----------------------------------------------------
createReferenceSet <- function(targets, baseline_samples) {
  output_name <- gsub("([0-z]*).(.*)", "\\1", basename(targets))
  # Check if targets are provided; if not, generate exons.hg19 object
  if (missing(targets) || is.null(targets)) {
    data("exons.hg19")
    targets <- exons.hg19
  } else {
    targets <- read.table(targets, header = FALSE, col.names = c("chrom", "start", "end", "info"))
  }
  
  Reference_Counts <- getBamCounts(bed.frame = targets,
                         bam.files = baseline_samples,
                         include.chr = FALSE) 

  save(Reference_Counts, file = file.path(paste0(output_name, ".ExomeDepth.RData")))
}

callCNVs <- function(targets, annotation, test_sample, baseline_samples, reference_counts_file) {

  # Check if targets are provided; if not, generate exons.hg19 object
  if (missing(targets) || is.null(targets)) {
    data("exons.hg19")
    targets <- exons.hg19
  } else {
    targets <- read.table(targets, header = FALSE, col.names = c("chrom", "start", "end", "info"))
  }
  
  # Check if annotations are provided; if not, generate genes.hg19 object
  if (missing(annotation) || is.null(annotation)) {
    data("genes.hg19")
    annotation <- genes.hg19 %>%
  dplyr::rename(gene_name = name) %>%
  mutate(chromosome = paste0("chr", chromosome)) %>%
  GRanges()                                                 
  } else {
    annotation <- rtracklayer::import(opt$annotation) %>% .[.$type == "gene"] %>% unique() # This needs to have "chr" within seqnames
  }

  if (exits(reference_counts_file)) {
    load(reference_counts_file)
  } else {
    error("Reference counts file not defined. Please generate a reference file before running this command.")
  }

  Reference_Counts <- getBamCounts(bed.frame = targets,
                         bam.files = baseline_samples,
                         include.chr = FALSE) 


  Test_Counts <- getBamCounts(bed.frame = targets,
                         bam.files = c(test_sample, baseline_samples),
                         include.chr = TRUE) 
                         
                    cbind()     %>%
    setNames(gsub("^X(\\d+)", "\\1", names(.))) # Remove 'X' from column names starting with a number. 
  
  Counts.df <- as.data.frame(Counts)
  
  my.reference.set <- as.matrix(Counts.df[, basename(baseline_samples)])
  
  my.test <- Counts[, basename(test_sample)]
    
  my.choice <- select.reference.set(test.counts = my.test,
                                    reference.counts = my.reference.set,
                                    bin.length = (Counts.df$end - Counts.df$start) / 1000,
                                    n.bins.reduced = 10000)
    
  my.matrix <- as.matrix(Counts.df[, my.choice$reference.choice, drop = FALSE])
  my.reference.selected <- apply(X = my.matrix, MAR = 1, FUN = sum)
    
  all.exons <- new('ExomeDepth', 
                   test = my.test,
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')

  all.exons <- CallCNVs(x = all.exons,
                        transition.probability = 10^-4,
                        chromosome = Counts$chromosome,
                        start = Counts$start,
                        end = Counts$end,
                        name = Counts$exon)

  # sort by BF value and annotate
  CNV_calls <- all.exons@CNV.calls %>% arrange(desc(BF)) %>% GRanges()

  # Find overlaps using "any" method to handle partial overlaps
overlap_hits <- findOverlaps(CNV_calls, annotation, type = "any")

# Combine gene names for overlapping ranges
combine_gene_names <- function(gene_names) {
  return(paste(unique(gene_names), collapse = ","))
}

# Initialize a list to store gene names for each CNV_calls entry
gene_names_list <- vector("list", length(CNV_calls))

# Loop through the overlaps and update the gene_names_list
for (i in seq_along(gene_names_list)) {
  # Find all overlaps for the current CNV_calls entry
  current_overlap_indices <- subjectHits(overlap_hits)[queryHits(overlap_hits) == i]
  
  if (length(current_overlap_indices) > 0) {
    # Extract the gene names from annotation for the current overlaps
    overlapping_genes <- mcols(annotation[current_overlap_indices])$gene_name
    
    # Combine gene names with commas and store them in the gene_names_list
    gene_names_list[[i]] <- combine_gene_names(overlapping_genes)
  }
}

# Update the "CNV_calls" object with the combined gene names, preserving NA values
CNV_calls$gene_name <- CNV_calls$gene_name <- unlist(lapply(gene_names_list, function(x) if (is.null(x)) NA else x))
  
  # Generate the output filename based on the test sample name
  sample_name <- gsub("\\.bam$", "", basename(test_sample))
  output_file <- file.path(output_directory, paste0(sample_name, "_CNV.csv"))
    
  write.csv(file = output_file,
            x = CNV_calls,
            row.names = FALSE)
  
  # Print completion message for the test sample
  cat("Analysis completed for", sample_name, "\n")
}

# Parse command-line options
option_list <- list(
  make_option("--test-samples", dest="test_samples", type="character"),
  make_option("--baseline-samples", dest="baseline_samples", type="character"),
  make_option("--targets", dest="targets", type="character"),
  make_option("--annotation", dest="annotation", type="character"),
  make_option("--output-directory", dest="output_directory", type="character"),
  make_option("--create-reference-set", dest="create_reference_set", type="boolean", default=FALSE) 
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Redirect output to a log file to suppress warnings and messages
log_file <- file.path(opt$output_directory, "exomedepth_log.txt")
sink(log_file, append = FALSE)

# Read test samples from TSV
test_samples <- read_tsv(opt$test_samples, col_names = "test_sample_path", show_col_types = F)

# Read baseline samples from TSV
baseline_samples <- read_tsv(opt$baseline_samples, col_names = "baseline_sample_path", show_col_types = F)

# Run the analysis for each test sample
if(create_reference_set){
    createReferenceSet(
        targets = opt$targets,
        baseline_samples = opt$baseline_sample_path 
    )
} else {
    for (test_sample_path in test_samples$test_sample_path) {
        callCNVs(
            targets = opt$targets,
            annotation = opt$annotation,
            test_sample = test_sample_path,
            baseline_samples = baseline_samples$baseline_sample_path,
            output_directory = opt$output_directory  # Updated argument name
            )
    }
}

# Close the sink to restore the standard output
sink()

#
targets = '/cmg1scratch/cromwell/FastqToVCFPipelineDatasets/ENRICHMENTS//nexteraExome_nochrYM.bed'
test_sample = '/mnt/dataSeq/DATA_REPOSITORY/BAMS_HG19/PX13069.marked.bam'
baseline_samples = c('/mnt/dataSeq/DATA_REPOSITORY/BAMS_HG19/PX13871.marked.bam', 
                     '/mnt/dataSeq/DATA_REPOSITORY/BAMS_HG19/PX13872.marked.bam', 
                     '/mnt/dataSeq/DATA_REPOSITORY/BAMS_HG19/PX14310.marked.bam')
annotation = NULL

opt=c()
opt$targets = targets
opt$test_samples = test_sample
opt$baseline_sample_path = baseline_samples