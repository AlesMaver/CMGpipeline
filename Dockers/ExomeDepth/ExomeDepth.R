# Script    : ExomeDepth.R
# Objective : To call CNVs using ExomeDepth
# Written by: egustavsson

# docker run --rm -it -v /mnt:/mnt -v /cmg1scratch:/cmg1scratch alesmaver/exome_depth R

library(optparse)
library(ExomeDepth)
library(GenomicRanges)
library(tidyverse)
library(rtracklayer)
library(doParallel)
library(RcppRoll)

# Function Definitions -----------------------------------------------------
createReferenceSet <- function(targets, baseline_samples) { 
  # Check if targets are provided; if not, generate exons.hg19 object
  if (missing(targets) || is.null(targets)) {
    data("exons.hg19")
    targets <- exons.hg19
  } else {
    targets <- read.table(targets, header = FALSE, sep="\t")
    if (ncol(targets) == 3) {
      colnames(targets) <- c("chrom", "start", "end")
    } 
    if (ncol(targets) >= 4) {
      colnames(targets)[1:4] <- c("chrom", "start", "end", "info")
    } 
  }
  
  Reference_Counts <- getBamCounts(bed.frame = targets,
                         bam.files = baseline_samples,
                         include.chr = FALSE) 

  return(Reference_Counts)
}

# Parallel creation of counts files
createReferenceCounts = function(samples, targets, working_directory, ncores=10) {
# Create a cluster with the specified number of cores
  registerDoParallel(cores = ncores)

  # Use parallel::foreach to parallelize the for loop
  foreach(sample = samples) %dopar% {
    counts <- createReferenceSet(
      targets = opt$targets,
      baseline_samples = sample
    )
    sample_name = gsub("([A-z0-9]*)\\..*", "\\1", basename(sample))
    write.table(counts, file = paste0(working_directory, sample_name, "_ExomeDepth_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

get_reference_ratio_sd <- function(reference_logratios_files){
  reference_logratios_files <- read.table(opt$reference_logratios_file_list , col.names = "baseline_sample_path")[,1]
  reference_logratios<-list()
  for (i in seq_along(reference_logratios_files)) {
      ratios <- read.table(reference_logratios_files[i], header = TRUE, sep = " ")[,4]
      reference_logratios[[basename(reference_logratios_files[i])]] <- ratios
  }

  lengths <- lapply(reference_logratios, function(x) length(x))
  acceptable_length <- max(unlist(lengths))
  reference_logratios[lengths != acceptable_length] <- NULL

  reference_logratios_df <- as.data.frame(reference_logratios)
  my.reference.selected.sd <- apply(X = reference_logratios_df, MAR = 1, FUN = sd)
  return(my.reference.selected.sd)
}

# Function to convert, sort, gzip, and tabix index a file
compress_and_index <- function(file_path) {
  # Convert space-delimited file to tab-delimited in-place
  system(paste("sed -i 's/ /\\t/g' ", file_path))
  cat("File converted to tab-delimited:", file_path, "\n")
  
  # Sort the file
  gz_file <- paste0(file_path, ".gz")
  system(paste("sort -k1,1 -k2,2n", file_path, " | bgzip > ", gz_file))
  cat("File sorted:", gz_file, "\n")

  # Check if index file already exists
  index_file <- paste0(gz_file, ".tbi")
  if (file.exists(index_file)) {
    cat("Index file already exists:", index_file, "\n")
  } else {
    # Index the compressed file using tabix
    system(paste("tabix -p bed", gz_file))
    cat("Index created for:", gz_file, "\n")
  }
}

# Function to process WIG files and call compress_and_index
process_wig_files <- function(working_directory, sample_name) {
  # List of WIG files
  wig_files <- list.files(working_directory, pattern = paste0(sample_name, "_ExomeDepth", "_.*\\.wig$"), full.names = TRUE)
  
  # Iterate through each WIG file
  for (wig_file in wig_files) {
    # Read the WIG file
    wig <- read.table(wig_file, header = FALSE, quote = "", sep = "\t", stringsAsFactors = FALSE)
    
    # Compress and index the WIG file
    compress_and_index(wig_file)
  }
}

# CNV calling function
callCNVs <- function(targets, annotation, test_counts_file, reference_counts_files, probes_variance = NULL, probes_variance_cutoff=1) {
  sample_name = gsub("_(.*)", "", basename(test_counts_file))

  ########################
  # TARGETS
  ########################
  # Check if targets are provided; if not, generate exons.hg19 object
  print("Reading targets file")
  if (missing(targets) || is.null(targets)) {
    data("exons.hg19")
    targets_df <- exons.hg19
  } else {
    targets_df <- read.table(targets, header = FALSE, sep="\t")

    targets_ncol = ncol(targets_df)

    if (ncol(targets_df) == 3) {    
      targets_df$info <- "NA"
    } 
    if (ncol(targets_df) > 4) {    
      targets_df <- targets_df[,1:4]
      targets_ncol = ncol(targets_df)
    } 
    
  }
  colnames(targets_df)[1:4] <- c("chrom", "start", "end", "info")

  ########################
  # ANNOTATION
  ########################
  print("Reading annotation file")
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

  ########################
  # REFERENCE COUNTS
  ########################
  print("Reading reference counts files")
  for (i in seq_along(reference_counts_files)) {
    print(paste0("Reading reference counts file ", reference_counts_files[i]))
    if (i == 1) {
      df <- read.table(reference_counts_files[i], header = TRUE, sep = "\t")
    } else {
      df_current <- read.table(reference_counts_files[i], header = TRUE, sep = "\t")
      
      df_current[1:targets_ncol]<-NULL # Skip first four columns that repeat in each file
      print(head(df_current))

      # Prevent errors due to merging of data frames with unequeal row numbers
      if(nrow(df) == nrow(df_current)){
      	df <- cbind(df, df_current)
	    }

    }
  }
  reference_counts = df
  ###

  ########################
  # TEST COUNTS
  ########################
  print("Reading test counts file")
  test_counts <- read.table(test_counts_file, header = TRUE, sep = "\t")

  ########################
  # REFERENCE AND TEST VECTORS
  ########################
  print("Creating reference and test vectors")
  my.reference.set <- as.matrix(reference_counts[,(targets_ncol+1):ncol(reference_counts)])
  head(my.reference.set)
  # Remove any samples from the references if they are in the test set
  my.reference.set <- my.reference.set[, !colnames(my.reference.set) %in% colnames(test_counts)]
  
  my.test <- test_counts[,(targets_ncol+1)]

  ########################
  # CHOOSE THE APPROPRIATE SAMPLES FOR THE REFERENCE
  ########################
  print("Choosing the appropriate samples for the reference")
  my.choice <- select.reference.set(test.counts = my.test,
                                    reference.counts = my.reference.set,
                                    bin.length = (df$end - df$start) / 1000,
                                    n.bins.reduced = 10000)
    
  my.matrix <- as.matrix(df[, my.choice$reference.choice, drop = FALSE])
  my.reference.selected <- apply(X = my.matrix, MAR = 1, FUN = sum)
  my.reference.selected.variance <- apply(X = my.matrix, MAR = 1, FUN = var)
  
  ########################
  # PERFORM EXOME DEPTH ANALYSIS
  ########################
  all.exons <- new('ExomeDepth', 
                   test = my.test,
                   reference = my.reference.selected,
                   formula = 'cbind(test, reference) ~ 1')

  if (targets_ncol==3){
    all.exons <- CallCNVs(x = all.exons,
                          transition.probability = 10^-4,
                          chromosome = reference_counts$chromosome,
                          start = reference_counts$start,
                          end = reference_counts$end,
                          name = rep("NA", nrow(reference_counts))
                          )
    } else {
    all.exons <- CallCNVs(x = all.exons,
                    transition.probability = 10^-4,
                    chromosome = reference_counts$chromosome,
                    start = reference_counts$start,
                    end = reference_counts$end,
                    name = reference_counts[,4]
                    )
    }

  if (nrow(all.exons@CNV.calls) == 0) {
    all.exons@CNV.calls <- data.frame(chromosome = c("chr1"),
                      start = c(1),
                      end = c(2),
                      BF = c(1),
                      type = c("DUP")
                      )
  }
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
  
  all.exons@CNV.calls = all.exons@CNV.calls[ order ( all.exons@CNV.calls$BF, decreasing = TRUE),]

  # Order columns in the desired manner
  CNV_table <- as.data.frame(CNV_calls)
  col_names_primary <- c("seqnames", "start", "end", "type", "BF", "gene_name")
  col_names_other <- setdiff(colnames(CNV_table), col_names_primary)
  CNV_table <- CNV_table[, c(col_names_primary, col_names_other)]
  CNV_table$type <- gsub("deletion", "DEL", CNV_table$type)
  CNV_table$type <- gsub("duplication", "DUP", CNV_table$type)
  # Write the CNV calls to a CSV file
  output_file_csv <- file.path(working_directory, paste0(sample_name, "_ExomeDepth_CNV_annotSV.bed"))
  write.table(file = output_file_csv,
            x = CNV_table,
            quote = FALSE,
            row.names = FALSE,
            sep = "\t")

  # Write the CNV calls to a BED file
  output_file_bed <- file.path(working_directory, paste0(sample_name, "_ExomeDepth_CNV.bed"))
  CNV_bed = as.data.frame(CNV_calls)[,c("seqnames", "start", "end", "type", "gene_name")]
  CNV_bed$type <- gsub("deletion", "DEL", CNV_bed$type)
  CNV_bed$type <- gsub("duplication", "DUP", CNV_bed$type)
  write.table(file = output_file_bed,
              x = CNV_bed,
              quote = FALSE,
              col.names = FALSE,
              row.names = FALSE,
              sep = "\t")
  
  # WIG FILE - TEST SAMPLE
  all.exons_freq = all.exons@test/(all.exons@test + all.exons@reference)
  all.exons_oe = all.exons_freq/all.exons@expected
  wig = cbind(reference_counts[,1:3], all.exons_oe)
  # wig = cbind(reference_counts[,1:3], all.exons@likelihood)
  # wig[,4] <- format(wig[,4], scientific = FALSE)
  # Save to a wig file without quotes and without column names
  write.table(wig, file = file.path(working_directory, paste0(sample_name, "_ExomeDepth", "_ratios_all.wig")), quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  if (!is.null(probes_variance)) {
      wig_clean <- wig[probes_variance<probes_variance_cutoff,]
  } else {
      wig_clean <- wig
  }
  wig_clean <- subset(wig_clean, !is.na(wig_clean[,4]))
  write.table(wig_clean, file = file.path(working_directory, paste0(sample_name, "_ExomeDepth", "_ratios_clean.wig")), quote = FALSE, col.names = FALSE, row.names = FALSE)

  tryCatch({
      wig[which(is.nan(wig[,4])),]<-NULL
  }, error = function(e) {
    # Handle the exception here
    print(paste("An error occurred:", e$message))
  })
  write.table(wig, file = file.path(working_directory, paste0(sample_name, "_ExomeDepth", "_ratios_nomissing.wig")), quote = FALSE, col.names = FALSE, row.names = FALSE)

  rolling_probes = 10
  wig_rolling = cbind(reference_counts[,1:3], all.exons_oe)
  wig_rolling <- subset(wig_rolling, !is.nan(wig_rolling[,4]))
  # We want to roll over cleaned probes to enhance the signal
  if (!is.null(probes_variance)) {
    wig_rolling <- wig_rolling[probes_variance<probes_variance_cutoff,]
  }
  wig_rolling[,4] <- roll_mean(wig_rolling[,4], n = rolling_probes, align = "center", fill = NA)
  wig_rolling <- subset(wig_rolling, !is.na(wig_rolling[,4]))
  wig_rolling[,4] <- format(wig_rolling[,4], digits = 3, scientific=FALSE)
  write.table(wig_rolling, file = file.path(working_directory, paste0(sample_name, "_ExomeDepth", "_rolling_ratios.wig")), quote = FALSE, col.names = FALSE, row.names = FALSE)

  # # WIG FILE - REFERENCE SAMPLES
  # all.exons_freq = all.exons@test/(all.exons@test + all.exons@reference)
  # all.exons_oe = all.exons_freq/all.exons@expected
  # wig = cbind(df[,1:3], all.exons_oe)
  # wig = wig[-which(is.nan(wig[,4])),]
  # # Save to a wig file without quotes and without column names
  # write.table(wig, file = file.path(working_directory, paste0("ExomeDepth", "_CNV.wig")), quote = FALSE, col.names = FALSE, row.names = FALSE)

  process_wig_files(working_directory=working_directory, sample_name=sample_name)

  # Print completion message for the test sample
  cat("Analysis completed for", sample_name, "\n")
}

# Parse command-line options
option_list <- list(
  make_option("--test-sample-bam", dest="test_sample_bam", type="character"),
  make_option("--baseline-sample-bam_list", dest="baseline_samples_bam_list", type="character", default=NULL),
  make_option("--targets", dest="targets", type="character"),
  make_option("--annotation", dest="annotation", type="character"),
  make_option("--output-directory", dest="output_directory", type="character", default="./"),
  make_option("--working-directory", dest="working_directory", type="character", default="./"),
  make_option("--test-counts-file", dest="test_counts_file", type="character", default=NULL),
  make_option("--reference-counts-file-list", dest="reference_counts_file_list", type="character", default=NULL),
  make_option("--reference-logratios-file-list", dest="reference_logratios_file_list", type="character", default=NULL)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

working_directory = opt$working_directory
print(opt)

# Redirect output to a log file to suppress warnings and messages
log_file <- file.path(opt$output_directory, "exomedepth_log.txt")
# sink(log_file, append = FALSE)

# OPTION ONE - a single BAM file is provided
# DO: Create a counts file for downstream use
if (!is.null(opt$test_sample_bam)){
  print("You have provided a single test sample bam, so we will generate counts file for downstream use. ")
  createReferenceCounts(samples=opt$test_sample_bam, opt$targets, working_directory=opt$working_directory, ncores=1)
  print("Counts file generated for the test sample.")
}

# OPTION TWO - a list of locations for multiple bam files is provided
# DO: Create multiple counts files for use in a reference set
if (!is.null(opt$baseline_samples_bam_list)){
  print("You have provided multiple baseline samples - will generate count files for all now... ")
  bam_files <- read.table(opt$baseline_samples_bam_list, col.names = "baseline_sample_path")[,1]

  createReferenceCounts(samples=bam_files, opt$targets, working_directory=opt$working_directory, ncores=30)
  print("Counts files generated for all baseline samples.")
}

# OPTION THREE - a file with test counts and a list of files with reference counts is provided
# DO: Perform Exome Depth analysis
if (!is.null(opt$test_counts_file) & !is.null(opt$reference_counts_file_list)){

  # Find probes with a high rate of variation in the reference population of ratio values
  if(!is.null(opt$reference_logratios_file_list)){
    reference_sd <- get_reference_ratio_sd(opt$reference_logratios_file_list)
  } else {
    reference_sd <- NULL
  }

  reference_counts_files <- read.table(opt$reference_counts_file_list, col.names = "baseline_sample_path")[,1]
  callCNVs(targets=opt$targets, annotation= NULL, test_counts_file=opt$test_counts_file, reference_counts_files = reference_counts_files, probes_variance = reference_sd, probes_variance_cutoff=0.2)
  print("Exome Depth analysis completed.")
}

# Close the sink to restore the standard output
# sink()

