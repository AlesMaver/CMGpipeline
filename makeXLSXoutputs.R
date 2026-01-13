#!/usr/bin/env Rscript
library("optparse")
library("openxlsx")
options(stringsAsFactors = F)

# Parse input options
option_list = list(
  make_option(c("--sample_basename"), type="character", default=NULL, help="Exome id, for example PX5000", metavar="character"),
  make_option(c("--RARE_FUNCTIONAL"), type="character", default=NULL, help="RARE_FUNCTIONAL output", metavar="character"),
  make_option(c("--HET_DOMINANT"), type="character", default=NULL, help="HET_DOMINANT output", metavar="character"),
  make_option(c("--COMPHET_RECESSIVE"), type="character", default=NULL, help="COMPHET_RECESSIVE output", metavar="character"),
  make_option(c("--HOM_RECESSIVE"), type="character", default=NULL, help="HOM_RECESSIVE output", metavar="character"),
  make_option(c("--CLINVAR_PATHOGENIC"), type="character", default=NULL, help="CLINVAR_PATHOGENIC output", metavar="character"),
  make_option(c("--CLINVAR_FILTERED"), type="character", default=NULL, help="CLINVAR_FILTERED output", metavar="character"),
  make_option(c("--CLINVAR_ALL"), type="character", default=NULL, help="CLINVAR_ALL output", metavar="character"),
  make_option(c("--PANEL_FILTERED"), type="character", default=NULL, help="PANEL_FILTERED output", metavar="character"),
  make_option(c("--PANEL_ALL"), type="character", default=NULL, help="PANEL_ALL output", metavar="character"),
  make_option(c("--MITOMAP"), type="character", default=NULL, help="MITOMAP output", metavar="character"),
  make_option(c("--XLSX_OUTPUT"), type="character", default=NULL, help="Output XLSX", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


if (is.null(opt$RARE_FUNCTIONAL) & 
    is.null(opt$HET_DOMINANT) & 
    is.null(opt$COMPHET_RECESSIVE) & 
    is.null(opt$HOM_RECESSIVE) & 
    is.null(opt$CLINVAR_PATHOGENIC) & 
    is.null(opt$CLINVAR_FILTERED) &
    is.null(opt$CLINVAR_ALL) &
    is.null(opt$PANEL_FILTERED) &
    is.null(opt$PANEL_ALL) &
    is.null(opt$MITOMAP)
  ){
  print_help(opt_parser)
  stop("At least one input must be supplied. Exiting with error...", call.=FALSE)
}

# #Debugging code - import outputs from a directory
# opt <- list()
# SAMPLE="SAMPLE"
# snpSift_outputDir<-"/home/ales/pipeline_debug/"
# opt$RARE_FUNCTIONAL<-paste0(snpSift_outputDir, SAMPLE, ".RARE_FUNCTIONAL.tab")
# opt$HET_DOMINANT<-paste0(snpSift_outputDir, SAMPLE, ".HET_DOMINANT.tab")
# opt$COMPHET_RECESSIVE<-paste0(snpSift_outputDir, SAMPLE, ".COMPHET_RECESSIVE.tab")
# opt$HOM_RECESSIVE<-paste0(snpSift_outputDir, SAMPLE, ".HOM_RECESSIVE.tab")
# opt$CLINVAR_PATHOGENIC<-paste0(snpSift_outputDir, SAMPLE, ".CLINVAR_PATHOGENIC.tab")
# opt$CLINVAR_FILTERED<-paste0(snpSift_outputDir, SAMPLE, ".CLINVAR_FILTERED.tab")
# opt$CLINVAR_ALL<-paste0(snpSift_outputDir, SAMPLE, ".CLINVAR_ALL.tab")
# opt$PANEL_FILTERED<-paste0(snpSift_outputDir, SAMPLE, ".PANEL_FILTERED.tab")
# opt$PANEL_ALL<-paste0(snpSift_outputDir, SAMPLE, ".PANEL_ALL.tab")
# opt$MITOMAP<-paste0(snpSift_outputDir, SAMPLE, ".MITOMAP.tab")


# Make a list of reports
reportList<-list()
if( !is.null(opt$RARE_FUNCTIONAL) & file.exists(opt$RARE_FUNCTIONAL) ) reportList$RARE_FUNCTIONAL <- read.table(opt$RARE_FUNCTIONAL, sep="\t", header=T, quote="", dec = ".", fill=NA)
if( !is.null(opt$HET_DOMINANT) & file.exists(opt$HET_DOMINANT) ) reportList$HET_DOMINANT <- read.table(opt$HET_DOMINANT, sep="\t", header=T, quote="", dec = ".", fill=NA)
if( !is.null(opt$COMPHET_RECESSIVE) & file.exists(opt$COMPHET_RECESSIVE) ) reportList$COMPHET_RECESSIVE <- read.table(opt$COMPHET_RECESSIVE, sep="\t", header=T, quote="", dec = ".", fill=NA)
if( !is.null(opt$HOM_RECESSIVE) & file.exists(opt$HOM_RECESSIVE) ) reportList$HOM_RECESSIVE <- read.table(opt$HOM_RECESSIVE, sep="\t", header=T, quote="", dec = ".", fill=NA)
if( !is.null(opt$CLINVAR_PATHOGENIC) & file.exists(opt$CLINVAR_PATHOGENIC) ) reportList$CLINVAR_PATHOGENIC <- read.table(opt$CLINVAR_PATHOGENIC, sep="\t", header=T, quote="", dec = ".", fill=NA)
if( !is.null(opt$CLINVAR_FILTERED) & file.exists(opt$CLINVAR_FILTERED) ) reportList$CLINVAR_FILTERED <- read.table(opt$CLINVAR_FILTERED, sep="\t", header=T, quote="", dec = ".", fill=NA)
if( !is.null(opt$CLINVAR_ALL) & file.exists(opt$CLINVAR_ALL) ) reportList$CLINVAR_ALL <- read.table(opt$CLINVAR_ALL, sep="\t", header=T, quote="", dec = ".", fill=NA)
if( !is.null(opt$PANEL_FILTERED) ) if ( file.exists(opt$PANEL_FILTERED) ) reportList$PANEL_FILTERED <- read.table(opt$PANEL_FILTERED, sep="\t", header=T, quote="", dec = ".", fill=NA)
if( !is.null(opt$PANEL_ALL) ) if ( file.exists(opt$PANEL_ALL) ) reportList$PANEL_ALL <- read.table(opt$PANEL_ALL, sep="\t", header=T, quote="", dec = ".", fill=NA)
if( !is.null(opt$MITOMAP) ) {
  if ( file.exists(opt$MITOMAP) ) {
    tryCatch({reportList$MITOMAP <- read.table(opt$MITOMAP, sep="\t", header=T, quote="", dec = ".", fill=NA)},
             error=function(e){print("Problem reading MITOMAP file - most likely the file is empty!")})
  }
}

# Convert logical columns to single letters - this is to disable T letters being read as logicals
convertLogicalsBackToCharacters <- function(df){
  df[,sapply(df,class) == "logical"] <- sapply(df[,sapply(df,class) == "logical"],
                                               function(i) substr(as.character(i),1,1))
  return(df)
}

for(sheetName in names(reportList)){
  if(nrow(reportList[[sheetName]])==0) next
  reportList[[sheetName]] <- convertLogicalsBackToCharacters(reportList[[sheetName]])
}

# Find compound heterozygous candidates
if("COMPHET_RECESSIVE" %in% names(reportList)){
  if( nrow(reportList$COMPHET_RECESSIVE) > 1) {
  COMPHET_RECESSIVE_gene_entries_list <- lapply(strsplit(reportList$COMPHET_RECESSIVE$ANN....GENE, split = ","), FUN=unique)
  COMPHET_RECESSIVE_gene_list <- unlist(COMPHET_RECESSIVE_gene_entries_list)
  COMPHET_RECESSIVE_comphet_genes <- COMPHET_RECESSIVE_gene_list[duplicated( COMPHET_RECESSIVE_gene_list )]
  COMPHET_RECESSIVE_overlaps<-lapply(COMPHET_RECESSIVE_gene_entries_list, intersect, COMPHET_RECESSIVE_comphet_genes)
  COMPHET_RECESSIVE_overlaps<-lapply(COMPHET_RECESSIVE_overlaps, length)>0
  reportList$COMPHET_RECESSIVE <- reportList$COMPHET_RECESSIVE[COMPHET_RECESSIVE_overlaps,]
  }
}

# Create a concatenated genome coordinate entry
for(sheetName in names(reportList)){
  if( !"ANN....GENE" %in% colnames(reportList[[sheetName]]) ) next # Skip processing of sheets that are not generated by snpEff
  VariantID <- apply(reportList[[sheetName]][,c("CHROM", "POS", "REF", "ALT")], 1, paste0, collapse="-")
  reportList[[sheetName]]$VariantID <- gsub(" ", "", VariantID)
}

# Select only the highest, first snpEff effect
selectFirstEffect<-function(String){
  String <- gsub(",$", ", ", String)
  String <- gsub("^,", " ,", String)
  String <- gsub(",,", ", ,", String)
  AllEffects<-strsplit(String, split=",")[[1]]
  FirstEffect<-AllEffects[1]
  return(FirstEffect)
}
# Test
# String<-c("c.123A>G,c.3345A>G,,")
# selectFirstEffect(String)

for(sheetName in names(reportList)){
  if( !"ANN....GENE" %in% colnames(reportList[[sheetName]]) ) next # Skip processing of sheets that are not generated by snpEff
  if(nrow(reportList[[sheetName]])==0) next
  
  reportList[[sheetName]]$ANN....GENE_SELECTED<-sapply(reportList[[sheetName]]$ANN....GENE, FUN=selectFirstEffect)
  reportList[[sheetName]]$ANN....FEATUREID_SELECTED<-sapply(reportList[[sheetName]]$ANN....FEATUREID, FUN=selectFirstEffect)
  reportList[[sheetName]]$ANN....HGVS_C_SELECTED<-sapply(reportList[[sheetName]]$ANN....HGVS_C, FUN=selectFirstEffect)
  reportList[[sheetName]]$ANN....HGVS_P_SELECTED<-sapply(reportList[[sheetName]]$ANN....HGVS_P, FUN=selectFirstEffect)
  reportList[[sheetName]]$ANN....RANK_SELECTED<-sapply(reportList[[sheetName]]$ANN....RANK, FUN=selectFirstEffect)
  reportList[[sheetName]]$ANN....IMPACT_SELECTED<-sapply(reportList[[sheetName]]$ANN....IMPACT, FUN=selectFirstEffect)
  reportList[[sheetName]]$ANN....EFFECT_SELECTED<-sapply(reportList[[sheetName]]$ANN....EFFECT, FUN=selectFirstEffect)
}

# Variant HGVS nomenclature
makeHGVS<-function(Transcript, Gene, HGVSc){
  #HGVS<-paste0(Transcript, "(", Gene, ")",":", HGVSc)
  HGVS<-paste0(Transcript, ":", HGVSc)
  return(HGVS)
}

# Create a concatenated genome coordinate entry
for(sheetName in names(reportList)){
  if( !"ANN....GENE" %in% colnames(reportList[[sheetName]]) ) next # Skip processing of sheets that are not generated by snpEff
  if(nrow(reportList[[sheetName]])==0) next
  
  reportList[[sheetName]]$HGVS <- mapply(makeHGVS,
                      Transcript=reportList[[sheetName]]$ANN....FEATUREID_SELECTED,
                      Gene=reportList[[sheetName]]$ANN....GENE_SELECTED,
                      HGVSc=reportList[[sheetName]]$ANN....HGVS_C_SELECTED)
}

# Make exon annotation readable
for(sheetName in names(reportList)){
  if( !"ANN....GENE" %in% colnames(reportList[[sheetName]]) ) next # Skip processing of sheets that are not generated by snpEff
  if(nrow(reportList[[sheetName]])==0) next
  
  reportList[[sheetName]]$ANN....RANK_SELECTED <- paste0("Exon ", reportList[[sheetName]]$ANN....RANK_SELECTED)
}


# SpliceAI field processor
reportList[[sheetName]]$SpliceAI.SpliceAI
SpliceAIprocessor <- function(SpliceAI){
  if ( !exists("SpliceAI") ) return("")
  if ( is.na(SpliceAI)) return("")
  if ( SpliceAI == "" ) return("")
  if ( !is.character(SpliceAI) ) return("")
  
  SpliceAIs <- strsplit(SpliceAI, split=",")[[1]]
  SPLICEAI_PREDICTION<-c()
  
  for (SpliceAI in SpliceAIs) {
    if (SpliceAI == "") next
    if (!is.character(SpliceAI)) next
    
    SpliceAI_preds <- as.vector(strsplit(SpliceAI, split="\\|")[[1]])
    names(SpliceAI_preds) <- c("Nucleotide_substitution", 
                               "Gene", 
                               "Acceptor_gain", 
                               "Acceptor_loss", 
                               "Donor_gain", 
                               "Donor_loss",
                               "Acceptor_gain_distance", 
                               "Acceptor_loss_distance", 
                               "Donor_gain_distance", 
                               "Donor_loss_distance"
                               )
    
    
    if( SpliceAI_preds["Acceptor_gain"]>=0.2 & SpliceAI_preds["Acceptor_gain"]<0.5 ) SPLICEAI_PREDICTION <- c(SPLICEAI_PREDICTION, "SPLICEAI_ACCEPTOR_GAIN_SUSPECTED")
    if( SpliceAI_preds["Acceptor_loss"]>=0.2 & SpliceAI_preds["Acceptor_loss"]<0.5 ) SPLICEAI_PREDICTION <- c(SPLICEAI_PREDICTION, "SPLICEAI_ACCEPTOR_LOSS_SUSPECTED")
    if( SpliceAI_preds["Donor_gain"]>=0.2 & SpliceAI_preds["Donor_gain"]<0.5 ) SPLICEAI_PREDICTION <- c(SPLICEAI_PREDICTION, "SPLICEAI_DONOR_GAIN_SUSPECTED")
    if( SpliceAI_preds["Donor_loss"]>=0.2 & SpliceAI_preds["Donor_loss"]<0.5) SPLICEAI_PREDICTION <- c(SPLICEAI_PREDICTION, "SPLICEAI_DONOR_LOSS_SUSPECTED")

    if( SpliceAI_preds["Acceptor_gain"]>=0.5 & SpliceAI_preds["Acceptor_gain"]<0.8 ) SPLICEAI_PREDICTION <- c(SPLICEAI_PREDICTION, "SPLICEAI_ACCEPTOR_GAIN_POSSIBLE")
    if( SpliceAI_preds["Acceptor_loss"]>=0.5 & SpliceAI_preds["Acceptor_loss"]<0.8 ) SPLICEAI_PREDICTION <- c(SPLICEAI_PREDICTION, "SPLICEAI_ACCEPTOR_LOSS_POSSIBLE")
    if( SpliceAI_preds["Donor_gain"]>=0.5 & SpliceAI_preds["Donor_gain"]<0.8 ) SPLICEAI_PREDICTION <- c(SPLICEAI_PREDICTION, "SPLICEAI_DONOR_GAIN_POSSIBLE")
    if( SpliceAI_preds["Donor_loss"]>=0.5 & SpliceAI_preds["Donor_loss"]<0.8) SPLICEAI_PREDICTION <- c(SPLICEAI_PREDICTION, "SPLICEAI_DONOR_LOSS_POSSIBLE")
    
    if( SpliceAI_preds["Acceptor_gain"]>=0.8 ) SPLICEAI_PREDICTION <- c(SPLICEAI_PREDICTION, "SPLICEAI_ACCEPTOR_GAIN_PROBABLE")
    if( SpliceAI_preds["Acceptor_loss"]>=0.8 ) SPLICEAI_PREDICTION <- c(SPLICEAI_PREDICTION, "SPLICEAI_ACCEPTOR_LOSS_PROBABLE")
    if( SpliceAI_preds["Donor_gain"]>=0.8 ) SPLICEAI_PREDICTION <- c(SPLICEAI_PREDICTION, "SPLICEAI_DONOR_GAIN_PROBABLE")
    if( SpliceAI_preds["Donor_loss"]>=0.8 ) SPLICEAI_PREDICTION <- c(SPLICEAI_PREDICTION, "SPLICEAI_DONOR_LOSS_PROBABLE")
    
    }
  
  SPLICEAI_PREDICTION <- paste0(SPLICEAI_PREDICTION, collapse=", ")
  return (SPLICEAI_PREDICTION)
  }

# Test
# SpliceAI="A|IGLL1|0.00|0.00|0.00|0.00|31|-16|31|-26"
# SpliceAI="A|IGLL1|0.80|0.00|0.00|0.00|31|-16|31|-26"
# SpliceAI="A|IGLL1|0.80|0.00|0.00|0.00|31|-16|31|-26,A|IGLL1|0.80|0.00|0.00|0.00|31|-16|31|-26"
# SpliceAI="A|IGLL1|0.80|0.00|0.00|0.00|31|-16|31|-26,"
# SpliceAIprocessor(SpliceAI)


for(sheetName in names(reportList)){
  if( !"ANN....GENE" %in% colnames(reportList[[sheetName]]) ) next # Skip processing of sheets that are not generated by snpEff
  if(nrow(reportList[[sheetName]])==0) next
  
  # Create an empty data frame that will contain prediction flags
  PREDICTIONS <- data.frame( row.names = 1:nrow(reportList[[sheetName]]) )
  
  # Get filtration tags for SpliceAI
  if("SpliceAI.SpliceAI" %in% colnames(reportList[[sheetName]]) ){
    PREDICTIONS$SPLICEAI_PREDICTION <- ""
    PREDICTIONS$SPLICEAI_PREDICTION<-sapply(reportList[[sheetName]]$SpliceAI.SpliceAI, FUN=SpliceAIprocessor)
  }
  
  # Get filtration tags for dbscSNV - ADA algorithm
  if("dbscSNV.ada_score" %in% colnames(reportList[[sheetName]]) ){
    PREDICTIONS$SPLICE_ADA_PATHOGENIC <- ""
    PREDICTIONS$SPLICE_ADA_PATHOGENIC [ reportList[[sheetName]]$dbscSNV.ada_score >= 0.6 ] <- "SPLICE_ADA_PATHOGENIC"
  }
  
  # Get filtration tags for dbscSNV - RF algorithm
  if("dbscSNV.rf_score" %in% colnames(reportList[[sheetName]]) ){
    PREDICTIONS$SPLICE_RF_PATHOGENIC <- ""
    PREDICTIONS$SPLICE_RF_PATHOGENIC [ reportList[[sheetName]]$dbscSNV.rf_score >= 0.6 ] <- "SPLICE_RF_PATHOGENIC"
  }
  
  if("dbNSFP_REVEL_rankscore" %in% colnames(reportList[[sheetName]]) ){
    # REVEL prediction
    # We choose the cutoff as detailed here http://m.ensembl.org/info/genome/variation/prediction/protein_function.html
    # REVEL does not provide a descriptive prediction but for convenience, we display scores above 0.5, as 'likely disease causing' and display scores below 0.5 as 'likely benign'. It was estimated that 75.4% of disease mutations but only 10.9% of neutral variants have a score above 0.5
    PREDICTIONS$REVEL_PREDICTION <- ""
    PREDICTIONS$REVEL_PREDICTION [ reportList[[sheetName]]$dbNSFP_REVEL_rankscore >= 0.5 ] <- "REVEL_PATHOGENIC"
  }
  
  if("dbNSFP_MetaSVM_pred" %in% colnames(reportList[[sheetName]]) ){
    # MetSVM prediction
    PREDICTIONS$METASVM_PREDICTION <- ""
    PREDICTIONS$METASVM_PREDICTION [ reportList[[sheetName]]$dbNSFP_MetaSVM_pred == "D" ] <- "METASVM_PATHOGENIC"
  }
  
  if("dbNSFP_CADD_phred" %in% colnames(reportList[[sheetName]]) ){
    # CADD prediction
    # We choose the cutoff as detailed here http://m.ensembl.org/info/genome/variation/prediction/protein_function.html
    # For convenience, in our transcript tables we display scores above 30 as 'likely deleterious' and scores below as 'likely benign'. Variants with scores over 30 are predicted to be the 0.1% most deleterious possible substitutions in the human genome.
    # For our use I recommend following the evidence on CADD official site which recommend cutoff equal to 15
    PREDICTIONS$CADD_PREDICTION <- ""
    PREDICTIONS$CADD_PREDICTION [ reportList[[sheetName]]$dbNSFP_CADD_phred >= 15 ] <- "CADD_PATHOGENIC"
  }
  
  if("dbNSFP_SIFT_pred" %in% colnames(reportList[[sheetName]]) ){
    # SIFT prediction
    # We choose the cutoff as detailed here http://m.ensembl.org/info/genome/variation/prediction/protein_function.html
    # For convenience, in our transcript tables we display scores above 30 as 'likely deleterious' and scores below as 'likely benign'. Variants with scores over 30 are predicted to be the 0.1% most deleterious possible substitutions in the human genome.
    PREDICTIONS$SIFT_PREDICTION <- ""
    PREDICTIONS$SIFT_PREDICTION [ grep("D", reportList[[sheetName]]$dbNSFP_SIFT_pred) ] <- "SIFT_PATHOGENIC"
  }
  
  if("dbNSFP_Polyphen2_HDIV_pred" %in% colnames(reportList[[sheetName]]) ){
    PREDICTIONS$POLYPHEN2_PREDICTION <- ""
    PREDICTIONS$POLYPHEN2_PREDICTION [ grep("D|P", reportList[[sheetName]]$dbNSFP_Polyphen2_HDIV_pred) ] <- "POLYPHEN_PATHOGENIC"
  }
  
  if("dbNSFP_MutationTaster_pred" %in% colnames(reportList[[sheetName]]) ){
    PREDICTIONS$MUTATIONTASTER_PREDICTION <- ""
    PREDICTIONS$MUTATIONTASTER_PREDICTION [ grep("D|P", reportList[[sheetName]]$dbNSFP_MutationTaster_pred) ] <- "MUTATIONTASTER_PATHOGENIC"
  }
  
  if("dbNSFP_GERP___RS" %in% colnames(reportList[[sheetName]]) ){
    PREDICTIONS$GERP_PREDICTION <- ""
    PREDICTIONS$GERP_PREDICTION [ reportList[[sheetName]]$dbNSFP_GERP___RS > 4 ] <- "GERP_CONSERVED"
  }
  
  PREDICTIONS_TAGS <- apply(PREDICTIONS, 1, paste0, collapse=",")
  #Beautify the final output, remove multiple commas, trailing commas, etc.
  PREDICTIONS_TAGS <- gsub("^,*|(?<=,),|,*$", "", PREDICTIONS_TAGS, perl = T)
  
  # Return a comma separated list of tags for predictions
  reportList[[sheetName]]$Predictions<-PREDICTIONS_TAGS
}

# Genotype file processor 
genotypeProcessor <- function(genotype){
  genotype <- gsub("0/0", "REF", genotype)
  genotype <- gsub("1/0", "HET", genotype)
  genotype <- gsub("0/1", "HET", genotype)
  genotype <- gsub("1/1", "HOM", genotype)
  genotype <- gsub("NA", "NA", genotype)
  #genotype <- gsub("\\.\\/\\.", "NA", genotype)
  genotype[genotype==0]<-"NA"
  return(genotype)
}

qualityProcessor <- function(AD, DP, GT, QUAL, GQ) {
  QUAL_TAGS<-c()

  # If the the GQ field is missing, this is most likely due to the Mutect2 output, therefore quality tags assignment should be skipped
  if ( is.na(GQ) ) { return(paste0(QUAL_TAGS, collapse=",")) }
  if ( is.na(QUAL) ) { return(paste0(QUAL_TAGS, collapse=",")) }
  
  # Assign quality tags
  if ( QUAL<100 ) QUAL_TAGS <- c(QUAL_TAGS, "BAD_VARIANT_QUALITY")
  if ( GQ<90 ) QUAL_TAGS <- c(QUAL_TAGS, "BAD_GENOTYPE_QUALITY")
  if ( DP <10 ) QUAL_TAGS <- c(QUAL_TAGS, "BAD_COVERAGE")
  
  AD <- strsplit(AD, split=",")[[1]]
  AD <- as.numeric(AD)
  if(length(AD)>2) { QUAL_TAGS<-c(QUAL_TAGS, "MULTIALLELIC"); return(paste0(QUAL_TAGS, collapse=",")) }
  
  
  if(length(AD)==2) {
    if ( (AD[1]+AD[2])>0 ) {
      AD_RATIO <- AD[2]/(AD[1]+AD[2])
    } else {
      AD_RATIO <- 0
    }
    if( (GT == "HET" | GT == "0/1" | GT == "1/0") & (AD_RATIO<0.25 || AD_RATIO>0.75) ) QUAL_TAGS<-c(QUAL_TAGS, "BAD_HET_RATIO")
  } 
  return(paste0(QUAL_TAGS, collapse=","))
}
# Test
# qualityProcessor("114,30", "144", "HET", "90", "89")

for(sheetName in names(reportList)){
  if( !"ANN....GENE" %in% colnames(reportList[[sheetName]]) ) next # Skip processing of sheets that are not generated by snpEff
  if(nrow(reportList[[sheetName]])==0) next
  
  # Select colnames that contain genotype and quality data
  genotypeColumns <- grep("GEN.*..GT", colnames(reportList[[sheetName]]))
  ADColumns <- grep("GEN.*..AD", colnames(reportList[[sheetName]]))
  GQColumns <- grep("GEN.*..GQ", colnames(reportList[[sheetName]]))
  DPColumns <- grep("GEN.*..DP", colnames(reportList[[sheetName]]))
  QUALColumns <- grep("QUAL", colnames(reportList[[sheetName]]))
  
  # Rename genotypes for easier readsbility
  for(genotypeColumn in genotypeColumns){
    reportList[[sheetName]][,genotypeColumn] <- sapply(reportList[[sheetName]][,genotypeColumn], FUN=genotypeProcessor)
  }
  
  # Create an empty data frame that will contain prediction flags
  reportList[[sheetName]]$QUALITY <- mapply (qualityProcessor, 
                     AD=reportList[[sheetName]][,ADColumns[1]],
                     DP=reportList[[sheetName]][,DPColumns[1]],
                     GT=reportList[[sheetName]][,genotypeColumns[1]],
                     QUAL=reportList[[sheetName]][,QUALColumns[1]],
                     GQ=reportList[[sheetName]][,GQColumns[1]])
}

# Add functional candidate tags to simplify interpretation
functionTags <- function(pLI, misZ, impact, effect, GT, gnomADexomes.AC, gnomAD.AC, SLOpopulation.AC_Het, VariantPredictions="", Quality="") {
  FUNC_TAGS<-c()
  if( is.na(pLI) ) pLI = 0
  if( is.na(misZ) ) misZ = 0
  if( is.na(impact) ) impact = "NONE_RETURNED"
  if( is.na(effect) ) effect = "NONE_RETURNED"
  
  isQualityOK <- function(Quality=Quality){
    return( !grepl("BAD_VARIANT_QUALITY|BAD_HET_RATIO", Quality, ignore.case = T) )
  }
  
  Metapredictor_score <- 0
  for( Predictor in c("REVEL_PATHOGENIC", "METASVM_PATHOGENIC", "CADD_PATHOGENIC") ) {
      if( grepl(Predictor, VariantPredictions) ) Metapredictor_score <- Metapredictor_score + 1; print(Metapredictor_score)
  }
  
  Predictor_score <- 0
  for( Predictor in c("SIFT_PATHOGENIC", "POLYPHEN_PATHOGENIC", "MUTATIONTASTER_PATHOGENIC") ) {
    if( grepl(Predictor, VariantPredictions) ) Predictor_score <- Predictor_score + 1
  }
  
  if( Metapredictor_score==3 ) {
    FUNC_TAGS <-c(FUNC_TAGS, "CONSENSUS_METAPREDICTOR_PATHOGENIC")
  } else if (Metapredictor_score + Predictor_score >= 4) {
    FUNC_TAGS <-c(FUNC_TAGS, "MAJORITY_PREDICTOR_PATHOGENIC")
  }
  
  if(  isQualityOK(Quality) & (gnomADexomes.AC < 3 || is.na(gnomADexomes.AC)) & (gnomAD.AC < 3 || is.na(gnomAD.AC)) & (SLOpopulation.AC_Het < 3 || is.na(SLOpopulation.AC_Het)) & grepl("Splice", VariantPredictions, ignore.case = T) ) FUNC_TAGS <-c(FUNC_TAGS, "SPLICE_CANDIDATE")
  
  if( grepl("GERP", VariantPredictions, ignore.case = T) ) FUNC_TAGS <-c(FUNC_TAGS, "CONSERVED")
  
  if( isQualityOK(Quality) & (gnomADexomes.AC < 1 || is.na(gnomADexomes.AC)) & (gnomAD.AC < 1 || is.na(gnomAD.AC)) & (SLOpopulation.AC_Het < 3 || is.na(SLOpopulation.AC_Het)) ) FUNC_TAGS <-c(FUNC_TAGS, "ULTRARARE_CANDIDATE")
  if( isQualityOK(Quality) & GT=="HOM" & (gnomADexomes.AC < 20 || is.na(gnomADexomes.AC)) & (gnomAD.AC < 20 || is.na(gnomAD.AC)) & (SLOpopulation.AC_Het < 5 || is.na(SLOpopulation.AC_Het)) ) FUNC_TAGS <-c(FUNC_TAGS, "ULTRARARE_HOM_CANDIDATE")
  if( isQualityOK(Quality) & as.numeric(pLI) > 0.9 & impact=="HIGH" & (GT == "HET" | GT == "0/1" | GT == "1/0") & (gnomADexomes.AC < 3 || is.na(gnomADexomes.AC)) ) FUNC_TAGS <- c(FUNC_TAGS, "LOF_CANDIDATE")
  if( isQualityOK(Quality) & as.numeric(misZ) > 2 & effect=="missense_variant" & (GT == "HET" | GT == "0/1" | GT == "1/0") & (gnomADexomes.AC < 3 || is.na(gnomADexomes.AC)) ) FUNC_TAGS <- c(FUNC_TAGS, "MISSENSE_CANDIDATE")
  if( isQualityOK(Quality) & (gnomADexomes.AC < 3 || is.na(gnomADexomes.AC)) & (gnomAD.AC < 3 || is.na(gnomAD.AC)) & (SLOpopulation.AC_Het < 3 || is.na(SLOpopulation.AC_Het)) & impact=="HIGH" ) FUNC_TAGS <- c(FUNC_TAGS, "RARE_LOF_CANDIDATE")
  return(paste0(FUNC_TAGS, collapse=","))
}
# Tests
# functionTags(pLI=0.99, misZ = 0, GT = "HET", impact = "HIGH", effect = "stopgain_variant", gnomADexomes.AC = NA, gnomAD.AC = NA, SLOpopulation.AC_Het = NA)
# functionTags(pLI=0.99, misZ = 2.3, GT = "HET", impact = "MODERATE", effect = "missense_variant", gnomADexomes.AC = NA, gnomAD.AC = NA, SLOpopulation.AC_Het = NA)
# functionTags(pLI=0.99, misZ = 0, GT = "HET", impact = "HIGH", effect = "stopgain_variant", gnomADexomes.AC = NA, gnomAD.AC = NA, SLOpopulation.AC_Het = NA, VariantPredictions = "REVEL_PATHOGENIC,CADD_PATHOGENIC,SIFT_PATHOGENIC,POLYPHEN_PATHOGENIC,MUTATIONTASTER_PATHOGENIC")

for(sheetName in names(reportList)){
  if( !"ANN....GENE" %in% colnames(reportList[[sheetName]]) ) next # Skip processing of sheets that are not generated by snpEff
  if(nrow(reportList[[sheetName]])==0) next
  
  # Select colnames that contain genotype and quality data
  genotypeColumns <- grep("GEN.*..GT", colnames(reportList[[sheetName]]))
  
  # Create an empty data frame that will contain prediction flags
  reportList[[sheetName]]$FUNCTIONAL_TAGS <- mapply (functionTags, 
                                             pLI = reportList[[sheetName]]$pLI,
                                             misZ = reportList[[sheetName]]$oe_mis,
                                             effect = reportList[[sheetName]]$ANN....EFFECT_SELECTED,
                                             impact = reportList[[sheetName]]$ANN....IMPACT_SELECTED,
                                             gnomADexomes.AC = reportList[[sheetName]]$gnomADexomes.AC,
                                             gnomAD.AC = reportList[[sheetName]]$gnomAD.AC,
                                             SLOpopulation.AC_Het = reportList[[sheetName]]$SLOpopulation.AC_Het,
                                             VariantPredictions = reportList[[sheetName]]$Predictions,
                                             Quality = reportList[[sheetName]]$QUALITY,
                                             GT = reportList[[sheetName]][,genotypeColumns[1]])
}


for(sheetName in names(reportList)){
  if( !"ANN....GENE" %in% colnames(reportList[[sheetName]]) ) next # Skip processing of sheets that are not generated by snpEff
  if(nrow(reportList[[sheetName]])==0) next
  
  allColumns <- colnames(reportList[[sheetName]])
  genotypeColumns <- grep("GEN\\.*", colnames(reportList[[sheetName]]))
  genotypeColumnNames <- colnames(reportList[[sheetName]])[genotypeColumns]
  
  columnsOrder<-c("VariantID", 
                genotypeColumnNames,
                "QUALITY",
                "HGVS", 
                "ANN....GENE_SELECTED", 
                "ANN....GENE",
                "ANN....RANK_SELECTED",
                "ANN....RANK",
                "ANN....HGVS_P_SELECTED",
                "ANN....HGVS_P",
                "ANN....IMPACT_SELECTED",
                "ANN....IMPACT",
                "ANN....EFFECT_SELECTED",
                "ANN....EFFECT",
                "Predictions",
                "FUNCTIONAL_TAGS",
                "Disease_name", "OMIM", "HPO", "Categorization", "Inheritance", "Age",
                "clinvar.ID", "clinvar.CLNSIG", "clinvar.CLNDN", "clinvar.CLNSIGCONF", "clinvar.CLNHGVS",
                "SLOpopulation.AC_Het", "SLOpopulation.AC_Hom", "SLOpopulation.AC_Hemi",
                "gnomAD.AC", "gnomAD.AF", "gnomAD.nhomalt",
                "gnomADexomes.AC", "gnomADexomes.AF", "gnomADexomes.nhomalt",
                "dbNSFP_REVEL_rankscore", "dbNSFP_MetaSVM_pred", "dbNSFP_CADD_phred", "dbNSFP_DANN_rankscore", "dbNSFP_SIFT_pred", "dbNSFP_SIFT4G_pred", "dbNSFP_Polyphen2_HDIV_pred", "dbNSFP_MutationTaster_pred", "dbNSFP_PrimateAI_pred", "dbNSFP_Polyphen2_HDIV_score",
                "SpliceAI.SpliceAI", "dbscSNV.ada_score", "dbscSNV.rf_score",
                "dbNSFP_GERP___NR", "dbNSFP_GERP___RS",
                "dbNSFP_Interpro_domain",
                "pLI", "oe_mis", "pRec",
                "LOF....GENE", "LOF....GENEID", "LOF....NUMTR", "LOF....PERC", 
                "NMD....GENE", "NMD....GENEID", "NMD....NUMTR", "NMD....PERC"
                )
  
  columnsOrder <- unique(columnsOrder[columnsOrder %in% allColumns])
  reportList[[sheetName]] <- reportList[[sheetName]][,columnsOrder]
}

# Add the interpretation fields for adding interpretations to the Excel file - sheet Panel filtered
INTERPRETATION_FIELDS<-c("Classification", "ACMG", "Report", "InterpretationSI", "InterpretationEN", "Condition", "Origin", "Interpretation_date")
if( !is.null(reportList[["PANEL_FILTERED"]]) ) if( nrow(reportList[["PANEL_FILTERED"]])>0 ) for (i in INTERPRETATION_FIELDS) reportList$PANEL_FILTERED[[i]]<-""

# Create a novel work book 
wb <- openxlsx::createWorkbook()

# Create sheets for each filtration set
for(sheetName in names(reportList)){
  openxlsx::addWorksheet(wb, sheetName, zoom=90)
  openxlsx::writeData(wb, sheetName, reportList[[sheetName]])
  openxlsx::addFilter(wb, sheetName, row = 1, cols = 1:ncol(reportList[[sheetName]]))
  
  # Set the widths of columns and hide selected columns
  ColNames<-colnames(reportList[[sheetName]])
  
  # Define column widths in the final report
  ColWidths<-rep(8, ncol(reportList[[sheetName]]))
  ColWidths[grepl("VariantID",ColNames)]<-20
  ColWidths[grepl("CHR",ColNames)]<-6
  ColWidths[grepl("POS",ColNames)]<-12
  ColWidths[grepl("REF",ColNames)]<-5
  ColWidths[grepl("ALT",ColNames)]<-5
  ColWidths[grepl("QUAL",ColNames)]<-5
  ColWidths[grepl("QUALITY",ColNames)]<-16
  ColWidths[grepl("..GT",ColNames)]<-5
  ColWidths[grepl("..AD",ColNames)]<-7
  ColWidths[grepl("..DP",ColNames)]<-4
  ColWidths[grepl("..GQ",ColNames)]<-4
  ColWidths[grepl("ANN....GENE",ColNames)]<-10
  ColWidths[grepl("Disease_name",ColNames)]<-24
  ColWidths[grepl("Categorization",ColNames)]<-14
  ColWidths[grepl("Inheritance",ColNames)]<-6
  ColWidths[grepl("Age",ColNames)]<-10
  ColWidths[grepl("HPO",ColNames)]<-20
  ColWidths[grepl("OMIM",ColNames)]<-20
  ColWidths[grepl("FEATUREID",ColNames)]<-12
  ColWidths[grepl("HGVS_C",ColNames)]<-12
  ColWidths[grepl("HGVS_P",ColNames)]<-12
  ColWidths[ColNames=="HGVS"]<-30
  ColWidths[grepl("EFFECT",ColNames)]<-18
  ColWidths[grepl("IMPACT",ColNames)]<-12
  ColWidths[grepl("FUNCTIONAL_TAGS",ColNames)]<-12
  ColWidths[grepl("Predictions",ColNames)]<-24
  ColWidths[grepl("SLOpopulation",ColNames)]<-5
  ColWidths[grepl("gnomAD",ColNames)]<-5
  ColWidths[grepl("CLNSIG",ColNames)]<-30
  ColWidths[grepl("CLNDN",ColNames)]<-30
  ColWidths[grepl("SpliceAI",ColNames)]<-18
  ColWidths[grepl("dbNSFP_Interpro_domain",ColNames)]<-40

  
  # Define which columns to hide
  colsToHide <- which(ColNames %in% c( 
                        "ANN....GENE",
                        "ANN....RANK",
                        "ANN....HGVS_P",
                        "ANN....IMPACT",
                        "ANN....EFFECT",
                        "OMIM", "HPO",
                        "clinvar.CLNHGVS",
                        "dbNSFP_REVEL_rankscore", "dbNSFP_MetaSVM_pred", "dbNSFP_CADD_phred", "dbNSFP_DANN_rankscore", "dbNSFP_SIFT_pred", "dbNSFP_SIFT4G_pred", "dbNSFP_Polyphen2_HDIV_pred", "dbNSFP_MutationTaster_pred", "dbNSFP_PrimateAI_pred", "dbNSFP_Polyphen2_HDIV_score",
                        "SpliceAI.SpliceAI", "dbscSNV.ada_score", "dbscSNV.rf_score",
                        "dbNSFP_GERP___NR", "dbNSFP_GERP___RS",
                        "LOF....GENE", "LOF....GENEID", "LOF....NUMTR", "LOF....PERC", 
                        "NMD....GENE", "NMD....GENEID", "NMD....NUMTR", "NMD....PERC"
  ))
  colsToHide <- unique(c(grep("*..DP", ColNames), colsToHide))
  colsToHide <- unique(c(grep("*..GQ", ColNames), colsToHide))
  ColHide<-rep(FALSE, ncol(reportList[[sheetName]]))
  ColHide[colsToHide]<-TRUE

  # Apply columns widths and hide the defined columns
  openxlsx::setColWidths(wb, sheetName, cols=1:length(ColWidths), widths = ColWidths, hidden = ColHide, ignoreMergedCells = FALSE)
  
  ## Styling the sheet
  headerStyle <- createStyle(fontColour = "#636363", fgFill = "#9ecae1",
                             halign = "center", valign = "center", textDecoration = "Bold",
                             border = "TopBottom", textRotation = 25)
  
  addStyle(wb, sheet = sheetName, headerStyle, rows = 1, cols = 1:length(ColWidths), gridExpand = TRUE)
  #bodyStyle <- createStyle(border = "TopBottom", borderColour = "#4F81BD") # This results in a corrupted XLSX file if the output only contains one row, therefore currently disabling
  #addStyle(wb, sheet = sheetName, bodyStyle, rows = 2:nrow(reportList[[sheetName]]), cols = 1:ncol(reportList[[sheetName]]), gridExpand = TRUE)
}

openxlsx::saveWorkbook(wb, file = opt$XLSX_OUTPUT, overwrite = TRUE)



