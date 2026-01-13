# This is a script that converts the RARE_FUNCTIONAL variants in our XLSX file to an input ready for SimulConsult tool

#!/usr/bin/env Rscript
library("optparse")
library("openxlsx")
options(stringsAsFactors = F)

# Parse input options
option_list = list(
  make_option(c("--XLSX_INPUT"), type="character", default=NULL, help="Annotated XLSX, containing the RARE_FUNCTIONAL sheet", metavar="character"),
  make_option(c("--sample_name"), type="character", default=NULL, help="Sample name", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Testing 
# opt <- list()
# opt$XLSX_INPUT = "EXOME.FinalReportNew.xlsx"

# Exit if no input file has been provided
if(is.null(opt$XLSX_INPUT)) {print("No input XLSX file defined! Exiting..."); quit(save = "no", status = 1, runLast = FALSE)}

# Get EXOME id from the file name
# Currently disabling this option in favor of a variable passed to the R script
# EXOME = gsub(".*\\/([A-Za-z0-9]*)\\..*", "\\1", opt$XLSX_INPUT)

# Headers required by the SimulConsult software
headers<-"hgncSymbol	geneNameLong	chrPos	cSeqAnnotation	cPosition	cRef	cAlt	pSeqAnnotation	pPosition	pRef	pAlt	rsid	zygProband	zygMother	zygFather	effect	freq1	freq2	homoShares	heteroShares	omimNumber	omimDiseaseNames	variantAccession	variantPathogenicity	polyPhen	mutationTaster	sift	gerp	grantham	phat	phast	phyloP	strandBias	knownSplice	totDepthProband	varDepthProband	qualProband	totDepthMother	varDepthMother	qualMother	totDepthFather	varDepthFather	qualFather"
headers<-strsplit(headers, split = "\t")[[1]]

# Read the rare functional variants from the annotated XLSX file
RARE_FUNCTIONAL <- openxlsx::read.xlsx(opt$XLSX_INPUT, sheet = "RARE_FUNCTIONAL")
RARE_FUNCTIONAL <- RARE_FUNCTIONAL[grep("BAD_HET_RATIO|BAD_VARIANT_QUALITY|BAD_COVERAGE", RARE_FUNCTIONAL$QUALITY, invert = T),]

# Create a data.frame, with columns in compatibility with SimulConsult
df <- data.frame(matrix(nrow = nrow(RARE_FUNCTIONAL), ncol=length(headers)))
colnames(df)<- headers

# Fill in the data
df$hgncSymbol<-RARE_FUNCTIONAL$ANN....GENE_SELECTED
df$chrPos<-gsub("^(chr[0-9XY]*)-([0-9]*).*", "\\1:\\2", RARE_FUNCTIONAL$VariantID)
df$cSeqAnnotation <- RARE_FUNCTIONAL$HGVS
GENOTYPE_PROBAND_COLUMN_NUMBER <- grep("GEN.*GT", colnames(RARE_FUNCTIONAL))[1]
df$zygProband <- RARE_FUNCTIONAL[,GENOTYPE_PROBAND_COLUMN_NUMBER]
df$effect <- RARE_FUNCTIONAL$ANN....EFFECT_SELECTED
df$freq1 <- format(RARE_FUNCTIONAL$gnomAD.AF, scientific = F)

# Write out the text file containing the SimulConsult inputs
write.table(df, file=paste0(opt$sample_name, ".SimulConsult.input.txt"), sep="\t", quote = F, row.names = F)
