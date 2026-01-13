# Load required libraries
library("optparse")

# Parse command line options
option_list = list(
  make_option(c("-i", "--annotatedIntervalsFile"), type="character", default=NULL,
              help="Outout of bedtools annotated intervals file", metavar="character")
)
opt_parser = optparse::OptionParser(option_list=option_list); opt = optparse::parse_args(opt_parser)

# This is function that will take a list of annotations and collapse them into a comma delimited format. 
# It will also perform some cleaning - ie. replace the underscores with spaces
collapseAnnotations=function(x, annotation_type){
  # Return unique and non-empty annotations
  x = unique(na.omit(x))
  x = x[x!=""]
  
  # If a single annotation is present, just return it
  if(length(x)==1) { return(x) }
  
  # Collapse additional annotations in a single field
  annotations = paste(x[2:length(x)], collapse = ", ")
  
  # Remove ugly underscores in annotations
  annotations = gsub("_", " ", annotations)
  
  # Make the pretty annotation format if more than 1 annotations are present, the first field is usually gene name and in brackets, other annotations will be included
  paste(c(x[1], " (", annotations , ")"), collapse="")
}

# Input annotation file 
annotated_intervals_file= opt$annotatedIntervalsFile

# Determine number of fields in the interval file
annotated_intervals_file_numfields = count.fields(annotated_intervals_file, sep = "\t")

if (length(annotated_intervals_file_numfields) > 0) {
  ncol=max(na.omit(annotated_intervals_file_numfields))

  annotated_intervals = read.table(annotated_intervals_file, fill=TRUE, sep = "\t", na.strings = "N/A", col.names = c("chr", "start", "stop", 4:ncol))
  annotation_types=unique(annotated_intervals$X5)
  annotation_type = annotation_types[1]
  createIntervalID=function(x){paste(trimws(x[1:3]),sep="-", collapse="-")}
  annotated_intervals=cbind(uniqueID=apply(annotated_intervals, 1, createIntervalID), annotated_intervals)

  interval_list=data.frame(uniqueID = unique(annotated_intervals$uniqueID), chr="", start="", stop="")
  for(annotation_type in annotation_types) {
    interval_list[,annotation_type]=""
  }
  interval_list[,"omim_hgncs"]=""

  for (interval in unique(annotated_intervals$uniqueID)){
    interval_list[interval_list$uniqueID == interval,c("chr", "start", "stop")]=annotated_intervals[annotated_intervals$uniqueID==interval,c("chr", "start", "stop")][1,]

    for (annotation_type in annotation_types) {
      DF = annotated_intervals[annotated_intervals$uniqueID == interval & annotated_intervals$X5 == annotation_type,]
      DF$collapsedAnnotations = apply(DF[10:ncol], 1, collapseAnnotations)
      interval_list[interval_list$uniqueID == interval,annotation_type] = paste(unique(DF$collapsedAnnotations, annotation_type=annotation_type), collapse="; ")
    }

    OMIM_DF=annotated_intervals[annotated_intervals$uniqueID==interval & annotated_intervals$X5=="OMIM",]
    if(nrow(OMIM_DF)>0){
      interval_list[interval_list$uniqueID == interval,"omim_hgncs"] = paste( unique(OMIM_DF[,10]), collapse = ", ")
    }
  }
} else {
  interval_list = ""
}
  
write.table(interval_list, file="output.txt", sep="\t", quote=F, row.names = F)
