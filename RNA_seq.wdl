version 1.0

## orig: import "https://raw.githubusercontent.com/AlesMaver/bioWDL_RNAseq/develop/align-hisat2.wdl" as AlignHisat2
import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/RNA_seq/align-hisat2.wdl" as AlignHisat2
import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/GenerateCRAM.wdl" as GenerateCRAM

workflow RNA_seq {
  input {
        String sample_basename

        # data for AlignHisat2:
        ### Array[FastqPair]+ inputReads 
        File input_fq1
        File input_fq2

        String outputDir = "."
        String library
        Array[String] readgroups
        Array[File]+ indexFiles

        # references for CRAM generation:
        File reference_fa
        File reference_fai
  }

  #FastqPair inputReads =  { 
  #            "R1": input_fq1,
  #            "R2": input_fq2 
  #          }

  Array[FastqPair] inputReads =  [
    { "R1": input_fq1, "R2": input_fq2 }
  ]


  call AlignHisat2.AlignHisat2 as AlignHisat2 {
      input:
        inputReads = inputReads, 
        outputDir = outputDir,
        sample = sample_basename,
        library = library,
        readgroups = readgroups,
        ### platform = "illumina",
        indexFiles = indexFiles
  }

  call GenerateCRAM.GenerateCRAM as GenerateCRAM {
      input:       
        input_bam = AlignHisat2.bamFile["file"],
        reference_fa = reference_fa,
        reference_fai = reference_fai,
        sample_basename = sample_basename
  }

  output {
    #File output_bam = AlignHisat2.bamFile["file"] 
    #File output_bam_index = AlignHisat2.bamFile["index"]
    
    File output_cram = GenerateCRAM.output_cram 
    File output_cram_index = GenerateCRAM.output_cram_index 
  }
}

struct FastqPair {
    File R1
    String? R1_md5
    File? R2
    String? R2_md5
}
