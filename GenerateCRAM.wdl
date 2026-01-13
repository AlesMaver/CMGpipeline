version 1.0

import "./CRAM_conversions.wdl" as CramConversions 

# WORKFLOW DEFINITION 
workflow GenerateCRAM {
  input {
    File input_bam
    File reference_fa
    File reference_fai
    String sample_basename
  }
  
call CramConversions.ConvertToCram as ConvertToCram {
    input:
      input_bam = input_bam,
      ref_fasta = reference_fa,
      ref_fasta_index = reference_fai,
      sample_basename = sample_basename
  }
  
  
  output {
    File output_cram = ConvertToCram.output_cram
    File output_cram_index = ConvertToCram.output_cram_index
  }
  
  # end-workflow
}
