version 1.0
## Copyright CMG@KIGM, Ales Maver, Tanja Majnik

# Wrapping WF 
# using old  wf with BAM file input, but having a CRAM file for the input
# and calling CRAM to BAM file conversion WF

import "../CRAM_conversions.wdl" as CramConversions 
import "../DeepVariant.wdl" as DeepVariant

# WORKFLOW DEFINITION 
workflow DeepVariant_CRAM_wrapper {
  input {
    String sample_basename
    File? input_cram
    File? input_cram_index
    File reference_fa
    File reference_fai
    File reference_dict
    String modelType
    Int? numShards
  } 

  # String sample_basename = sub(basename(input_cram), "[\_,\.].*", "" )

  call CramConversions.CramToBam as CramToBam {
      input:
        sample_name = sample_basename,
        input_cram = input_cram,
        ref_fasta = reference_fa,
        ref_fasta_index = reference_fai,
        ref_dict = reference_dict,
        docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817",
        samtools_path = "samtools"
  }

  call DeepVariant.DeepVariant as DeepVariant {
      input:
            sample_basename = sample_basename,
            inputBam = CramToBam.output_bam,
            inputBamIndex = CramToBam.output_bai,
            referenceFasta = reference_fa,
            referenceFastaIndex = reference_fai,
            modelType = modelType,
            numShards = numShards
  }

  output {
      File outputVCF = DeepVariant.outputVCF
      File outputVCFIndex = DeepVariant.outputVCFIndex
      File? outputVCFStatsReport = DeepVariant.outputVCFStatsReport
      File? outputGVCF = DeepVariant.outputGVCF
      File? outputGVCFIndex = DeepVariant.outputGVCFIndex
  }
}
