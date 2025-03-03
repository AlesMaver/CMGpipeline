version 1.0
## Copyright CMG@KIGM, Ales Maver, Tanja Majnik

# Wrapping WF 
# using old wf with BAM file input, but having a CRAM file for the input
# and calling CRAM to BAM file conversion WF

import "../CRAM_conversions.wdl" as CramConversions 
import "../Conifer_OnlyRPKM.wdl" as Conifer_OnlyRPKM

# WORKFLOW DEFINITION 
workflow Conifer_OnlyRPKM_CRAM_wrapper {
  input {
    String sample_basename
    File? input_cram
    File? input_cram_index
    File reference_fa
    File reference_fai
    File reference_dict
    String? enrichment
    File? enrichment_bed
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

  call Conifer_OnlyRPKM.Conifer_OnlyRPKM as Conifer_OnlyRPKM {
      input:
        input_bam = CramToBam.output_bam,
        input_bam_index = CramToBam.output_bai,
        enrichment = enrichment,
        enrichment_bed = enrichment_bed
  }

  output {
    File output_rpkm = Conifer_OnlyRPKM.output_rpkm
  }
}
