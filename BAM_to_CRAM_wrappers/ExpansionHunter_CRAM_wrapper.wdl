version 1.0
## Copyright CMG@KIGM, Ales Maver, Tanja Majnik

# Wrapping WF 
# using old wf with BAM file input, but having a CRAM file for the input
# and calling CRAM to BAM file conversion WF

import "../CRAM_conversions.wdl" as CramConversions 
import "../exp_hunter.wdl" as ExpansionHunter 

# WORKFLOW DEFINITION 
workflow ExpansionHunter_CRAM_wrapper {
  input {
    String sample_basename
    File? input_cram
    File? input_cram_index
    File reference_fa
    File reference_fai
    File reference_dict
    String expansion_hunter_docker
  } 

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

  call ExpansionHunter.ExpansionHunter as ExpansionHunter {
      input:
        sample_id = sample_basename,
        bam_file = CramToBam.output_bam,
        bai_file = CramToBam.output_bai,
        reference_fasta = reference_fa,
        expansion_hunter_docker = expansion_hunter_docker
  }

  output {
    File? expansion_hunter_vcf_annotated = ExpansionHunter.expansion_hunter_vcf_annotated
  }
}
