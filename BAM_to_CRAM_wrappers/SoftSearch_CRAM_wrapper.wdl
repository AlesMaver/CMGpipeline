version 1.0
## Copyright CMG@KIGM, Ales Maver, Tanja Majnik

# Wrapping WF 
# using old Conifer wf with BAM file input, but having a CRAM file for the input
# and calling CRAM to BAM file conversion WF

import "../CRAM_conversions.wdl" as CramConversions 
import "../softsearch/softsearch_v2.wdl" as SoftSearch

# WORKFLOW DEFINITION 
workflow SoftSearch_CRAM_wrapper {
  input {
    String sample_basename
    File? input_cram
    File? input_cram_index
    File reference_fa
    File reference_fai
    File reference_dict
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

  call SoftSearch.SoftSearchWF as SoftSearch {
      input:
        input_bam = CramToBam.output_bam,
        input_bam_index = CramToBam.output_bai,
        reference_fa = reference_fa,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        sample_basename = sample_basename
  }

  output {
    File output_complete_vcf = SoftSearch.output_complete_vcf
    File output_complete_vcf_index = SoftSearch.output_complete_vcf_index
    File output_vcf = SoftSearch.output_vcf
    File? output_tsv_name = SoftSearch.output_tsv_name
  }
}
