version 1.0
## Copyright CMG@KIGM, Ales Maver, Tanja Majnik

# Wrapping WF 
# using old wf with BAM file input, but having a CRAM file for the input
# and calling CRAM to BAM file conversion WF

import "../CRAM_conversions.wdl" as CramConversions 
import "../manta/manta_workflow.wdl" as Manta

# WORKFLOW DEFINITION 
workflow Manta_CRAM_wrapper {
  input {
    String sample_basename
    File? input_cram
    File? input_cram_index
    File reference_fa
    File reference_fai
    File reference_dict
    Boolean exome = false
    Array[File] input_manta_reference_vcfs
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

  call Manta.SVcalling as Manta {
      input:
        bamFile = CramToBam.output_bam,
        bamIndex = CramToBam.output_bai,
        referenceFasta = reference_fa,
        referenceFastaFai = reference_fai,
        referenceFastaDict = reference_dict,  
        exome = false,
        sample = sample_basename, 
        input_manta_reference_vcfs = input_manta_reference_vcfs
  }

  output {
        File mantaVcf = Manta.mantaVcf
        File mantaVcfindex = Manta.mantaVcfindex
        File output_sv_table = Manta.output_sv_table
        File output_manta_filtered_vcf = Manta.output_manta_filtered_vcf
        File? annotSV_tsv = Manta.annotSV_tsv
  }
}
