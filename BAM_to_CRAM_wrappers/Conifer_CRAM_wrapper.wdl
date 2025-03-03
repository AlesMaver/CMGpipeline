version 1.0
## Copyright CMG@KIGM, Ales Maver, Tanja Majnik

# Wrapping WF 
# using old Conifer wf with BAM file input, but having a CRAM file for the input
# and calling CRAM to BAM file conversion WF

import "../CRAM_conversions.wdl" as CramConversions 
import "../Conifer.wdl" as Conifer

# WORKFLOW DEFINITION 
workflow Conifer_CRAM_wrapper {
  input {
    String sample_basename
    File? input_cram
    File? input_cram_index
    File reference_fa
    File reference_fai
    File reference_dict
    #Array[File]? input_reference_rpkms 
    Array[File] input_reference_rpkms = select_first([input_reference_rpkms, [""]])
    Int? CONIFER_svd
    Float? CONIFER_threshold
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

  call Conifer.Conifer as Conifer {
      input:
        input_bam = CramToBam.output_bam,
        input_bam_index = CramToBam.output_bai,
        input_reference_rpkms = input_reference_rpkms,
        CONIFER_svd = CONIFER_svd,
        CONIFER_threshold = CONIFER_threshold,
        enrichment = enrichment,
        enrichment_bed = enrichment_bed
  }

  output {
    File output_conifer_calls = Conifer.output_conifer_calls
    File output_conifer_calls_wig = Conifer.output_conifer_calls_wig
    # Array[File] output_plotcalls = Conifer.output_plotcalls
    File conifer_plots_tar = Conifer.conifer_plots_tar
    File CNV_bed = Conifer.CNV_bed
    File CNV_wig = Conifer.CNV_wig
    File output_rpkm = Conifer.output_rpkm
    File? annotSV_tsv = Conifer.annotSV_tsv
  }
}

