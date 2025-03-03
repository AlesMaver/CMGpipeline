version 1.0
## Copyright CMG@KIGM, Ales Maver, Tanja Majnik

# Wrapping WF 
# using old wf with BAM file input, but having a CRAM file for the input
# and calling CRAM to BAM file conversion WF

import "../CRAM_conversions.wdl" as CramConversions 
import "../ExomeDepth.wdl" as ExomeDepth

# WORKFLOW DEFINITION 
workflow ExomeDepth_CRAM_wrapper {
  input {
    String sample_basename
    File? input_cram
    File? input_cram_index
    File reference_fa
    File reference_fai
    File reference_dict
    File? target_bed
    File? exome_depth_counts_input
    Array[File]? reference_counts_files

    File? input_bam
  } 


  if(!defined(exome_depth_counts_input)) {
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
    call ExomeDepth.ExomeDepth as ExomeDepth {
      input:
        input_bam = CramToBam.output_bam,
        input_bam_index = CramToBam.output_bai,
        sample_name = sample_basename,
        target_bed = target_bed,
        exome_depth_counts_input = exome_depth_counts_input,
        reference_counts_files = reference_counts_files
    }
  }

  if(defined(exome_depth_counts_input)) {
    call ExomeDepth.ExomeDepth as ExomeDepthWithoutBAM {
      input:
        #input_bam = select_first([CramToBam.output_bam, input_bam]),
        #input_bam_index = select_first([CramToBam.output_bai,""]),
        sample_name = sample_basename,
        target_bed = target_bed,
        exome_depth_counts_input = exome_depth_counts_input,
        reference_counts_files = reference_counts_files
    }    
  }

  output {
    File? exome_depth_counts_1 = ExomeDepthWithoutBAM.exome_depth_counts
    File? exome_depth_cnv_calls_bed_1 = ExomeDepthWithoutBAM.exome_depth_cnv_calls_bed
    File? exome_depth_cnv_calls_csv_1 = ExomeDepthWithoutBAM.exome_depth_cnv_calls_csv
    File? exome_depth_ratios_all_wig_gz_1 = ExomeDepthWithoutBAM.exome_depth_ratios_all_wig_gz
    File? exome_depth_ratios_all_wig_gz_tbi_1 = ExomeDepthWithoutBAM.exome_depth_ratios_all_wig_gz_tbi
    File? exome_depth_rolling_ratios_wig_1 = ExomeDepthWithoutBAM.exome_depth_rolling_ratios_wig
    File? exome_depth_rolling_ratios_wig_gz_tbi_1 = ExomeDepthWithoutBAM.exome_depth_rolling_ratios_wig_gz_tbi
    File? exome_depth_ratios_clean_wig_gz_1 = ExomeDepthWithoutBAM.exome_depth_ratios_clean_wig_gz
    File? exome_depth_ratios_clean_wig_gz_tbi_1 = ExomeDepthWithoutBAM.exome_depth_ratios_clean_wig_gz_tbi
    File? exome_depth_ratios_nomissing_wig_gz_1 = ExomeDepthWithoutBAM.exome_depth_ratios_nomissing_wig_gz
    File? exome_depth_ratios_nomissing_wig_gz_tbi_1 = ExomeDepthWithoutBAM.exome_depth_ratios_nomissing_wig_gz_tbi
    File? exome_depth_annotSV_tsv_1 = ExomeDepthWithoutBAM.exome_depth_annotSV_tsv

    File? exome_depth_counts = ExomeDepth.exome_depth_counts
    File? exome_depth_cnv_calls_bed = ExomeDepth.exome_depth_cnv_calls_bed
    File? exome_depth_cnv_calls_csv = ExomeDepth.exome_depth_cnv_calls_csv
    File? exome_depth_ratios_all_wig_gz = ExomeDepth.exome_depth_ratios_all_wig_gz
    File? exome_depth_ratios_all_wig_gz_tbi = ExomeDepth.exome_depth_ratios_all_wig_gz_tbi
    File? exome_depth_rolling_ratios_wig = ExomeDepth.exome_depth_rolling_ratios_wig
    File? exome_depth_rolling_ratios_wig_gz_tbi = ExomeDepth.exome_depth_rolling_ratios_wig_gz_tbi
    File? exome_depth_ratios_clean_wig_gz = ExomeDepth.exome_depth_ratios_clean_wig_gz
    File? exome_depth_ratios_clean_wig_gz_tbi = ExomeDepth.exome_depth_ratios_clean_wig_gz_tbi
    File? exome_depth_ratios_nomissing_wig_gz = ExomeDepth.exome_depth_ratios_nomissing_wig_gz
    File? exome_depth_ratios_nomissing_wig_gz_tbi = ExomeDepth.exome_depth_ratios_nomissing_wig_gz_tbi
    File? exome_depth_annotSV_tsv = ExomeDepth.exome_depth_annotSV_tsv
  }
}
