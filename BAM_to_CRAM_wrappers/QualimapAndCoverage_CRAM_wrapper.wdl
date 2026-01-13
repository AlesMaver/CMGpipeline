version 1.0
## Copyright CMG@KIGM, Ales Maver, Tanja Majnik

# Wrapping WF 
# using old wf with BAM file input, but having a CRAM file for the input
# and calling CRAM to BAM file conversion WF

import "../CRAM_conversions.wdl" as CramConversions 
import "../Qualimap/QualimapAndCoverageWrapper.wdl" as QualimapAndCoverageWrapper

# WORKFLOW DEFINITION 
workflow QualimapAndCoverageWrapper_CRAM_wrapper {
  input {
    String sample_basename
    File? input_cram
    File? input_cram_index
    File reference_fa
    File reference_fai
    File reference_dict

    File reference_fixed_fa
    File reference_fixed_fai
    String? enrichment
    File? enrichment_bed
    File refSeqFile
    String? targetRegions
    Boolean? perform_masked_alignment
    Int threads

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

  call QualimapAndCoverageWrapper.QualimapAndCoverageWrapper as QualimapAndCoverageWrapper {
      input:
        sample_basename = sample_basename,
        input_bam = CramToBam.output_bam,
        input_bam_index = CramToBam.output_bai,
        reference_fa = reference_fa,
        reference_fai = reference_fai,
        reference_dict = reference_dict,
        reference_fixed_fa = reference_fixed_fa,        
        reference_fixed_fai = reference_fixed_fai,
        enrichment = enrichment,
        enrichment_bed = enrichment_bed,
        refSeqFile = refSeqFile,
        targetRegions = targetRegions,
        perform_masked_alignment = perform_masked_alignment,
        threads = threads
  }

  output {
    File? Qualimap_results = QualimapAndCoverageWrapper.Qualimap_results
    File? QualimapWGS_results = QualimapAndCoverageWrapper.QualimapWGS_results
    File? DepthOfCoverage_output = QualimapAndCoverageWrapper.DepthOfCoverage_output
    File? DepthOfCoverageWGS_output = QualimapAndCoverageWrapper.DepthOfCoverageWGS_output
  }
}
