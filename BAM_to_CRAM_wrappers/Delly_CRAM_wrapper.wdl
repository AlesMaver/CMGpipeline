version 1.0
## Copyright CMG@KIGM, Ales Maver, Tanja Majnik

# Wrapping WF 
# using old wf with BAM file input, but having a CRAM file for the input
# and calling CRAM to BAM file conversion WF

import "../CRAM_conversions.wdl" as CramConversions 
import "../Delly/DELLY_single3.wdl" as Delly

# WORKFLOW DEFINITION 
workflow Delly_CRAM_wrapper {
  input {
    String sample_basename
    File? input_cram
    File? input_cram_index
    File reference_fa
    File reference_fai
    File reference_dict
    File population_bcf
    File population_bcf_index
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

  call Delly.DELLY as Delly {
      input:
        input_bam = CramToBam.output_bam,
        input_bai = CramToBam.output_bai,
        population_bcf = population_bcf,
        population_bcf_index = population_bcf_index,
        reference_fasta = reference_fa
  }

  output {
    File call_bcf_file = Delly.call_bcf_file
    File sample_bcfs = Delly.sample_bcfs
    File sample_bcf_indices = Delly.sample_bcf_indices
    File unfiltered_population_bcf = Delly.unfiltered_population_bcf
    File filtered_population_bcf = Delly.filtered_population_bcf
    File? delly_annotSV = Delly.delly_annotSV
  }
}
