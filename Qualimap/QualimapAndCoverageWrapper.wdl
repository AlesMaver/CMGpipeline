version 1.0
## Copyright CMG@KIGM, Ales Maver, TM


# Subworkflows
import "../CRAM_conversions.wdl" as CramConversions
import "../PrepareMaskedGenomeFiles.wdl" as PrepareMaskedGenomeFileswdl
import "./QualimapAndCoverage.wdl" as QualimapAndCoverageWdl

# WORKFLOW DEFINITION 
workflow QualimapAndCoverageWrapper {
  input {
    String sample_basename
    File? input_bam
    File? input_bam_index
    File? input_cram
    File? input_cram_index

    File reference_fixed_fa
    File reference_fixed_fai

    File reference_fa
    File reference_fai
    File reference_dict

    String? enrichment
    File? enrichment_bed

    File refSeqFile

    String? targetRegions
    Boolean? perform_masked_alignment

    Int threads
  }

  # START

  if (defined(input_cram)) {
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
  }

  
  # 
  if ( defined(targetRegions) && select_first([perform_masked_alignment, false]) ) {
    call PrepareMaskedGenomeFileswdl.PrepareMaskedGenomeFasta as PrepareMaskedGenomeFasta {
      input:
        reference_fixed_fa=reference_fixed_fa,
        reference_fixed_fai=reference_fixed_fai,
        targetRegions=targetRegions,

        # Runtime 
        docker = "pegi3s/bedtools"
    }
  }
  
  call QualimapAndCoverageWdl.QualimapAndCoverage as QualimapAndCoverage {
       input:
          sample_basename = sample_basename,
          input_bam = select_first([input_bam, CramToBam.output_bam]),
          input_bam_index = select_first([input_bam_index, CramToBam.output_bai]),
          reference_fa = reference_fa,
          reference_fai = reference_fai,
          reference_dict = reference_dict,
          enrichment = enrichment,
          enrichment_bed = enrichment_bed,
          refSeqFile = refSeqFile,
          targetRegions = targetRegions,
          maskedGenomeFastaTargetRegions_bed = PrepareMaskedGenomeFasta.targetRegions_bed,
          threads = threads
  }

  output {
    File? Qualimap_results = QualimapAndCoverage.Qualimap_results
    File? QualimapWGS_results = QualimapAndCoverage.QualimapWGS_results
    File? DepthOfCoverage_output = QualimapAndCoverage.DepthOfCoverage_output
    File? DepthOfCoverageWGS_output = QualimapAndCoverage.DepthOfCoverageWGS_output
  }

}

