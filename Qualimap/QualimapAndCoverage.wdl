version 1.0
## Copyright CMG@KIGM, Ales Maver, TM


# Subworkflows
import "../Qualimap.wdl" as Qualimap

# WORKFLOW DEFINITION 
workflow QualimapAndCoverage {
  input {
    String sample_basename
    File input_bam
    File input_bam_index

    File reference_fa
    File reference_fai
    File reference_dict

    String? enrichment
    File? enrichment_bed

    File refSeqFile

    String? targetRegions
    File? maskedGenomeFastaTargetRegions_bed

    Int threads
  }

  # START
  
  if ( defined(targetRegions) ){
    call RegionsToBed {
      input:
        targetRegions=select_first([targetRegions,""]),

        # Runtime 
        docker = "pegi3s/bedtools"

    }
  }

  if( defined(enrichment_bed) || defined(maskedGenomeFastaTargetRegions_bed) || defined(RegionsToBed.targetRegions_bed) ){
    call Qualimap.bamqc as Qualimap {
    input:
      bam = input_bam,
      sample_basename = sample_basename,
      
      enrichment_bed = select_first([maskedGenomeFastaTargetRegions_bed, RegionsToBed.targetRegions_bed, enrichment_bed]),

      ncpu = 8
    }
  }

  if ( defined(enrichment) ){
    if( enrichment=="WGS1Mb" ){
       call Qualimap.DownsampleBED as DownsampleBED {
        input:
           bed_file = enrichment_bed,
           reference_fai = reference_fai
      }
    }
  }

  # Depth of coverage
  if( defined(enrichment_bed) || defined(maskedGenomeFastaTargetRegions_bed) || defined(RegionsToBed.targetRegions_bed) ){
    call Qualimap.DepthOfCoverage34 as DepthOfCoverage {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        sample_basename = sample_basename,

        reference_fa=reference_fa,
        reference_fai=reference_fai,
        reference_dict=reference_dict,

        enrichment_bed = select_first([maskedGenomeFastaTargetRegions_bed, RegionsToBed.targetRegions_bed, DownsampleBED.downsampled_bed_file, enrichment_bed]),

        refSeqFile = refSeqFile,

        threads = threads,
        docker = "broadinstitute/gatk3:3.8-1",
        gatk_path = "/usr/GenomeAnalysisTK.jar"
    }
  }
  
  # Calculate WGS coverage if neither enrichment_bed nor target regions parameter is defined
  if( !defined(enrichment_bed) && !defined(maskedGenomeFastaTargetRegions_bed) && !defined(RegionsToBed.targetRegions_bed) ){
  
    call Qualimap.bamqc as QualimapWGS {
      input:
        bam = input_bam,
        sample_basename = sample_basename,
        ncpu = 8
    }
  
    call Qualimap.DepthOfCoverage34 as DepthOfCoverageWGS {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        sample_basename = sample_basename,

        reference_fa=reference_fa,
        reference_fai=reference_fai,
        reference_dict=reference_dict,

        refSeqFile = refSeqFile,

        threads = threads,
        docker = "broadinstitute/gatk3:3.8-1",
        gatk_path = "/usr/GenomeAnalysisTK.jar"
    }
  }

  output {
    File? Qualimap_results = Qualimap.results
    File? QualimapWGS_results = QualimapWGS.results

    File? DepthOfCoverage_output = DepthOfCoverage.DepthOfCoverage_output
    File? DepthOfCoverageWGS_output = DepthOfCoverageWGS.DepthOfCoverage_output
  }

}

# -----
# TASKS
# -----

task RegionsToBed {
  input {
    # Command parameters
    String? targetRegions # FORMAT "chr1:123033-130000;chrX:1-1000"

    # Runtime parameters
    String docker
  }
  
  command <<<
    set -e

    echo "~{targetRegions}"  | tr ';' '\n' | tr ':' '\t' | tr '-' '\t' > targetRegions.bed
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 500
    cpu: 1
    runtime_minutes: 10
  }
  output {
    File targetRegions_bed = "targetRegions.bed"
  }
}






