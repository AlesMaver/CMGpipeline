version 1.0

import "./CRAM_conversions.wdl" as CramConversions
import "./manta/manta_workflow.wdl" as Manta

workflow ExomeDepth {
  input {
    String sample_name
    File? target_bed
    File? input_bam
    File? input_bam_index
    File? input_cram
    File? input_cram_index
    File? reference_fa
    File? reference_fai
    File? reference_dict
    File? exome_depth_counts_input
    Array[File]? reference_counts_files
  }

  if(!defined(exome_depth_counts_input)) {

    if (defined(input_cram)) {
      call CramConversions.CramToBam as CramToBam {
        input:
          sample_name = sample_name,
          input_cram = input_cram,
          ref_fasta = reference_fa,
          ref_fasta_index = reference_fai,
          ref_dict = reference_dict,
          docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817",
          samtools_path = "samtools"
      }
    }

    call GetCounts {
      input:
          sample_name = sample_name,
          target_bed = target_bed,
          # input_bam = select_first([input_bam]),
          # input_bam_index = select_first([input_bam_index])
          input_bam = select_first([input_bam, CramToBam.output_bam]),
          input_bam_index = select_first([input_bam_index, CramToBam.output_bai])
    }
  }

  if(defined(reference_counts_files)) {
    call ExomeDepth {
      input:
        sample_name = sample_name,
        target_bed = target_bed,
        test_counts_file = select_first([exome_depth_counts_input, GetCounts.exome_depth_counts]),
        reference_counts_files = select_first([reference_counts_files])
    }

    call Manta.annotSV as annotSV {
      input:
        genome_build = "GRCh37",
        input_vcf = ExomeDepth.exome_depth_cnv_calls_csv,
        output_tsv_name = sample_name + ".ExomeDepth.annotSV.tsv"
    }
  }

  output {
    File? exome_depth_counts = GetCounts.exome_depth_counts
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
    File? exome_depth_annotSV_tsv = annotSV.sv_variants_tsv
  }
}


task GetCounts {
  input {
    String sample_name
    File? target_bed
    File input_bam
    File input_bam_index
  }

  command <<<
    wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/Dockers/ExomeDepth/ExomeDepth.R

    Rscript ExomeDepth.R \
        --targets ~{target_bed} \
        --test-sample-bam ~{input_bam}
  >>>

  runtime {
        docker: "alesmaver/exome_depth"
        maxRetries: 3
        # requested_memory_mb_per_core: 3000
        # cpu: 6
        requested_memory_mb_per_core: 1000
        cpu: 16
        runtime_minutes: 180
    }

  output {
        File exome_depth_counts = "~{sample_name}_ExomeDepth_counts.tsv"
  }
}

task ExomeDepth {
  input {
    String sample_name
    File? target_bed
    File test_counts_file
    Array[File] reference_counts_files
  }

  command <<<
    wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/Dockers/ExomeDepth/ExomeDepth.R

    Rscript ExomeDepth.R \
        --test-counts-file ~{test_counts_file} \
        --reference-counts-file-list ~{write_lines(reference_counts_files)} \
        --targets ~{target_bed}
    >>>

    runtime {
        docker: "alesmaver/exome_depth"
        maxRetries: 3
        requested_memory_mb_per_core: 1000
        cpu: 6
        runtime_minutes: 180
        }
        
    output {
      File exome_depth_cnv_calls_bed = "~{sample_name}_ExomeDepth_CNV.bed"
      File exome_depth_cnv_calls_csv = "~{sample_name}_ExomeDepth_CNV_annotSV.bed"

      File exome_depth_ratios_all_wig_gz = "~{sample_name}_ExomeDepth_ratios_all.wig.gz"
      File exome_depth_ratios_all_wig_gz_tbi = "~{sample_name}_ExomeDepth_ratios_all.wig.gz.tbi"

      File exome_depth_rolling_ratios_wig = "~{sample_name}_ExomeDepth_rolling_ratios.wig.gz"
      File exome_depth_rolling_ratios_wig_gz_tbi = "~{sample_name}_ExomeDepth_rolling_ratios.wig.gz.tbi"

      File exome_depth_ratios_clean_wig_gz = "~{sample_name}_ExomeDepth_ratios_clean.wig.gz"
      File exome_depth_ratios_clean_wig_gz_tbi = "~{sample_name}_ExomeDepth_ratios_clean.wig.gz.tbi"

      File exome_depth_ratios_nomissing_wig_gz = "~{sample_name}_ExomeDepth_ratios_nomissing.wig.gz"
      File exome_depth_ratios_nomissing_wig_gz_tbi = "~{sample_name}_ExomeDepth_ratios_nomissing.wig.gz.tbi"
    }
}
