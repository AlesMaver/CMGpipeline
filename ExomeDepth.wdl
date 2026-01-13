version 1.0

import "./CRAM_conversions.wdl" as CramConversions
import "./manta/manta_workflow.wdl" as Manta

workflow ExomeDepth {
  input {
    String sample_name
    String test_type = "WES"
    String? enrichment
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

  if (defined(enrichment)) {
    String test_type = if (enrichment == "WGS1Mb") then "WGS" else "WES"
  }

  # If the input has not yet been counted, convert the CRAM to BAM and count the reads
  if(!defined(exome_depth_counts_input) && test_type == "WES") {

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

  if(!defined(exome_depth_counts_input) && test_type == "WGS") {
    call GetCounts_Bedtools {
      input:
          sample_name = sample_name,
          target_bed = target_bed,
          input_bam = select_first([input_cram, input_bam]),
          input_bam_index = select_first([input_cram_index, input_bam_index]),
          reference_fa = reference_fa,
          reference_fai = reference_fai,
          reference_dict = reference_dict
      }
    }

  ## Select the output from either GetCounts or GetCounts_Bedtools
  # old: File exome_depth_counts_calculated = select_first([GetCounts.exome_depth_counts, GetCounts_Bedtools.exome_depth_counts])

  if(defined(reference_counts_files)) {
    call ExomeDepth {
      input:
        sample_name = sample_name,
        target_bed = target_bed,
        # old: test_counts_file = select_first([exome_depth_counts_input, exome_depth_counts_calculated]),
        test_counts_file = select_first([exome_depth_counts_input, GetCounts.exome_depth_counts, GetCounts_Bedtools.exome_depth_counts]),
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
    ## File? exome_depth_counts = exome_depth_counts_calculated
    File? exome_depth_counts = select_first([exome_depth_counts_input, GetCounts.exome_depth_counts, GetCounts_Bedtools.exome_depth_counts])
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

task GetCounts_Bedtools {
  input {
    String sample_name
    File? target_bed
    File input_bam
    File input_bam_index
    File? reference_fa
    File? reference_fai
    File? reference_dict
    Int samtools_threads = 30
  }

  command <<<
    # Get BAM chromosome order
    samtools view -H ~{input_bam} | grep '@SQ' | awk '{print $2}' | cut -d':' -f2 > bam_chrom_order.txt

    # Get BED chromosome order
    cut -f1 ~{target_bed} | uniq > bed_chrom_order.txt

    # Sort BED file by chromosome order in BAM (so that bedtools will work)
    awk 'NR==FNR {order[$1]=NR; next} $1 in order {print order[$1], $0}' bam_chrom_order.txt ~{target_bed} | sort -k1,1n -k3,3n | cut -d' ' -f2- > bam_sorted.bed

    # Calculate coverage using bedtools
    samtools view -@ ~{samtools_threads} -h -F 3852 -f 3 -q 30 ~{if defined(reference_fa) then "-T " + reference_fa else ""} ~{input_bam} | \
    awk '$1 ~ /^@/ || $9 > 0' | \
    samtools view -@ ~{samtools_threads} -O BAM -h | \
    bedtools bamtobed -i stdin | \
    bedtools coverage -a bam_sorted.bed -b stdin -counts -sorted > bam_sorted_counts.bed

    awk 'NR==FNR {order[$1]=NR; next} $1 in order {print order[$1], $0}' bed_chrom_order.txt bam_sorted_counts.bed | sort -k1,1n -k3,3n | cut -d' ' -f2- > bed_sorted_counts.bed
    (echo -e "chromosome\tstart\tend\t~{sample_name}.bam"; cat bed_sorted_counts.bed) > ~{sample_name}_ExomeDepth_counts.tsv
  >>>

  runtime {
    docker: "danhumassmed/samtools-bedtools:1.0.2"
    maxRetries: 3
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
        cpu: 16
        runtime_minutes: 180
        # continue only when rc=0 or rc=1 (know problem inside the R script)
        continueOnReturnCode: [0, 1]
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
