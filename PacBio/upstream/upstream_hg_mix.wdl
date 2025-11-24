version 1.0

# T.M.
# hg mix means calling PB upstream workflow twice with different HG reference data:
# calling pbmm2, mosdepth, deep variant and sawfish with hg19 reference data using 
# calling pbmm2, paraphase, mitorsaw and trgt tasks with hg38 reference data
# cherry picking the result files
# converting bam files to cram

import "./upstream.wdl" as upstream_hg38
import "./upstream_my_version.wdl" as upstream_hg19
# importing structure "RuntimeAttributes"
import "https://raw.githubusercontent.com/PacificBiosciences/wdl-common/1b8bbbcaf6f8783189c1ca1421f5ea94ca0f10c4/wdl/structs.wdl"

import "../../VEP/Vep2.wdl" as VEP
import "../../CRAM_conversions.wdl" as CramConversions


workflow PB_upstream {
  input {
    String sample_id
    String? sex
    Array[File] hifi_reads
    # Array[File]? fail_reads

    File hg19_ref_map_file
    File hg38_ref_map_file

    Int max_reads_per_alignment_chunk

    File trgt_catalog

    #File? fail_reads_bed
    #File? fail_reads_bait_fasta
    #File? fail_reads_bait_index

    #Boolean gpu = false

    RuntimeAttributes default_runtime_attributes
  }


  # Part HG38 
  call upstream_hg38.upstream as upstream_hg38 {
        input:
          sample_id = sample_id,
          sex = sex,
          hifi_reads = hifi_reads,
          ref_map_file = hg38_ref_map_file,
          max_reads_per_alignment_chunk = max_reads_per_alignment_chunk,
          trgt_catalog = trgt_catalog,
          single_sample = true,
          gpu = false,
          default_runtime_attributes = default_runtime_attributes
  }

  # Part HG19 
  call upstream_hg19.upstream as upstream_hg19 {
        input:
          sample_id = sample_id,
          sex = sex,
          hifi_reads = hifi_reads,
          ref_map_file = hg19_ref_map_file,
          max_reads_per_alignment_chunk = max_reads_per_alignment_chunk,
          trgt_catalog = trgt_catalog,
          single_sample = true,
          gpu = false,
          default_runtime_attributes = default_runtime_attributes
  }

  # let's rename the output files to our like
  call Rename_files {
        input:
          sample_id = sample_id,
          aligned_bam              = upstream_hg19.out_bam,
          aligned_bam_index        = upstream_hg19.out_bam_index,
          small_variant_vcf        = upstream_hg19.small_variant_vcf,
          small_variant_vcf_index  = upstream_hg19.small_variant_vcf_index,
          small_variant_gvcf       = upstream_hg19.small_variant_gvcf,
          small_variant_gvcf_index = upstream_hg19.small_variant_gvcf_index          
  }

  call VEP.VEP as VEPDeepVariant {
      input:
        sample_basename = sample_id,
        #input_vcf = RunDeepVariant.outputVCF,
        #filename_infix = ".DeepVariant"
        input_vcf = Rename_files.output_small_variant_vcf,
        filename_infix = ".DeepVariant"
  }

  Map[String, String] hg19_ref_map = read_map(hg19_ref_map_file)
  call CramConversions.ConvertToCram as ConvertToCram {
	    input:
	      input_bam = Rename_files.output_bam,
	      ref_fasta = hg19_ref_map["fasta"],
	      ref_fasta_index = hg19_ref_map["fasta_index"],
	      sample_basename = sample_id
  }

  output {
    # alignments
    File out_bam       = upstream_hg19.out_bam
    File out_bam_index = upstream_hg19.out_bam_index

    # mosdepth outputs
    File mosdepth_summary                 = upstream_hg19.mosdepth_summary
    File mosdepth_region_bed              = upstream_hg19.mosdepth_region_bed
    File mosdepth_region_bed_index        = upstream_hg19.mosdepth_region_bed_index
    File mosdepth_depth_distribution_plot = upstream_hg19.mosdepth_depth_distribution_plot
    ##String inferred_sex                 = upstream_hg19.inferred_sex
    ##String stat_depth_mean              = upstream_hg19.stat_depth_mean

    # per sample sv signatures
    ### not needed: File discover_tar = upstream_hg19.discover_tar

    # sawfish outputs for single sample
    File? sv_vcf                        = upstream_hg19.sv_vcf
    File? sv_vcf_index                  = upstream_hg19.sv_vcf_index
    File? sv_supporting_reads           = upstream_hg19.sv_supporting_reads
    File? sv_copynum_bedgraph           = upstream_hg19.sv_copynum_bedgraph
    File? sv_depth_bw                   = upstream_hg19.sv_depth_bw
    File? sv_gc_bias_corrected_depth_bw = upstream_hg19.sv_gc_bias_corrected_depth_bw
    File? sv_maf_bw                     = upstream_hg19.sv_maf_bw
    File? sv_copynum_summary            = upstream_hg19.sv_copynum_summary

    # small variant outputs
    File small_variant_vcf        = upstream_hg19.small_variant_vcf
    File small_variant_vcf_index  = upstream_hg19.small_variant_vcf_index
    File small_variant_gvcf       = upstream_hg19.small_variant_gvcf
    File small_variant_gvcf_index = upstream_hg19.small_variant_gvcf_index

    # trgt outputs
    File   trgt_vcf                  = upstream_hg38.trgt_vcf
    File   trgt_vcf_index            = upstream_hg38.trgt_vcf_index
    File   trgt_spanning_reads       = upstream_hg38.trgt_spanning_reads
    File   trgt_spanning_reads_index = upstream_hg38.trgt_spanning_reads_index
    #String stat_trgt_genotyped_count = upstream_hg38.stat_trgt_genotyped_count
    #String stat_trgt_uncalled_count  = upstream_hg38.stat_trgt_uncalled_count

    # paraphase outputs
    File? paraphase_output_json         = upstream_hg38.paraphase_output_json
    File? paraphase_realigned_bam       = upstream_hg38.paraphase_realigned_bam
    File? paraphase_realigned_bam_index = upstream_hg38.paraphase_realigned_bam_index
    File? paraphase_vcfs                = upstream_hg38.paraphase_vcfs

    # per sample mitorsaw outputs
    File mitorsaw_vcf       = upstream_hg38.mitorsaw_vcf
    File mitorsaw_vcf_index = upstream_hg38.mitorsaw_vcf_index
    File mitorsaw_hap_stats = upstream_hg38.mitorsaw_hap_stats

  }
}

# renaming some of the upstream output files so that they match our naming rules:
# PX15843.DeepVariant.vcf.gz
# PX15843.DeepVariant.vcf.gz.tbi
# PX15843.DeepVariant.VEP.hg19.annotated.vcf.gz
# PX15843.DeepVariant.VEP.hg19.annotated.vcf.gz.tbi
# bam: PX19220.PX19220.reset.hg19.aligned.bam --> PX19220.hg19.aligned.bam

task Rename_files {
  input {
    String sample_id
    File aligned_bam
    File aligned_bam_index
    File small_variant_vcf
    File small_variant_vcf_index
    File small_variant_gvcf
    File small_variant_gvcf_index
  }

  command {
    cp  ~{aligned_bam} ~{sample_id}.aligned.bam
    cp  ~{aligned_bam_index} ~{sample_id}.aligned.bam.bai
    cp ~{small_variant_vcf} ~{sample_id}.DeepVariant.vcf.gz
    cp ~{small_variant_vcf_index} ~{sample_id}.DeepVariant.vcf.gz.tbi
    cp ~{small_variant_gvcf} ~{sample_id}.DeepVariant.gvcf.gz
    cp ~{small_variant_gvcf_index} ~{sample_id}.DeepVariant.gvcf.tbi
  }

  runtime {
    docker: "alpine:latest"
    cpu: 1
    requested_memory_mb_per_core: 1000
    runtime_minutes: 5
  }
  output {
    File output_bam = "~{sample_id}.aligned.bam"
    File output_bam_index = "~{sample_id}.aligned.bam.bai"
    File output_small_variant_vcf "~{sample_id}.DeepVariant.vcf.gz"
    File output_small_variant_vcf_index "~{sample_id}.DeepVariant.vcf.gz.tbi"
    File output_small_variant_gvcf "~{sample_id}.DeepVariant.gvcf.gz"
    File output_small_variant_gvcf_index "~{sample_id}.DeepVariant.gvcf.tbi"
  }
}

