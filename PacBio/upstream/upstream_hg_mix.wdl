version 1.0

# T.M.
# hg mix means calling PB upstream workflow twice with different HG reference data:
# calling pbmm2, mosdepth, deep variant and sawfish with hg19 reference data using 
# calling pbmm2, paraphase, mitorsaw and trgt tasks with hg38 reference data
# cherry picking the result files
# converting bam files to cram

import "./upstream.wdl" as upstream_hg38
import "./upstream_my_version.wdl" as upstream_hg19


workflow PB_upstream {
  meta {
    description: "Given a set of HiFi reads for a human sample, run steps upstream of phasing."
  }

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
          runtime_attributes = default_runtime_attributes
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
          runtime_attributes = default_runtime_attributes
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
    File discover_tar = upstream_hg19.discover_tar

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
    String stat_trgt_genotyped_count = upstream_hg38.stat_trgt_genotyped_count
    String stat_trgt_uncalled_count  = upstream_hg38.stat_trgt_uncalled_count

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
