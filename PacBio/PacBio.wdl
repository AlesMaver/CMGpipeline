version 1.0

#import "https://raw.githubusercontent.com/AlesMaver/HiFi-human-WGS-WDL/refs/heads/main/workflows/upstream/upstream.wdl" as Upstream
import "http://10.3.248.96:3500/WDL/HiFi-human-WGS-WDL/workflows/upstream/upstream.wdl" as Upstream

import "../VEP/Vep2.wdl" as VEP
import "../CRAM_conversions.wdl" as CramConversions
import "../AnnotationPipeline.wdl" as Annotation

workflow PacBioWorkflow {
  meta {
    description: "PacBio HiFi upstream workflow with VEP annotation and custom annotation pipeline"
    author: "Ales Maver"
  }

  input {
    # ========================================
    # PacBio Upstream Workflow Inputs
    # ========================================
    String sample_id
    String? sex
    Array[File] hifi_reads
    Array[File]? fail_reads
    File ref_map_file
    Int max_reads_per_alignment_chunk = 5000000
    File trgt_catalog
    File? fail_reads_bed
    File? fail_reads_bait_fasta
    File? fail_reads_bait_index
    File? regions_bed
    Boolean single_sample = true
    Boolean gpu = false
    RuntimeAttributes default_runtime_attributes
    Boolean run_sawfish = true

    # ========================================
    # Post-Processing Options
    # ========================================
    Boolean run_vep_annotation = true
    Boolean run_custom_annotation = true
    Boolean convert_to_cram = true

    # ========================================
    # VEP Annotation Inputs
    # ========================================
    String vep_filename_infix = ""

    # ========================================
    # Custom Annotation Pipeline Inputs
    # ========================================
    File? chromosome_list
    File? gnomAD_vcf
    File? gnomAD_vcf_index
    File? gnomADexomes_vcf
    File? gnomADexomes_vcf_index
    File? TopMed_vcf
    File? TopMed_vcf_index
    File? Regeneron_vcf
    File? Regeneron_vcf_index
    File? GnomAD4_exomes_vcf
    File? GnomAD4_exomes_vcf_index
    File? GnomAD4_genomes_vcf
    File? GnomAD4_genomes_vcf_index
    File? blacklisted_regions_bed
    File? blacklisted_regions_bed_index
    File? GERP_bed
    File? GERP_bed_index
    File? SLOpopulation_vcf
    File? SLOpopulation_vcf_index
    File? ClinVar_vcf
    File? ClinVar_vcf_index
    File? SpliceAI
    File? SpliceAI_index
    File? dbscSNV
    File? dbscSNV_index
    File? HPO
    File? HPO_index
    File? OMIM
    File? OMIM_index
    File? gnomadConstraints
    File? gnomadConstraints_index
    File? CGD
    File? CGD_index
    File? bcftools_annotation_header
    File? gnomAD_pLi_bed
    File? gnomAD_pLi_bed_index
    File? gnomAD_misz_bed
    File? gnomAD_misz_bed_index
    File? metaDome_bed
    File? metaDome_bed_index
    File? pext_bed
    File? pext_bed_index
    File? dbNSFP
    File? dbNSFP_index
  }

  # Read reference map file for CRAM and annotation
  Map[String, String] ref_map = read_map(ref_map_file)

  # ========================================
  # Run PacBio Upstream Workflow
  # ========================================
  call Upstream.upstream as pacbio_upstream {
    input:
      sample_id                      = sample_id,
      sex                            = sex,
      hifi_reads                     = hifi_reads,
      fail_reads                     = fail_reads,
      ref_map_file                   = ref_map_file,
      max_reads_per_alignment_chunk  = max_reads_per_alignment_chunk,
      trgt_catalog                   = trgt_catalog,
      fail_reads_bed                 = fail_reads_bed,
      fail_reads_bait_fasta          = fail_reads_bait_fasta,
      fail_reads_bait_index          = fail_reads_bait_index,
      regions_bed                    = regions_bed,
      single_sample                  = single_sample,
      gpu                            = gpu,
      run_sawfish                    = run_sawfish,
      default_runtime_attributes     = default_runtime_attributes
  }

  # ========================================
  # Convert BAM to CRAM (Optional)
  # ========================================
  if (convert_to_cram) {
    call CramConversions.ConvertToCram as ConvertToCram {
      input:
        input_bam             = pacbio_upstream.out_bam,
        ref_fasta       = ref_map["fasta"],
        ref_fasta_index = ref_map["fasta_index"],
        sample_basename       = sample_id
    }
  }

  # ========================================
  # VEP Annotation on Small Variants (Optional)
  # ========================================
  if (run_vep_annotation) {
    call VEP.VEP as VEPDeepVariant {
      input:
        sample_basename = sample_id,
        input_vcf       = pacbio_upstream.small_variant_vcf,
        filename_infix  = vep_filename_infix
    }
  }

  # ========================================
  # Custom Annotation on Small Variants (Optional)
  # ========================================
  if (run_custom_annotation) {
    
    call Annotation.AnnotateVCF as AnnotateSmallVariants {
      input:
        input_vcf                      = pacbio_upstream.small_variant_vcf,
        sample_basename                = sample_id,
        chromosome_list                = select_first([chromosome_list]),
        gnomAD_vcf                     = select_first([gnomAD_vcf]),
        gnomAD_vcf_index               = select_first([gnomAD_vcf_index]),
        gnomADexomes_vcf               = select_first([gnomADexomes_vcf]),
        gnomADexomes_vcf_index         = select_first([gnomADexomes_vcf_index]),
        TopMed_vcf                     = select_first([TopMed_vcf]),
        TopMed_vcf_index               = select_first([TopMed_vcf_index]),
        Regeneron_vcf                  = select_first([Regeneron_vcf]),
        Regeneron_vcf_index            = select_first([Regeneron_vcf_index]),
        GnomAD4_exomes_vcf             = select_first([GnomAD4_exomes_vcf]),
        GnomAD4_exomes_vcf_index       = select_first([GnomAD4_exomes_vcf_index]),
        GnomAD4_genomes_vcf            = select_first([GnomAD4_genomes_vcf]),
        GnomAD4_genomes_vcf_index      = select_first([GnomAD4_genomes_vcf_index]),
        blacklisted_regions_bed        = select_first([blacklisted_regions_bed]),
        blacklisted_regions_bed_index  = select_first([blacklisted_regions_bed_index]),
        GERP_bed                       = select_first([GERP_bed]),
        GERP_bed_index                 = select_first([GERP_bed_index]),
        SLOpopulation_vcf              = select_first([SLOpopulation_vcf]),
        SLOpopulation_vcf_index        = select_first([SLOpopulation_vcf_index]),
        ClinVar_vcf                    = select_first([ClinVar_vcf]),
        ClinVar_vcf_index              = select_first([ClinVar_vcf_index]),
        SpliceAI                       = select_first([SpliceAI]),
        SpliceAI_index                 = select_first([SpliceAI_index]),
        dbscSNV                        = select_first([dbscSNV]),
        dbscSNV_index                  = select_first([dbscSNV_index]),
        HPO                            = select_first([HPO]),
        HPO_index                      = select_first([HPO_index]),
        OMIM                           = select_first([OMIM]),
        OMIM_index                     = select_first([OMIM_index]),
        gnomadConstraints              = select_first([gnomadConstraints]),
        gnomadConstraints_index        = select_first([gnomadConstraints_index]),
        CGD                            = select_first([CGD]),
        CGD_index                      = select_first([CGD_index]),
        bcftools_annotation_header     = select_first([bcftools_annotation_header]),
        gnomAD_pLi_bed                 = select_first([gnomAD_pLi_bed]),
        gnomAD_pLi_bed_index           = select_first([gnomAD_pLi_bed_index]),
        gnomAD_misz_bed                = select_first([gnomAD_misz_bed]),
        gnomAD_misz_bed_index          = select_first([gnomAD_misz_bed_index]),
        metaDome_bed                   = select_first([metaDome_bed]),
        metaDome_bed_index             = select_first([metaDome_bed_index]),
        pext_bed                       = select_first([pext_bed]),
        pext_bed_index                 = select_first([pext_bed_index]),
        fasta_reference                = ref_map["fasta"],
        fasta_reference_index          = ref_map["fasta_index"],
        fasta_reference_dict           = ref_map["fasta_dict"],
        dbNSFP                         = select_first([dbNSFP]),
        dbNSFP_index                   = select_first([dbNSFP_index])
    }
  }

  output {
    # ========================================
    # Alignment Outputs
    # ========================================
    File out_bam                = pacbio_upstream.out_bam
    File out_bam_index          = pacbio_upstream.out_bam_index
    File? out_cram              = ConvertToCram.output_cram
    File? out_cram_index        = ConvertToCram.output_cram_index

    # ========================================
    # Quality Control Outputs
    # ========================================
    File   mosdepth_summary                 = pacbio_upstream.mosdepth_summary
    File   mosdepth_region_bed              = pacbio_upstream.mosdepth_region_bed
    File   mosdepth_region_bed_index        = pacbio_upstream.mosdepth_region_bed_index
    File   mosdepth_depth_distribution_plot = pacbio_upstream.mosdepth_depth_distribution_plot
    String inferred_sex                     = pacbio_upstream.inferred_sex
    String stat_depth_mean                  = pacbio_upstream.stat_depth_mean

    # ========================================
    # Small Variant Outputs
    # ========================================
    File small_variant_vcf              = pacbio_upstream.small_variant_vcf
    File small_variant_vcf_index        = pacbio_upstream.small_variant_vcf_index
    File small_variant_gvcf             = pacbio_upstream.small_variant_gvcf
    File small_variant_gvcf_index       = pacbio_upstream.small_variant_gvcf_index
    File? small_variant_vcf_vep         = VEPDeepVariant.output_vcf
    File? small_variant_vcf_vep_index   = VEPDeepVariant.output_vcf_index
    File? small_variant_vcf_annotated   = AnnotateSmallVariants.output_vcf

    # ========================================
    # Structural Variant Outputs
    # ========================================
    File  discover_tar                  = pacbio_upstream.discover_tar
    File? sv_vcf                        = pacbio_upstream.sv_vcf
    File? sv_vcf_index                  = pacbio_upstream.sv_vcf_index
    File? sv_supporting_reads           = pacbio_upstream.sv_supporting_reads
    File? sv_copynum_bedgraph           = pacbio_upstream.sv_copynum_bedgraph
    File? sv_depth_bw                   = pacbio_upstream.sv_depth_bw
    File? sv_gc_bias_corrected_depth_bw = pacbio_upstream.sv_gc_bias_corrected_depth_bw
    File? sv_maf_bw                     = pacbio_upstream.sv_maf_bw
    File? sv_copynum_summary            = pacbio_upstream.sv_copynum_summary

    # ========================================
    # Tandem Repeat Genotyping (TRGT) Outputs
    # ========================================
    File   trgt_vcf                  = pacbio_upstream.trgt_vcf
    File   trgt_vcf_index            = pacbio_upstream.trgt_vcf_index
    File   trgt_spanning_reads       = pacbio_upstream.trgt_spanning_reads
    File   trgt_spanning_reads_index = pacbio_upstream.trgt_spanning_reads_index
    String stat_trgt_genotyped_count = pacbio_upstream.stat_trgt_genotyped_count
    String stat_trgt_uncalled_count  = pacbio_upstream.stat_trgt_uncalled_count

    # ========================================
    # Paraphase Outputs (Pharmacogenes)
    # ========================================
    File? paraphase_output_json         = pacbio_upstream.paraphase_output_json
    File? paraphase_realigned_bam       = pacbio_upstream.paraphase_realigned_bam
    File? paraphase_realigned_bam_index = pacbio_upstream.paraphase_realigned_bam_index
    File? paraphase_vcfs                = pacbio_upstream.paraphase_vcfs

    # ========================================
    # Mitochondrial Variant Outputs
    # ========================================
    File mitorsaw_vcf       = pacbio_upstream.mitorsaw_vcf
    File mitorsaw_vcf_index = pacbio_upstream.mitorsaw_vcf_index
    File mitorsaw_hap_stats = pacbio_upstream.mitorsaw_hap_stats

    # ========================================
    # QC Messages
    # ========================================
    Array[String] qc_messages = pacbio_upstream.msg
  }
}
