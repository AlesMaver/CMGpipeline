version 1.0

# T.M.
# hg mix means calling PB upstream workflow twice with different HG reference data:
# calling pbmm2, mosdepth, deep variant and sawfish with hg19 reference data using 
# calling pbmm2, paraphase, mitorsaw and trgt tasks with hg38 reference data
# cherry picking the result files
# converting bam files to cram
# annotating deep variant vcf with VEP tool 
# annotating deep variant vcf with AnnotateVCF tool
# annotating sawfish sv_vcf with annotSV tool

import "https://raw.githubusercontent.com/tanmaj/PacBio-HiFi-human-WGS-WDL/refs/heads/our-general-changes/workflows/upstream/upstream_kigm.wdl" as upstream
# importing structure "RuntimeAttributes". This structure is hardly needed, taking along due to compatibility with original workflow.
import "https://raw.githubusercontent.com/tanmaj/wdl-common/refs/heads/19ca392-branch/wdl/structs.wdl"

import "../../VEP/Vep2.wdl" as VEP
import "../../CRAM_conversions.wdl" as CramConversions
import "../../manta/manta_workflow.wdl" as manta
import "../../AnnotationPipeline.wdl" as Annotation


workflow PacBioUpstream {
  input {
    # --------------------------------#
    #  inputs for PacBio upstream 
    # --------------------------------#
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

    # --------------------------------#
    #  inputs for annotateVCF 
    # --------------------------------#
    File chromosome_list
    
    File gnomAD_vcf
    File gnomAD_vcf_index

    File gnomADexomes_vcf
    File gnomADexomes_vcf_index

    File TopMed_vcf
    File TopMed_vcf_index

    File Regeneron_vcf
    File Regeneron_vcf_index

    File GnomAD4_exomes_vcf
    File GnomAD4_exomes_vcf_index

    File GnomAD4_genomes_vcf
    File GnomAD4_genomes_vcf_index

    File blacklisted_regions_bed
    File blacklisted_regions_bed_index

    File GERP_bed
    File GERP_bed_index

    File SLOpopulation_vcf
    File SLOpopulation_vcf_index

    File ClinVar_vcf
    File ClinVar_vcf_index

    File SpliceAI
    File SpliceAI_index

    File dbscSNV
    File dbscSNV_index

    File HPO
    File HPO_index
    File OMIM
    File OMIM_index
    File gnomadConstraints
    File gnomadConstraints_index
    File CGD
    File CGD_index
    File bcftools_annotation_header

    File gnomAD_pLi_bed
    File gnomAD_pLi_bed_index

    File gnomAD_misz_bed
    File gnomAD_misz_bed_index

    File metaDome_bed
    File metaDome_bed_index

    File pext_bed
    File pext_bed_index

    File fasta_reference
    File fasta_reference_index
    File fasta_reference_dict

    File dbNSFP
    File dbNSFP_index

  }


  # Part HG38 
  call upstream.upstream as upstream_hg38 {
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
  call upstream.upstream as upstream_hg19 {
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

  if ( defined(upstream_hg19.sv_vcf) ) {
    # annotSV for sawfish sv vcf file
    call manta.annotSV as annotSV {
        input:
            genome_build = "GRCh37",
            input_vcf = select_first([upstream_hg19.sv_vcf, ""]),
            output_tsv_name = sample_id + ".sawfish.AnnotSV.tsv"
    }
  }

  # let's rename the output files to our like
  call Rename_files {
        input:
          sample_id = sample_id,
          #aligned_bam              = select_first([upstream_hg19.out_bam, ""]),
          #aligned_bam_index        = select_first([upstream_hg19.out_bam_index, ""]),
          small_variant_vcf        = select_first([upstream_hg19.small_variant_vcf, ""]),
          small_variant_vcf_index  = select_first([upstream_hg19.small_variant_vcf_index, ""]),
          small_variant_gvcf       = select_first([upstream_hg19.small_variant_gvcf, ""]),
          small_variant_gvcf_index = select_first([upstream_hg19.small_variant_gvcf_index, ""]),   
          sv_vcf                        = select_first([upstream_hg19.sv_vcf, ""]),
          sv_vcf_index                  = select_first([upstream_hg19.sv_vcf_index, ""]),
          sv_supporting_reads           = select_first([upstream_hg19.sv_supporting_reads, ""]),
          sv_copynum_bedgraph           = select_first([upstream_hg19.sv_copynum_bedgraph, ""]),
          sv_depth_bw                   = select_first([upstream_hg19.sv_depth_bw, ""]),
          sv_gc_bias_corrected_depth_bw = select_first([upstream_hg19.sv_gc_bias_corrected_depth_bw, ""]),
          sv_maf_bw                     = select_first([upstream_hg19.sv_maf_bw, ""]),
          sv_copynum_summary            = select_first([upstream_hg19.sv_copynum_summary, ""]),
          mosdepth_summary                 = select_first([upstream_hg19.mosdepth_summary, ""]),
          mosdepth_region_bed              = select_first([upstream_hg19.mosdepth_region_bed, ""]),
          mosdepth_region_bed_index        = select_first([upstream_hg19.mosdepth_region_bed_index, ""]),
          mosdepth_depth_distribution_plot = select_first([upstream_hg19.mosdepth_depth_distribution_plot, ""])
  }

  call VEP.VEP as VEPDeepVariant {
      input:
        sample_basename = sample_id,
        input_vcf = Rename_files.output_small_variant_vcf,
        filename_infix = ".DeepVariant"
  }

  # transformation of upstream_hg19.small_variant_vcf file to match AnnotateVCF workflow: no mitochondria data, transcribing GRCh37 notation into hg19 notation
  call TransformVcfFile {
      input:
        output_vcf_name = sample_id + ".DeepVariant.vcf.gz",
        input_vcf = select_first([upstream_hg19.small_variant_vcf, ""])
  }

  Map[String, String] hg19_ref_map = read_map(hg19_ref_map_file)
  call CramConversions.ConvertToCram as ConvertToCram {
	    input:
	      input_bam = Rename_files.output_bam,
	      ref_fasta = hg19_ref_map["fasta"],
	      ref_fasta_index = hg19_ref_map["fasta_index"],
	      sample_basename = sample_id
  }

  call Annotation.AnnotateVCF as AnnotateVCF_DeepVariant{
      input:
        sample_basename = sample_id,
        output_filename = sample_id + ".DeepVariant.annotated.vcf",
        input_vcf = select_first([TransformVcfFile.output_vcf, ""]),

        chromosome_list = chromosome_list,
        gnomAD_vcf = gnomAD_vcf,
        gnomAD_vcf_index = gnomAD_vcf_index,

        gnomADexomes_vcf = gnomADexomes_vcf,
        gnomADexomes_vcf_index = gnomADexomes_vcf_index,

        TopMed_vcf = TopMed_vcf,
        TopMed_vcf_index = TopMed_vcf_index,

        GnomAD4_exomes_vcf = GnomAD4_exomes_vcf,
        GnomAD4_exomes_vcf_index = GnomAD4_exomes_vcf_index,

        GnomAD4_genomes_vcf = GnomAD4_genomes_vcf,
        GnomAD4_genomes_vcf_index = GnomAD4_genomes_vcf_index,

        blacklisted_regions_bed = blacklisted_regions_bed,
        blacklisted_regions_bed_index = blacklisted_regions_bed_index,

        GERP_bed = GERP_bed,
        GERP_bed_index = GERP_bed_index,

        Regeneron_vcf = Regeneron_vcf,
        Regeneron_vcf_index = Regeneron_vcf_index,

        SLOpopulation_vcf = SLOpopulation_vcf,
        SLOpopulation_vcf_index = SLOpopulation_vcf_index,

        ClinVar_vcf = ClinVar_vcf,
        ClinVar_vcf_index = ClinVar_vcf_index,

        SpliceAI = SpliceAI,
        SpliceAI_index = SpliceAI_index,

        dbscSNV = dbscSNV,
        dbscSNV_index = dbscSNV_index,

        HPO = HPO,
        HPO_index = HPO_index,
        OMIM = OMIM,
        OMIM_index = OMIM_index,
        gnomadConstraints = gnomadConstraints,
        gnomadConstraints_index = gnomadConstraints_index,
        CGD = CGD,
        CGD_index = CGD_index,
        bcftools_annotation_header = bcftools_annotation_header,

        pext_bed = pext_bed,
        pext_bed_index = pext_bed_index,

        gnomAD_pLi_bed = gnomAD_pLi_bed,
        gnomAD_pLi_bed_index = gnomAD_pLi_bed_index,

        gnomAD_misz_bed = gnomAD_misz_bed,
        gnomAD_misz_bed_index = gnomAD_misz_bed_index,

        metaDome_bed = metaDome_bed,
        metaDome_bed_index = metaDome_bed_index,

        fasta_reference = fasta_reference,
        fasta_reference_index = fasta_reference_index,
        fasta_reference_dict = fasta_reference_dict,

        dbNSFP = dbNSFP,
        dbNSFP_index = dbNSFP_index

  }

  # ------------
  # OUTPUT FILES
  # ------------
  output {
    # alignments
    File out_bam        = upstream_hg19.out_bam
    File out_bam_index  = upstream_hg19.out_bam_index
    File out_cram       = ConvertToCram.output_cram
    File out_cram_index = ConvertToCram.output_cram_index

    # mosdepth outputs
    File? mosdepth_summary                 = upstream_hg19.mosdepth_summary
    #File mosdepth_region_bed              = upstream_hg19.mosdepth_region_bed
    #File mosdepth_region_bed_index        = upstream_hg19.mosdepth_region_bed_index
    #File mosdepth_depth_distribution_plot = upstream_hg19.mosdepth_depth_distribution_plot
    ##String inferred_sex                 = upstream_hg19.inferred_sex
    ##String stat_depth_mean              = upstream_hg19.stat_depth_mean
    #File output_mosdepth_summary                 = Rename_files.output_mosdepth_summary
    File? output_mosdepth_region_bed              = Rename_files.output_mosdepth_region_bed
    File? output_mosdepth_region_bed_index        = Rename_files.output_mosdepth_region_bed_index
    File? output_mosdepth_depth_distribution_plot = Rename_files.output_mosdepth_depth_distribution_plot

    # per sample sv signatures
    ### not needed: File discover_tar = upstream_hg19.discover_tar

    # sawfish outputs for single sample
    ##File? sv_vcf                        = upstream_hg19.sv_vcf
    ##File? sv_vcf_index                  = upstream_hg19.sv_vcf_index
    ##File? sv_supporting_reads           = upstream_hg19.sv_supporting_reads
    ##File? sv_copynum_bedgraph           = upstream_hg19.sv_copynum_bedgraph
    ##File? sv_depth_bw                   = upstream_hg19.sv_depth_bw
    ##File? sv_gc_bias_corrected_depth_bw = upstream_hg19.sv_gc_bias_corrected_depth_bw
    ##File? sv_maf_bw                     = upstream_hg19.sv_maf_bw
    ##File? sv_copynum_summary            = upstream_hg19.sv_copynum_summary
    File? output_sv_vcf = Rename_files.output_sv_vcf
    File? output_sv_vcf_index = Rename_files.output_sv_vcf_index
    File? output_sv_supporting_reads = Rename_files.output_sv_supporting_reads
    File? output_sv_copynum_bedgraph = Rename_files.output_sv_copynum_bedgraph
    File? output_sv_depth_bw = Rename_files.output_sv_depth_bw
    File? output_sv_gc_bias_corrected_depth_bw = Rename_files.output_sv_gc_bias_corrected_depth_bw
    File? output_sv_maf_bw = Rename_files.output_sv_maf_bw
    File? output_sv_copynum_summary = Rename_files.output_sv_copynum_summary
    
    File? sv_annotsv                    = annotSV.sv_variants_tsv

    # deep variant: small variant outputs
    ##File small_variant_vcf        = upstream_hg19.small_variant_vcf
    ##File small_variant_vcf_index  = upstream_hg19.small_variant_vcf_index
    ##File small_variant_gvcf       = upstream_hg19.small_variant_gvcf
    ##File small_variant_gvcf_index = upstream_hg19.small_variant_gvcf_index
    File? output_small_variant_vcf        = Rename_files.output_small_variant_vcf
    File? output_small_variant_vcf_index  = Rename_files.output_small_variant_vcf_index
    File? output_small_variant_gvcf       = Rename_files.output_small_variant_gvcf
    File? output_small_variant_gvcf_index = Rename_files.output_small_variant_gvcf_index

    File? deep_variant_vcf_modified           = TransformVcfFile.output_vcf
    File? deep_variant_vcf_modified_index     = TransformVcfFile.output_vcf_index
    File deep_variant_annotated_vcf           = AnnotateVCF_DeepVariant.output_vcf
    File deep_variant_annotated_vcf_index     = AnnotateVCF_DeepVariant.output_vcf_index
    File vep_deep_variant_annotated_vcf       = VEPDeepVariant.output_vcf
    File vep_deep_variant_annotated_vcf_index = VEPDeepVariant.output_vcf_index


    # trgt outputs
    File?   trgt_vcf                  = upstream_hg38.trgt_vcf
    File?   trgt_vcf_index            = upstream_hg38.trgt_vcf_index
    File?   trgt_spanning_reads       = upstream_hg38.trgt_spanning_reads
    File?   trgt_spanning_reads_index = upstream_hg38.trgt_spanning_reads_index
    #String stat_trgt_genotyped_count = upstream_hg38.stat_trgt_genotyped_count
    #String stat_trgt_uncalled_count  = upstream_hg38.stat_trgt_uncalled_count

    # paraphase outputs
    File? paraphase_output_json         = upstream_hg38.paraphase_output_json
    File? paraphase_realigned_bam       = upstream_hg38.paraphase_realigned_bam
    File? paraphase_realigned_bam_index = upstream_hg38.paraphase_realigned_bam_index
    File? paraphase_vcfs                = upstream_hg38.paraphase_vcfs

    # per sample mitorsaw outputs
    File? mitorsaw_vcf       = upstream_hg38.mitorsaw_vcf
    File? mitorsaw_vcf_index = upstream_hg38.mitorsaw_vcf_index
    File? mitorsaw_hap_stats = upstream_hg38.mitorsaw_hap_stats

  }
}

# ------
# TASKS 
# ------

# renaming some of the upstream output files so that they match our naming rules:
task Rename_files {
  input {
    String sample_id
    #File aligned_bam
    #File aligned_bam_index
    File small_variant_vcf
    File small_variant_vcf_index
    File small_variant_gvcf
    File small_variant_gvcf_index
    File? sv_vcf
    File? sv_vcf_index 
    File? sv_supporting_reads
    File? sv_copynum_bedgraph
    File? sv_depth_bw
    File? sv_gc_bias_corrected_depth_bw
    File? sv_maf_bw
    File? sv_copynum_summary
    File mosdepth_summary
    File mosdepth_region_bed
    File mosdepth_region_bed_index
    File mosdepth_depth_distribution_plot
  }

	#####cp  ~{aligned_bam} ~{sample_id}.aligned.bam
	#####cp  ~{aligned_bam_index} ~{sample_id}.aligned.bam.bai

  command {
	cp ~{small_variant_vcf} ~{sample_id}.DeepVariant.GRCh37.vcf.gz
	cp ~{small_variant_vcf_index} ~{sample_id}.DeepVariant.GRCh37.vcf.gz.tbi
	cp ~{small_variant_gvcf} ~{sample_id}.DeepVariant.gvcf.gz
	cp ~{small_variant_gvcf_index} ~{sample_id}.DeepVariant.gvcf.gz.tbi
	cp ~{sv_vcf} ~{sample_id}.hg19.sawfish.structural_variants.vcf.gz
	cp ~{sv_vcf_index} ~{sample_id}.hg19.sawfish.structural_variants.vcf.gz.tbi
	cp ~{sv_copynum_bedgraph} ~{sample_id}.hg19.sawfish.structural_variants.copynum.bedgraph
	cp ~{sv_copynum_summary} ~{sample_id}.hg19.sawfish.structural_variants.copynum.summary.json
	cp ~{sv_depth_bw} ~{sample_id}.hg19.sawfish.structural_variants.depth.bw
	cp ~{sv_gc_bias_corrected_depth_bw} ~{sample_id}.hg19.sawfish.structural_variants.gc_bias_corrected_depth.bw
	cp ~{sv_maf_bw} ~{sample_id}.hg19.sawfish.structural_variants.maf.bw
	cp ~{sv_supporting_reads} ~{sample_id}.hg19.sawfish.structural_variants.supporting_reads.json.gz
	cp ~{mosdepth_summary} ~{sample_id}.hg19.mosdepth.summary.txt.rename
	cp ~{mosdepth_region_bed} ~{sample_id}.hg19.mosdepth.regions.bed.gz
	cp ~{mosdepth_region_bed_index} ~{sample_id}.hg19.mosdepth.regions.bed.gz.csi
	cp ~{mosdepth_depth_distribution_plot} ~{sample_id}.hg19.mosdepth.depth_distribution.png
  }

  runtime {
    docker: "debian:buster-slim"
    cpu: 4
    requested_memory_mb_per_core: 1000
    runtime_minutes: 30
  }

  output {
    #File output_bam = "~{sample_id}.aligned.bam"
    #File output_bam_index = "~{sample_id}.aligned.bam.bai"
    File output_small_variant_vcf = "~{sample_id}.DeepVariant.GRCh37.vcf.gz"
    File output_small_variant_vcf_index = "~{sample_id}.DeepVariant.GRCh37.vcf.gz.tbi"
    File output_small_variant_gvcf = "~{sample_id}.DeepVariant.gvcf.gz"
    File output_small_variant_gvcf_index = "~{sample_id}.DeepVariant.gvcf.gz.tbi"
    File? output_sv_vcf = "~{sample_id}.hg19.sawfish.structural_variants.vcf.gz"
    File? output_sv_vcf_index = "~{sample_id}.hg19.sawfish.structural_variants.vcf.gz.tbi"
    File? output_sv_supporting_reads = "~{sample_id}.hg19.sawfish.structural_variants.supporting_reads.json.gz"
    File? output_sv_copynum_bedgraph = "~{sample_id}.hg19.sawfish.structural_variants.copynum.bedgraph"
    File? output_sv_depth_bw = "~{sample_id}.hg19.sawfish.structural_variants.depth.bw"
    File? output_sv_gc_bias_corrected_depth_bw = "~{sample_id}.hg19.sawfish.structural_variants.gc_bias_corrected_depth.bw"
    File? output_sv_maf_bw = "~{sample_id}.hg19.sawfish.structural_variants.maf.bw"
    File? output_sv_copynum_summary = "~{sample_id}.hg19.sawfish.structural_variants.copynum.summary.json"
	File  output_mosdepth_summary = "~{sample_id}.hg19.mosdepth.summary.txt.rename"
	File  output_mosdepth_region_bed = "~{sample_id}.hg19.mosdepth.regions.bed.gz"
	File  output_mosdepth_region_bed_index = "~{sample_id}.hg19.mosdepth.regions.bed.gz.csi"
	File  output_mosdepth_depth_distribution_plot = "~{sample_id}.hg19.mosdepth.depth_distribution.png"
  }
}

# prepare GRCh37 notated vcf file for use in hg19 notated environment of AnnotateVCF workflow; cut out all MT data
task TransformVcfFile {
  input {
    File input_vcf
	String output_vcf_name
  }

  command {
    zcat ~{input_vcf} | grep -v MT | sed -e 's/^/chr/' | sed -e 's/chr#/#/' | sed -e 's/contig=<ID=/contig=<ID=chr/' | bgzip -c > ~{output_vcf_name}
    tabix -p vcf ~{output_vcf_name}
  }

  output {
    File output_vcf = "~{output_vcf_name}"
    File output_vcf_index = "~{output_vcf_name}.tbi"
  }

  runtime {
    docker: "broadinstitute/gatk:4.2.0.0"
    cpu: 4
    requested_memory_mb_per_core: 1000
    runtime_minutes: 30
  }
}
