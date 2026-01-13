version 1.0
## Copyright CMG@KIGM, Ales Maver and Tanja Majnik

## Usage
# This workflow accepts three types of inputs: Illumina FASTQ files, a BAM file, or CRAM files
# Currently, the input CRAM files should be aligned to the hg19 reference genome assembly, we will implement support for other genome formats in the future
# The CRAM output is optional and disabled by default at the moment, until production switches to CRAM
# Manta is for genome analysis

# Subworkflows
import "./AnnotationPipeline.wdl" as Annotation
import "./CRAM_conversions.wdl" as CramConversions  
import "./Conifer.wdl" as Conifer
import "./Qualimap.wdl" as Qualimap
import "./ROH.wdl" as ROH
import "./CreateInterpretationTable.wdl" as CreateInterpretationTable
## import "./MitoMap.wdl" as MitoMap
import "./mitomap/mitomap_workflow.wdl" as MitoMap
import "./exp_hunter.wdl" as ExpansionHunter
import "./manta/manta_workflow.wdl" as Manta
import "./Delly/DELLY_single3.wdl" as Delly
import "./DeepVariant.wdl" as DeepVariant
import "./optimised_optitypeDNA" as Optitype
import "./SMN_caller/SMN_caller.wdl" as SMN
import "./SCRAMBLE/scramble_v2.wdl" as Scramble
import "./softsearch/softsearch_v2.wdl" as Softsearch
import "./bigWig/wigToBigWig_conversion" as BigWig
import "https://raw.githubusercontent.com/AlesMaver/gatk/master/scripts/mutect2_wdl/mutect2.wdl" as Mutect2
import "./MitochondriaPipeline/MitochondriaPipeline.wdl" as MitochondriaPipeline

# WORKFLOW DEFINITION 
workflow FastqToVCF {
  input {
    Array[File]? input_fq1
    Array[File]? input_fq2

    File? input_bam
    File? input_bam_index

    File? input_cram
    File? input_cram_index
    
    File? input_cram_hg38
    File? input_cram_hg38_index

    String sample_basename
    
    File illuminaAdapters

    File chromosome_list

    String? targetRegions
    Boolean? perform_masked_alignment
    Boolean call_somatic_variants = false

    Int bwa_threads
    Int threads
    Int split_reads_num

    File reference_fixed_fa
    File reference_fixed_fai
    File reference_fixed_amb
    File reference_fixed_ann
    File reference_fixed_bwt
    File reference_fixed_pac
    File reference_fixed_sa

    File reference_fa
    File reference_fai
    File reference_dict
        
    File reference_hg38_fa
    File reference_hg38_fai
    File reference_hg38_dict
    
    File gnomAD_vcf
    File gnomAD_vcf_index

    File gnomADexomes_vcf
    File gnomADexomes_vcf_index

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

    File dbNSFP
    File dbNSFP_index

    # Files containing frequency information for common SNPs used for ROH calculation
    File dbSNPcommon_bed
    File dbSNPcommon_bed_index
    File gnomAD_maf01_vcf
    File gnomAD_maf01_vcf_index
    File gnomAD_maf01_tab
    File gnomAD_maf01_tab_index

    File dbsnp_vcf
    File dbsnp_vcf_index
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices

    File refSeqFile
    
    Boolean make_bamout = false
    Boolean do_SoftSearch = true

    String? enrichment
    File? enrichment_bed

    Array[File]? input_reference_rpkms 
    Int? CONIFER_svd
    Float? CONIFER_threshold

    # Manta
    Array[File]? input_manta_reference_vcfs
    Boolean exome = false
    
    # DELLY
    File population_bcf
    File population_bcf_index

    ## Boolean GenerateCRAM = false
    Boolean? GenerateCRAM

    Boolean GVCFmode = false

    String? panel_gene_list
    Array[File]? relative_vcfs
    Array[File]? relative_vcf_indexes

    # TopMed
    File TopMed_vcf
    File TopMed_vcf_index

    File pext_bed
    File pext_bed_index

    File gnomAD_pLi_bed
    File gnomAD_pLi_bed_index

    File gnomAD_misz_bed
    File gnomAD_misz_bed_index

    File metaDome_bed
    File metaDome_bed_index

    File Regeneron_vcf
    File Regeneron_vcf_index

    File GnomAD4_exomes_vcf
    File GnomAD4_exomes_vcf_index

    # Mitochondria pipeline 
    Boolean do_MitochondriaPipeline = true
    File mt_dict
    File mt_fasta
    File mt_fasta_index
    File mt_amb
    File mt_ann
    File mt_bwt
    File mt_pac
    File mt_sa
    File blacklisted_sites
    File blacklisted_sites_index
    File mt_shifted_dict
    File mt_shifted_fasta
    File mt_shifted_fasta_index
    File mt_shifted_amb
    File mt_shifted_ann
    File mt_shifted_bwt
    File mt_shifted_pac
    File mt_shifted_sa
    File shift_back_chain
    File control_region_shifted_reference_interval_list
    File non_control_region_interval_list
    File gnomad_mito_sites_vcf
    File gnomad_mito_sites_vcf_index


    # Here are the global docker environment variables for tools used in this workflow
    # TO DO: Move the other task-specific docker definitions here for clarity, unless necessary
    String cutadapt_docker = "kfdrc/cutadapt:latest"
    # Uncommenting because the latest gatk docker might be causing problem with the HC genotyping step, forcing 4.2.0.0 temporarily
    # String gatk_docker = "broadinstitute/gatk:latest"
    String gatk_docker = "broadinstitute/gatk:4.2.0.0"
    String gatk_path = "/gatk/gatk"
    String gitc_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
    String samtools_path = "samtools" # Path to samtools command within GITC docker
    String vcfanno_docker = "clinicalgenomics/vcfanno:0.3.2"
    String bcftools_docker = "biocontainers/bcftools:v1.9-1-deb_cv1"
    String SnpEff_docker = "alesmaver/snpeff_v43:latest"
    String expansion_hunter_docker = "gbergant/expansionhunter:latest"
  }  

  # Terminate workflow in case neither input_fq1 or input_bam or input_cram is provided
  # Float fileSize = size(select_first([input_cram, input_bam, input_fq1, ""]))
  Float fileSize = size(select_first([input_cram_hg38, input_cram, input_bam, input_fq1, ""]))

  # Get sample name from either an input FASTQ R1 file or from the input BAM file - this causes problems with optional inputs, so it is left disabled and the input variable sample_basename is now a workflow input
  # String sample_basename = select_first([sub(basename(input_fq1), "[\_,\.].*", "" ), sub(basename(input_bam), "[\_,\.].*", "" )])

  # Get a list of chromosome contig names, containing one contig per line 
  Array[String] chromosomes = read_lines(chromosome_list)

  # Run adaptor trimming and create uBAM, but only when input FASTQ is provided
  if ( defined(input_fq1) ) {

    call MergeFastqFiles as MergeFastqFiles_fq1 {
      input:
        input_fq=input_fq1,
        sample_basename=sample_basename,
        read = "R1",
        docker = gatk_docker
    }

    call MergeFastqFiles as MergeFastqFiles_fq2 {
      input:
        input_fq=input_fq2,
        sample_basename=sample_basename,
        read = "R2",
        docker = gatk_docker
    }

    call CutAdapters as CutAdapters_fq1 {
      input:
        input_fq=MergeFastqFiles_fq1.output_fq_merged,
        sample_basename=sample_basename,
        illuminaAdapters=illuminaAdapters,

        # Runtime 
        docker = cutadapt_docker
    }

    call CutAdapters as CutAdapters_fq2 {
      input:
        input_fq=MergeFastqFiles_fq2.output_fq_merged,
        sample_basename=sample_basename,
        illuminaAdapters=illuminaAdapters,

        # Runtime 
        docker = cutadapt_docker
    }

    call PairedFastQsToUnmappedBAM {
        input:
          sample_name = sample_basename,
          fastq_1 = CutAdapters_fq1.output_fq_trimmed,
          fastq_2 = CutAdapters_fq2.output_fq_trimmed,
          readgroup_name = sample_basename,
          library_name = sample_basename,
          platform_unit = "NextSeq550",
          run_date = "2020-01-01",
          platform_name = "illumina",
          sequencing_center = "KIGM",
          gatk_path = gatk_path,
          docker = gatk_docker
      }
  }

  if ( defined(input_cram) ) {
    call CramConversions.CramToBam as Cram_hg19_ToBam {
      input:
        input_cram = input_cram,
        sample_name = sample_basename,
        ref_dict = reference_dict,
        ref_fasta = reference_fa,
        ref_fasta_index = reference_fai,
        docker = gitc_docker,
        samtools_path = samtools_path
    }
  }

  if ( defined(input_cram_hg38) ) {
    call CramConversions.CramToBam as Cram_hg38_ToBam {
      input:
        input_cram = input_cram_hg38,
        sample_name = sample_basename,
        ref_dict = reference_hg38_dict,
        ref_fasta = reference_hg38_fa,
        ref_fasta_index = reference_hg38_fai,
        docker = gitc_docker,
        samtools_path = samtools_path
    }
  }
  
  call SamSplitter {
    input :
      input_bam = select_first([Cram_hg38_ToBam.output_bam, Cram_hg19_ToBam.output_bam, input_bam, PairedFastQsToUnmappedBAM.output_unmapped_bam]),
      n_reads = split_reads_num,
      preemptible_tries = 3,
      compression_level = 2
  }

  if ( defined(targetRegions) && select_first([perform_masked_alignment, false]) ) {
    call PrepareMaskedGenomeFasta {
      input:
        reference_fixed_fa=reference_fixed_fa,
        reference_fixed_fai=reference_fixed_fai,
        targetRegions=targetRegions,

        # Runtime 
        docker = "pegi3s/bedtools"
    }

    call PrepareMaskedBWAIndex {
      input:
        reference_masked_fa=PrepareMaskedGenomeFasta.reference_masked_fa,

        # Runtime 
        docker="alesmaver/bwa_samtools_picard"
    }
  }

  scatter(unmapped_bam in SamSplitter.split_bams) {
    call Align {
      input:
       input_bam=unmapped_bam,
       sample_basename=sample_basename,

       reference_fixed_fa=select_first([PrepareMaskedBWAIndex.reference_masked_fasta, reference_fixed_fa]),
       reference_fixed_fai=reference_fixed_fai,
       reference_fixed_amb=select_first([PrepareMaskedBWAIndex.reference_masked_amb, reference_fixed_amb]),
       reference_fixed_ann=select_first([PrepareMaskedBWAIndex.reference_masked_ann, reference_fixed_ann]),
       reference_fixed_bwt=select_first([PrepareMaskedBWAIndex.reference_masked_bwt, reference_fixed_bwt]),
       reference_fixed_pac=select_first([PrepareMaskedBWAIndex.reference_masked_pac, reference_fixed_pac]),
       reference_fixed_sa=select_first([PrepareMaskedBWAIndex.reference_masked_sa, reference_fixed_sa]),

       reference_fa=reference_fa,
       reference_dict=reference_dict,

       picard_path="/usr/picard/picard.jar",

       threads=bwa_threads,

       docker="alesmaver/bwa_samtools_picard"
    }
  }

  call GatherUnsortedBamFiles {
    input:
      input_bams = Align.output_bam,
      output_bam_basename = sample_basename,
      preemptible_tries = 3,
      compression_level = 2
  }

  # Sort aggregated+deduped BAM file and fix tags
  call SortSam {
    input:
      input_bam = GatherUnsortedBamFiles.output_bam,
      output_bam_basename = sample_basename + ".marked",
      compression_level = 2,
      preemptible_tries = 3
  }

  Float agg_bam_size = size(SortSam.output_bam, "GiB")

  # Create list of sequences for scatter-gather parallelization
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = reference_dict,
      preemptible_tries = 3
  }

  # We need disk to localize the sharded input and output due to the scatter for BQSR.
  # If we take the number we are scattering by and reduce by 3 we will have enough disk space
  # to account for the fact that the data is not split evenly.
  Int num_of_bqsr_scatters = length(CreateSequenceGroupingTSV.sequence_grouping)
  Int potential_bqsr_divisor = num_of_bqsr_scatters - 10
  Int bqsr_divisor = if potential_bqsr_divisor > 1 then potential_bqsr_divisor else 1

  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        input_bam = SortSam.output_bam,
        input_bam_index = SortSam.output_bam_index,
        recalibration_report_filename = sample_basename + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        known_indels_sites_vcfs = known_indels_sites_vcfs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = reference_dict,
        ref_fasta = reference_fa,
        ref_fasta_index = reference_fai,
        bqsr_scatter = bqsr_divisor,
        preemptible_tries = 3
    }
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  # The reports are always the same size
  call GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = sample_basename + ".recal_data.csv",
      preemptible_tries = 3
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        input_bam = SortSam.output_bam,
        input_bam_index = SortSam.output_bam_index,
        output_bam_basename = sample_basename,
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = reference_dict,
        ref_fasta = reference_fa,
        ref_fasta_index = reference_fai,
        bqsr_scatter = bqsr_divisor,
        compression_level = 2,
        preemptible_tries = 3,
        bin_base_qualities = true,
        somatic = false
    }
  }

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherSortedBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = sample_basename,
      total_input_size = agg_bam_size,
      compression_level = 2,
      preemptible_tries = 3
  }

  if ( select_first([GenerateCRAM, false])) {
	  call CramConversions.ConvertToCram as ConvertToCram {
	    input:
	      input_bam = SortSam.output_bam,
	      ref_fasta = reference_fa,
	      ref_fasta_index = reference_fai,
	      sample_basename = sample_basename
  	}
  }
  
  
  # Calculate Optitype only if targetRegions are not present
  if ( !defined(targetRegions) ) {
    call Optitype.Optimised_OptitypeDna as Optitype {
      input:
        sample_basename = sample_basename,
        input_fq1 = CutAdapters_fq1.output_fq_trimmed,
        input_fq2 = CutAdapters_fq2.output_fq_trimmed,
	### input_bam = select_first([Cram_hg38_ToBam.output_bam, Cram_hg19_ToBam.output_bam, input_bam])
        input_bam = input_bam
    }
  }
  

  # Calculate SMN only if targetRegions are not present
  if( !defined(targetRegions) ) {
     call SMN.SMN_caller as SMN_caller {
     input:
       input_bam=SortSam.output_bam,
       input_bam_index=SortSam.output_bam_index,
       sample_basename=sample_basename
     }
  }
  
  
  # Calculate Scramble only if targetRegions are not present
  if( !defined(targetRegions) ) {
     call Scramble.SCRAMBLE_workflow as Scramble {
     input:
       input_bam = SortSam.output_bam,
       input_bam_index = SortSam.output_bam_index,
       reference_fa = reference_fa,
       reference_fai = reference_fai,
       reference_dict = reference_dict 
     }
  }
  
  
  # Calculate Softsearch only if targetRegions are not present
  if (do_SoftSearch) {
    if( !defined(targetRegions) ) {
      call Softsearch.SoftSearchWF as SoftsearchWF {
        input:
          input_bam = SortSam.output_bam,
          input_bam_index = SortSam.output_bam_index,
          reference_fa = reference_fa,
          reference_fai = reference_fai,
          reference_dict = reference_dict,
          sample_basename = sample_basename
      }
    }
  }

  # Mitochondria pipeline 
  if (do_MitochondriaPipeline) {
      call MitochondriaPipeline.MitochondriaPipeline as MitochondriaPipeline {
        input:
          sample_basename = sample_basename,
          wgs_aligned_input_bam_or_cram = SortSam.output_bam,
          wgs_aligned_input_bam_or_cram_index = SortSam.output_bam_index,
          ref_fasta = reference_fa,
          ref_fasta_index = reference_fai,
          ref_dict = reference_dict,
          mt_fasta = mt_fasta,
          mt_fasta_index = mt_fasta_index,
          mt_dict = mt_dict,
          mt_amb = mt_amb,
          mt_ann = mt_ann,
          mt_bwt = mt_bwt, 
          mt_pac = mt_pac, 
          mt_sa = mt_sa, 
          blacklisted_sites = blacklisted_sites,
          blacklisted_sites_index = blacklisted_sites_index, 
          mt_shifted_dict = mt_shifted_dict, 
          mt_shifted_fasta = mt_shifted_fasta, 
          mt_shifted_fasta_index = mt_shifted_fasta_index, 
          mt_shifted_amb = mt_shifted_amb, 
          mt_shifted_ann = mt_shifted_ann, 
          mt_shifted_bwt = mt_shifted_bwt, 
          mt_shifted_pac = mt_shifted_pac, 
          mt_shifted_sa = mt_shifted_sa, 
          shift_back_chain = shift_back_chain, 
          control_region_shifted_reference_interval_list = control_region_shifted_reference_interval_list, 
          non_control_region_interval_list = non_control_region_interval_list, 
          gnomad_mito_sites_vcf = gnomad_mito_sites_vcf, 
          gnomad_mito_sites_vcf_index = gnomad_mito_sites_vcf_index 
      }
  }


  if( defined(targetRegions) ) {
    call StringToArray {
      input:
        input_string = select_first([targetRegions, ""]),
        separator = ";"
    }
  }

  if( defined(targetRegions) && call_somatic_variants ) {
    call Mutect2.Mutect2 {
        input:
          ref_fasta = reference_fa,
          ref_fai = reference_fai,
          ref_dict = reference_dict, 
          scatter_count = 1,
          intervals = StringToArray.intervals_list,
          tumor_reads = GatherSortedBamFiles.output_bam,
          tumor_reads_index = GatherSortedBamFiles.output_bam_index,
          make_bamout = true, 
          # We mostly use Mutect2 for detection using PCR-based approaches, therefore disabling NotDuplicateReadFilter for testing purposes
          m2_extra_args = " --disable-read-filter NotDuplicateReadFilter ",
          gatk_docker = gatk_docker
    }
  }

  scatter (chromosome in select_first([StringToArray.values, chromosomes]) ) {
    call HaplotypeCaller {
      input:
        input_bam = GatherSortedBamFiles.output_bam,
        input_bai = GatherSortedBamFiles.output_bam_index,
        sample_basename=sample_basename,

        reference_fa=reference_fa,
        reference_fai=reference_fai,
        reference_dict=reference_dict,

        chromosome=chromosome,
	
	make_bamout = make_bamout,

        threads=threads,

        gatk_path=gatk_path,
        docker=gatk_docker
    }
  }

  # Merge per-interval GVCFs
  call MergeVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_vcf,
      input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
      sample_basename = sample_basename,
      docker = gatk_docker,
      gatk_path = gatk_path
  }
  
  if (make_bamout) {
    call MergeBamOuts {
      input:
	ref_fasta = reference_fa,
	ref_fai = reference_fai,
	ref_dict = reference_dict,
	bam_outs = HaplotypeCaller.output_bamout,
	sample_basename = sample_basename,
        docker = gatk_docker,
        gatk_path = gatk_path
    }
  }

  call SplitSNPindel {
  input:
    input_vcf = MergeVCFs.output_vcf,
    input_vcf_index = MergeVCFs.output_vcf_index,
    sample_basename=sample_basename,

    reference_fa=reference_fa,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    threads=threads,

    gatk_path=gatk_path,
    docker=gatk_docker
  }

  call VariantFiltrationSNP {
  input:
    input_vcf = SplitSNPindel.output_vcf_SNP,
    sample_basename=sample_basename,

    reference_fa=reference_fa,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    threads=threads,

    gatk_path=gatk_path,
    docker=gatk_docker
  }

  call VariantFiltrationINDEL {
  input:
    input_vcf = SplitSNPindel.output_vcf_INDEL,
    sample_basename=sample_basename,

    reference_fa=reference_fa,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    threads=threads,

    gatk_path=gatk_path,
    docker=gatk_docker
  }

  call CombineSNPindel {
  input:
    input_snp_vcf = VariantFiltrationSNP.output_vcf,
    input_indel_vcf = VariantFiltrationINDEL.output_vcf,
    sample_basename=sample_basename,

    reference_fa=reference_fa,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    threads=threads,

    picard_path="/usr/picard/picard.jar",
    docker="alesmaver/bwa_samtools_picard"
  }

  call SelectFinalVariants {
  input:
    input_vcf = CombineSNPindel.output_vcf,
    sample_basename=sample_basename,

    reference_fa=reference_fa,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    threads=threads,

    gatk_path=gatk_path,
    docker=gatk_docker
  }

  call Annotation.AnnotateVCF as AnnotateVCF{
    input:
      input_vcf = select_first([Mutect2.filtered_vcf, SelectFinalVariants.output_vcf]),
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

      fasta_reference = reference_fa,
      fasta_reference_index = reference_fai,
      fasta_reference_dict = reference_dict,

      dbNSFP = dbNSFP,
      dbNSFP_index = dbNSFP_index,

      targetRegions = targetRegions,

      #bcftools_docker = bcftools_docker,
      #SnpEff_docker = SnpEff_docker,
      gatk_docker = gatk_docker,
      gatk_path = gatk_path,
      vcfanno_docker = vcfanno_docker
  }

  #call MitoMap.CreateMitoFasta as CreateMitoFasta {
  #  input:
  #  input_vcf = SelectFinalVariants.output_vcf,
  #  sample_basename = sample_basename,

  #  reference_fa = reference_fa,
  #  reference_fai = reference_fai,
  #  reference_dict = reference_dict,

  #  docker = "broadinstitute/gatk3:3.8-1"
  #}

  #call MitoMap.MitoMap as MitoMap {
  #  input:
  #  mtDNA_fasta = CreateMitoFasta.mtDNA_fasta,
  #  sample_basename = sample_basename
  #}
  
  call MitoMap.MitoMapWorkflow as MitoMapWorkflow  {
    input:
      input_vcf = SelectFinalVariants.output_vcf,
      sample_basename = sample_basename,
      reference_fa = reference_fa,
      reference_fai = reference_fai,
      reference_dict = reference_dict,
  }

  call CreateInterpretationTable.CreateInterpretationTable as CreateInterpretationTable {
    input:
      input_vcf = AnnotateVCF.output_vcf,
      input_vcf_index = AnnotateVCF.output_vcf_index,
      relative_vcfs = relative_vcfs,
      relative_vcf_indexes = relative_vcf_indexes,
      panel_gene_list = panel_gene_list,
      mitoResults_txt = MitoMapWorkflow.mitoResults_txt
  }

  if( defined(input_reference_rpkms) ){
    call Conifer.Conifer as Conifer{
    input:
      input_bam = SortSam.output_bam,
      input_bam_index = SortSam.output_bam_index,

      input_reference_rpkms = input_reference_rpkms,
      CONIFER_svd = CONIFER_svd,
      CONIFER_threshold = CONIFER_threshold,

      enrichment = enrichment,
      enrichment_bed = enrichment_bed
    }

    call Conifer.CONIFER_Annotate as CONIFER_Annotate{
    input:
      conifer_calls_wig = Conifer.output_conifer_calls_wig,
      
      sample_basename = sample_basename,

      HPO = HPO,
      HPO_index = HPO_index,
      OMIM = OMIM,
      OMIM_index = OMIM_index,
      gnomadConstraints = gnomadConstraints,
      gnomadConstraints_index = gnomadConstraints_index,
      CGD = CGD,
      CGD_index = CGD_index
    }   

    call Conifer.GenerateAnnotatedIntervalFile as CONIFER_GenerateAnnotatedIntervalFile {
    input:
      bedtools_annotated_file = CONIFER_Annotate.conifer_bedtools
    }     
  }

  if( defined(input_manta_reference_vcfs) && !defined(targetRegions) ){
    call Manta.SVcalling as Manta{
    input:
      bamFile = SortSam.output_bam,
      bamIndex = SortSam.output_bam_index,

      referenceFasta=reference_fa,
      referenceFastaFai=reference_fai,
      referenceFastaDict=reference_dict,  

      exome = false,
      sample = sample_basename, 

      input_manta_reference_vcfs = input_manta_reference_vcfs, 
    }

    call Annotation.bcftoolsAnnotate as MantaAnnotation {
      input:
        input_vcf = Manta.output_manta_filtered_vcf,
        input_vcf_index = Manta.mantaVcfindex,

        sample_basename=sample_basename,

        HPO = HPO,
        HPO_index = HPO_index,
        OMIM = OMIM,
        OMIM_index = OMIM_index,
        gnomadConstraints = gnomadConstraints,
        gnomadConstraints_index = gnomadConstraints_index,
        CGD = CGD,
        CGD_index = CGD_index,

        bcftools_annotation_header = bcftools_annotation_header,
        
        output_filename = sample_basename + ".manta.annotated.vcf.gz",
        docker = "dceoy/bcftools"
    }
  }

  # DELLY
  if ( defined(enrichment) ){
    if( enrichment=="WGS1Mb" ){
       call Delly.DELLY as Delly {
       input:
   	 input_bam = SortSam.output_bam,
         input_bai = SortSam.output_bam_index,
         population_bcf = population_bcf,
         population_bcf_index = population_bcf_index,
         reference_fasta = reference_fa
       }
    }
  }

  # Deep Variant 
  if (defined(enrichment)) {
    String model_type = if (enrichment == "WGS1Mb") then "WGS" else "WES"
    call DeepVariant.DeepVariant as DeepVariant {
       input:
	  sample_basename = sample_basename,
          inputBam = SortSam.output_bam,
          inputBamIndex = SortSam.output_bam_index,
          referenceFasta = reference_fa,
          referenceFastaIndex = reference_fai,
	  modelType = model_type,
          # numShards = 128
    }
  }

  if ( defined(targetRegions) ){
    call RegionsToBed {
      input:
        targetRegions=select_first([targetRegions,""]),

        # Runtime 
        docker = "pegi3s/bedtools"

    }
  }

  if( defined(enrichment_bed) || defined(PrepareMaskedGenomeFasta.targetRegions_bed) || defined(RegionsToBed.targetRegions_bed) ){
    call Qualimap.bamqc as Qualimap {
    input:
      bam = SortSam.output_bam,
      sample_basename=sample_basename,
      
      enrichment_bed = select_first([PrepareMaskedGenomeFasta.targetRegions_bed, RegionsToBed.targetRegions_bed, enrichment_bed]),

      ncpu = 8
    }
  }

  if ( defined(enrichment) ){
    if( enrichment=="WGS1Mb" ){
       call Qualimap.DownsampleBED as DownsampleBED {
        input:
           bed_file = enrichment_bed,
           reference_fai=reference_fai
      }
    }
  }

  # Depth of coverage
  if( defined(enrichment_bed) || defined(PrepareMaskedGenomeFasta.targetRegions_bed) || defined(RegionsToBed.targetRegions_bed) ){
    call Qualimap.DepthOfCoverage34 as DepthOfCoverage {
      input:
        input_bam = SortSam.output_bam,
        input_bam_index = SortSam.output_bam_index,
        sample_basename = sample_basename,

        reference_fa=reference_fa,
        reference_fai=reference_fai,
        reference_dict=reference_dict,

        enrichment_bed = select_first([PrepareMaskedGenomeFasta.targetRegions_bed, RegionsToBed.targetRegions_bed, DownsampleBED.downsampled_bed_file, enrichment_bed]),

        refSeqFile = refSeqFile,

        threads = threads,
        docker = "broadinstitute/gatk3:3.8-1",
        gatk_path = "/usr/GenomeAnalysisTK.jar"
    }
  }
  
  # Calculate WGS coverage if neither enrichment_bed nor target regions parameter is defined
  if( !defined(enrichment_bed) && !defined(PrepareMaskedGenomeFasta.targetRegions_bed) && !defined(RegionsToBed.targetRegions_bed) ){
  
    call Qualimap.bamqc as QualimapWGS {
      input:
        bam = SortSam.output_bam,
        sample_basename=sample_basename,

        ncpu = 8
    }
  
    call Qualimap.DepthOfCoverage34 as DepthOfCoverageWGS {
      input:
        input_bam = SortSam.output_bam,
        input_bam_index = SortSam.output_bam_index,
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

  # Downsample dbSNP bed file if this is a WGS analysis
  if ( defined(enrichment) ){
    if( enrichment=="WGS1Mb" ){
       call ROH.Downsample_dbSNP as Downsample_dbSNP {
        input:
          dbSNPcommon_bed = dbSNPcommon_bed,
          dbSNPcommon_bed_index = dbSNPcommon_bed_index
      }
    }
  }

  # Do not perform ROH calling if target regions are defined - this means targeted sequencing
  # Currently ROHs are left on because they do not disrupt the WF and they generate the files, even if empty, and make it easier for copying, etc
  #if ( !defined(targetRegions) ){
    call ROH.calculateBAF as calculateBAF {
    input:
      input_bam = SortSam.output_bam,
      input_bam_index = SortSam.output_bam_index,
      sample_basename=sample_basename,

      reference_fa=reference_fa,

      dbSNPcommon_bed = select_first([Downsample_dbSNP.downsampled_dbSNPcommon_bed, dbSNPcommon_bed]),
      dbSNPcommon_bed_index = select_first([Downsample_dbSNP.downsampled_dbSNPcommon_bed_index, dbSNPcommon_bed_index]),

      docker = "alesmaver/bwa_samtools_picard"
    }

    call ROH.CallROH as CallROH {
    input:
      input_bam = SortSam.output_bam,
      input_bam_index = SortSam.output_bam_index,
      sample_basename=sample_basename,
    
      reference_fa=reference_fa,
    
      dbSNPcommon_bed = select_first([Downsample_dbSNP.downsampled_dbSNPcommon_bed, dbSNPcommon_bed]),
      dbSNPcommon_bed_index = select_first([Downsample_dbSNP.downsampled_dbSNPcommon_bed_index, dbSNPcommon_bed_index]),
    
      gnomAD_maf01_vcf = gnomAD_maf01_vcf,
      gnomAD_maf01_vcf_index = gnomAD_maf01_vcf_index,

      gnomAD_maf01_tab = gnomAD_maf01_tab,
      gnomAD_maf01_tab_index = gnomAD_maf01_tab_index,
    
      docker = bcftools_docker
    }

    call Manta.annotSV as ROH_annotSV {
        input:
          genome_build = "GRCh37",
          input_vcf = select_first([CallROH.ROH_calls_annotSV_input_bed, ""]),
          output_tsv_name = sample_basename + ".ROH.annotSV.tsv"
    }
  #}

  call ExpansionHunter.ExpansionHunter as ExpansionHunter {
    input:
      sample_id = sample_basename,
      bam_file = SortSam.output_bam,
      bai_file = SortSam.output_bam_index,
      reference_fasta = reference_fa,
      expansion_hunter_docker = expansion_hunter_docker
  }

  #call ROH.CallPlink as CallPlink {
  #input:
  #  input_vcf = CallROH.BAF_vcf,
  #  sample_basename=sample_basename,
  #
  #  docker = "asherkhb/plink"
  #}
  
  call BigWig.wigToBigWigConversion as BigWig {
  input:
    sample_basename = sample_basename,
    ref_fasta_index = reference_fai,
    baf_wig = calculateBAF.output_BAF,
    roh_calls_qual_wig = CallROH.ROH_calls_qual,
    cnv_genome_wig = Conifer.CNV_wig,
    cnv_wig = Conifer.output_conifer_calls_wig,
    coverage_neg_wig = DepthOfCoverage.coverage_neg_wig,
    coverage_wig = DepthOfCoverage.coverage_wig,
    coverage_mean_wig = DepthOfCoverage.coverage_mean_wig
  }

  output {
    File output_bam = SortSam.output_bam
    File output_bam_index = SortSam.output_bam_index

    File? output_cram = ConvertToCram.output_cram
    File? output_cram_index = ConvertToCram.output_cram_index

    File output_vcf_raw = MergeVCFs.output_vcf
    File output_vcf_raw_index = MergeVCFs.output_vcf_index

    File output_vcf = SelectFinalVariants.output_vcf
    File output_vcf_index = SelectFinalVariants.output_vcf_index

    File output_annotated_vcf = AnnotateVCF.output_vcf
    File output_annotated_vcf_index = AnnotateVCF.output_vcf_index
    File? XLSX_OUTPUT = CreateInterpretationTable.XLSX_OUTPUT

    File? output_rpkm = Conifer.output_rpkm 
    File? output_conifer_calls = Conifer.output_conifer_calls
    Array[File]? output_plotcalls = Conifer.output_plotcalls
    File? output_conifer_calls_wig = Conifer.output_conifer_calls_wig
    #File CNV_bed = Conifer.CNV_bed
    File? CNV_wig = Conifer.CNV_wig
    File? conifer_calls_annotated = CONIFER_Annotate.conifer_calls_annotated
    File? conifer_annotSV_tsv = Conifer.annotSV_tsv
    File? conifer_plots_tar = Conifer.conifer_plots_tar

    File? mantaVCF = Manta.mantaVcf
    File? mantaVCFindex = Manta.mantaVcfindex
    File? mantaSVs = Manta.output_sv_table
    File? mantaSVs_annotSV_tsv = Manta.annotSV_tsv
    
    # Delly
    File? call_bcf_file = Delly.call_bcf_file
    File? sample_bcfs = Delly.sample_bcfs
    File? sample_bcf_indices = Delly.sample_bcf_indices
    #File? unfiltered_population_bcf = Delly.sample_bcf_indices
    File? filtered_population_bcf = Delly.filtered_population_bcf
    File? delly_annotSV = Delly.delly_annotSV 

    # Deep Variant
    File? outputVCF = DeepVariant.outputVCF
    File? outputVCFIndex = DeepVariant.outputVCFIndex
    File? outputVCFStatsReport = DeepVariant.outputVCFStatsReport
    File? outputGVCF = DeepVariant.outputGVCF
    File? outputGVCFIndex = DeepVariant.outputGVCFIndex

    File? Qualimap_results = Qualimap.results
    File? QualimapWGS_results = QualimapWGS.results

    File? DepthOfCoverage_output = DepthOfCoverage.DepthOfCoverage_output
    File? DepthOfCoverageWGS_output = DepthOfCoverageWGS.DepthOfCoverage_output

    File? output_BAF = calculateBAF.output_BAF
    File? ROH_calls_qual = CallROH.ROH_calls_qual
    File? ROH_calls_size = CallROH.ROH_calls_size
    File? ROH_intervals_state = CallROH.ROH_intervals_state
    File? ROH_intervals_qual = CallROH.ROH_intervals_qual
    File? ROH_annotSV_tsv = ROH_annotSV.sv_variants_tsv
    #File ROHplink_calls = CallPlink.ROHplink_calls

    ##File? mitoResults_xls = MitoMap.mitoResults_xls
    ##File? mitoResults_txt = MitoMap.mitoResults_txt
    File? mitoResults_xls = MitoMapWorkflow.mitoResults_xls
    File? mitoResults_txt = MitoMapWorkflow.mitoResults_txt
    
    File? optitype_tsv = Optitype.optitype_tsv
    File? optitype_plot = Optitype.optitype_plot

    File? output_tsv = SMN_caller.output_tsv
    File? output_json = SMN_caller.output_json
    File? output_pdf = SMN_caller.output_pdf

    File? MEI_output = Scramble.MEI_output
    File? MEI_output_annotated = Scramble.MEI_output_annotated
    
    File? softsearch_complete_vcf = SoftsearchWF.output_complete_vcf
    File? softsearch_complete_vcf_index = SoftsearchWF.output_complete_vcf_index
    File? softsearch_filter_vcf = SoftsearchWF.output_vcf
    File? softsearch_annotSV = SoftsearchWF.output_tsv_name

    File? expansion_hunter_vcf_annotated = ExpansionHunter.expansion_hunter_vcf_annotated
    
    # merged haplotype caller bamout:
    File? bamout = MergeBamOuts.merged_bam_out
    File? bamout_index = MergeBamOuts.merged_bam_out_index
    
    ## Array[File] bigWig_files = BigWig.bigWig_files
    File? bigWig_tar_file = BigWig.bigWig_tar_file

    # Mitochondria pipeline
    File? mt_aligned_bam = MitochondriaPipeline.mt_aligned_bam
    File? mt_aligned_bai = MitochondriaPipeline.mt_aligned_bai
    File? split_vcf = MitochondriaPipeline.split_vcf
    File? split_vcf_index = MitochondriaPipeline.split_vcf_index
    File? input_vcf_for_haplochecker = MitochondriaPipeline.input_vcf_for_haplochecker
    File? duplicate_metrics = MitochondriaPipeline.duplicate_metrics
    File? coverage_metrics = MitochondriaPipeline.coverage_metrics
    File? coverage_mean_metrics = MitochondriaPipeline.coverage_mean_metrics
    File? coverage_median_metrics = MitochondriaPipeline.coverage_median_metrics
    File? major_haplogroup_file = MitochondriaPipeline.major_haplogroup_file
    File? base_level_coverage_metrics = MitochondriaPipeline.base_level_coverage_metrics
    File? filterVCF_output_vcf = MitochondriaPipeline.filterVCF_output_vcf
    File? filterVCF_mito_table = MitochondriaPipeline.filterVCF_mito_table

  }
}


##################
# TASK DEFINITIONS
##################

task MergeFastqFiles {
  input {
    Array[File]? input_fq
    String sample_basename
    String read

    # Runtime parameters
    String docker
  }

  command {
  set -e
     cat ${sep = ' ' input_fq} > ~{sample_basename}_~{read}.fq.gz
  }
  runtime {
    docker: docker
    requested_memory_mb_per_core: 1000
    cpu: 4
    runtime_minutes: 800
  }
  output {
    File output_fq_merged = "~{sample_basename}_~{read}.fq.gz"
  }
}

task CutAdapters {
  input {
    # Command parameters
    File? input_fq # The fastq file is given the optional flag to correspond to optional fastq input in the workflow
    String sample_basename
    
    File illuminaAdapters

    # Runtime parameters
    String docker
  }
  
  command {
  set -e
     cutadapt -j 4 -a file:~{illuminaAdapters} --mask-adapter -o ~{sample_basename}.trimmed.fq.gz ~{input_fq}
  }
  runtime {
    docker: docker
    requested_memory_mb_per_core: 1000
    cpu: 4
    runtime_minutes: 800
  }
  output {
    File output_fq_trimmed = "~{sample_basename}.trimmed.fq.gz"
  }
}

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

task PrepareMaskedGenomeFasta {
  input {
    # Command parameters
    File reference_fixed_fa
    File reference_fixed_fai

    String? targetRegions # FORMAT "chr1:123033-130000;chrX:1-1000"

    # Runtime parameters
    String docker
  }
  
  command <<<
    set -e

    echo "~{targetRegions}"  | tr ';' '\n' | tr ':' '\t' | tr '-' '\t' > targetRegions.bed

    awk -v OFS='\t' {'print $1,$2'} ~{reference_fixed_fai} > ~{reference_fixed_fa}.genome

    bedtools complement  -i targetRegions.bed -g ~{reference_fixed_fa}.genome > targetRegions.complement.bed

    bedtools maskfasta -fi ~{reference_fixed_fa} -bed targetRegions.complement.bed -fo reference.masked.fa
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 4000
    cpu: 4
    runtime_minutes: 800
  }
  output {
    File reference_masked_fa = "reference.masked.fa"
    File targetRegions_bed = "targetRegions.bed"
  }
}

task PrepareMaskedBWAIndex {
  input {
    # Command parameters
    File reference_masked_fa
    String reference_masked_fa_filename = basename(reference_masked_fa)
    # Runtime parameters
    String docker
  }
  
  command {
    set -e

    cp ~{reference_masked_fa} ./

    bwa index ~{reference_masked_fa_filename}
  }
  runtime {
    docker: docker
    requested_memory_mb_per_core: 5000
    cpu: 4
    runtime_minutes: 800
  }
  output {
    File reference_masked_fasta = "~{reference_masked_fa_filename}"
    File reference_masked_sa = "~{reference_masked_fa_filename}.sa"
    File reference_masked_amb = "~{reference_masked_fa_filename}.amb"
    File reference_masked_ann = "~{reference_masked_fa_filename}.ann"
    File reference_masked_bwt = "~{reference_masked_fa_filename}.bwt"
    File reference_masked_pac = "~{reference_masked_fa_filename}.pac"
  }
}

task Align {
  input {
  File input_bam

  File reference_fixed_fa
  File reference_fixed_fai
  File reference_fixed_amb
  File reference_fixed_ann
  File reference_fixed_bwt
  File reference_fixed_pac
  File reference_fixed_sa

  File reference_fa
  File reference_dict

  String sample_basename

  String picard_path="/usr/picard/picard.jar"

  # Runtime params
  Int threads
  String docker
  }
  
  command {
    java -Xms1000m -Xmx1000m -jar ~{picard_path} SamToFastq INPUT=~{input_bam} FASTQ=/dev/stdout INTERLEAVE=true | bwa mem -p -t ~{threads} ~{reference_fixed_fa} /dev/stdin | samtools sort -@ ~{threads} - | java -jar ~{picard_path} AddOrReplaceReadGroups I=/dev/stdin O=/dev/stdout RGID=4 RGLB=~{sample_basename} RGPL=illumina RGPU=unit1 RGSM=~{sample_basename} | java -jar ~{picard_path} ReorderSam I=/dev/stdin O=~{sample_basename}.sorted.bam REFERENCE_SEQUENCE=~{reference_fa} SEQUENCE_DICTIONARY=~{reference_dict}

    java -jar ~{picard_path} MarkDuplicates  I=~{sample_basename}.sorted.bam O=~{sample_basename}.marked.bam M=~{sample_basename}.metrics.txt CREATE_INDEX=true
  }

  output {
    File output_bam = "~{sample_basename}.marked.bam"
    File output_bam_index = "~{sample_basename}.marked.bai"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
    requested_memory_mb_per_core: 2000
    cpu: 8
    runtime_minutes: 300
  }
}

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  input {
    File ref_dict
    Int preemptible_tries
  }
  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("~{ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
    memory: "2 GiB"
    requested_memory_mb_per_core: 1000
    cpu: 1
    runtime_minutes: 20
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}


task BaseRecalibrator {
  input {
    File input_bam
    File input_bam_index
    String recalibration_report_filename
    Array[String] sequence_group_interval
    File dbsnp_vcf
    File dbsnp_vcf_index
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int bqsr_scatter
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.10.1"
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Float dbsnp_size = size(dbsnp_vcf, "GiB")
  Int disk_size = ceil((size(input_bam, "GiB") / bqsr_scatter) + ref_size + dbsnp_size) + 20

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command {
    gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Xms5g" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{recalibration_report_filename} \
      --known-sites ~{dbsnp_vcf} \
      --known-sites ~{sep=" -known-sites " known_indels_sites_vcfs} \
      -L ~{sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "6 GiB"
    disks: "local-disk " + disk_size + " HDD"
    maxRetries: 3
    requested_memory_mb_per_core: 6000
    cpu: 1
    runtime_minutes: 120
  }
  output {
    File recalibration_report = "~{recalibration_report_filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  input {
    File input_bam
    File input_bam_index
    String output_bam_basename
    File recalibration_report
    Array[String] sequence_group_interval
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int compression_level
    Int bqsr_scatter
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.10.1"
    Int memory_multiplier = 1
    Int additional_disk = 20
    Boolean bin_base_qualities = true
    Boolean somatic = false
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil((size(input_bam, "GiB") * 3 / bqsr_scatter) + ref_size) + additional_disk

  Int memory_size = ceil(3500 * memory_multiplier)

  Boolean bin_somatic_base_qualities = bin_base_qualities && somatic

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command {
    gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=~{compression_level} -Xms3000m" \
      ApplyBQSR \
      --create-output-bam-md5 \
      --add-output-sam-program-record \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{output_bam_basename}.bam \
      -bqsr ~{recalibration_report} \
      ~{true='--static-quantized-quals 10' false='' bin_base_qualities} \
      ~{true='--static-quantized-quals 20' false='' bin_base_qualities} \
      ~{true='--static-quantized-quals 30' false='' bin_base_qualities} \
      ~{true='--static-quantized-quals 40' false='' bin_somatic_base_qualities} \
      ~{true='--static-quantized-quals 50' false='' bin_somatic_base_qualities} \
      -L ~{sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    maxRetries: 3
    memory: "~{memory_size} MiB"
    disks: "local-disk " + disk_size + " HDD"
    requested_memory_mb_per_core: 3000
    cpu: 1
    runtime_minutes: 90
  }
  output {
    File recalibrated_bam = "~{output_bam_basename}.bam"
    File recalibrated_bam_checksum = "~{output_bam_basename}.bam.md5"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  input {
    Array[File] input_bqsr_reports
    String output_report_filename
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.10.1"
  }

  command {
    gatk --java-options "-Xms3000m" \
      GatherBQSRReports \
      -I ~{sep=' -I ' input_bqsr_reports} \
      -O ~{output_report_filename}
    }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "3500 MiB"
    disks: "local-disk 20 HDD"
    requested_memory_mb_per_core: 3500
    cpu: 1
    runtime_minutes: 60
  }
  output {
    File output_bqsr_report = "~{output_report_filename}"
  }
}

# Combine multiple *sorted* BAM files
task GatherSortedBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename
    Float total_input_size
    Int compression_level
    Int preemptible_tries
  }

  # Multiply the input bam size by two to account for the input and output
  Int disk_size = ceil(2 * total_input_size) + 20

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms2000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=true
      
    samtools index ~{output_bam_basename}.bam
    
    cp ~{output_bam_basename}.bam.bai ~{output_bam_basename}.bai
    }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    preemptible: preemptible_tries
    maxRetries: 3
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
    requested_memory_mb_per_core: 4000
    cpu: 1
    runtime_minutes: 120
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}


task RecalibrateBAM {
  input {
  File input_bam
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  File dbsnp_vcf
  File dbsnp_vcf_index
  Array[File] known_indels_sites_vcfs
  Array[File] known_indels_sites_indices

  # Runtime params
  String gatk_path
  Int threads
  String docker
  }
  
  command {
  set -e
    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      BaseRecalibrator \
      -R ~{reference_fa} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{sample_basename}.sample.recal_data.grp \
      --known-sites ~{dbsnp_vcf} \
      --known-sites ~{sep=" -known-sites " known_indels_sites_vcfs}

      ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      ApplyBQSR \
      --create-output-bam-md5 \
      --add-output-sam-program-record \
      -R ~{reference_fa} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{sample_basename}.recalibrated.bam \
      -bqsr ~{sample_basename}.sample.recal_data.grp 
  }

  output {
    File output_bam = "~{sample_basename}.recalibrated.bam"
    File output_bai = "~{sample_basename}.recalibrated.bai"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
  }
}

task HaplotypeCaller {
  input {
  File input_bam
  File input_bai
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  String chromosome

  Boolean? make_bamout

  # Runtime params
  String gatk_path
  Int threads
  String docker
  }

  String chromosome_intervals = sub(chromosome, ",", " -L ")
  
  command {
    set -e   
    # We need to create this file, even if it stays empty
    touch bamout.bam
    
    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}" \
      HaplotypeCaller \
      -R ~{reference_fa} \
      -stand-call-conf 50.0 \
      -L ~{chromosome_intervals} \
      -I ~{input_bam} \
      -O ~{sample_basename}.raw.vcf.gz \
      ~{true='-bamout bamout.bam' false='' make_bamout}
  }

  output {
    File output_vcf = "~{sample_basename}.raw.vcf.gz"
    File output_vcf_index = "~{sample_basename}.raw.vcf.gz.tbi"
    File output_bamout = "bamout.bam"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
    requested_memory_mb_per_core: 6000
    cpu: 1
    runtime_minutes: 1200
  }
}

# Merge GVCFs generated per-interval for the same sample
task MergeVCFs {
  input {
    # Command parameters
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String sample_basename

    # Runtime parameters
    String gatk_path
    String docker
  }
  
  command {
  set -e

    ~{gatk_path} --java-options -Xmx4G  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{sample_basename}.raw.vcf.gz
  }
  runtime {
    docker: docker
    maxRetries: 3
    requested_memory_mb_per_core: 4000
    cpu: 1
    runtime_minutes: 60
  }
  output {
    File output_vcf = "~{sample_basename}.raw.vcf.gz"
    File output_vcf_index = "~{sample_basename}.raw.vcf.gz.tbi"
  }
}

task MergeBamOuts {
  input {
    File ref_fasta
    File ref_fai
    File ref_dict
    # Array[File]+ bam_outs
    Array[File] bam_outs
    String sample_basename
    String gatk_path
    String docker
  }

  command <<<
	set -e
	# This command block assumes that there is at least one file in bam_outs.
	#  Do not call this task if len(bam_outs) == 0  
	~{gatk_path} --java-options -Xmx4G GatherBamFiles -I ~{sep=" -I " bam_outs} -O unsorted.out.bam -R ~{ref_fasta}

	# We must sort because adjacent scatters may have overlapping (padded) assembly regions, hence overlapping bamouts
	~{gatk_path} --java-options -Xmx4G SortSam -I unsorted.out.bam -O bamout.bam --SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT
	~{gatk_path} --java-options -Xmx4G BuildBamIndex -I bamout.bam -VALIDATION_STRINGENCY LENIENT

	mv bamout.bam  ~{sample_basename}.bamout.bam
	mv bamout.bai  ~{sample_basename}.bamout.bai
  >>>


  runtime {
	    docker: docker
	    maxRetries: 3
	    requested_memory_mb_per_core: 4000
	    cpu: 1
	    runtime_minutes: 60 
  }
  output {
        File merged_bam_out = "~{sample_basename}.bamout.bam"
        File merged_bam_out_index = "~{sample_basename}.bamout.bai"
  }
}


task SplitSNPindel {
  input {
  File input_vcf
  File input_vcf_index
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  # Runtime params
  String gatk_path
  Int threads
  String docker
  }
  
  command {
  set -e
    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      SelectVariants \
      -R ~{reference_fa} \
      -V ~{input_vcf} \
      -O ~{sample_basename}.raw.snp.vcf \
      -select-type SNP

      ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      SelectVariants \
      -R ~{reference_fa} \
      -V ~{input_vcf} \
      -O ~{sample_basename}.raw.indel.vcf \
      -select-type INDEL

  }

  output {
    File output_vcf_SNP = "~{sample_basename}.raw.snp.vcf"
    File output_vcf_INDEL = "~{sample_basename}.raw.indel.vcf"
  }

  runtime {
    docker: "~{docker}"
    requested_memory_mb_per_core: 6000
    cpu: 1
    runtime_minutes: 60
    maxRetries: 3
  }
}

task VariantFiltrationSNP {
  input {
  File input_vcf
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  # Runtime params
  String gatk_path
  Int threads
  String docker
  }
  
  command {
  set -e
    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      VariantFiltration \
      -R ~{reference_fa} \
      -V ~{input_vcf} \
      -O ~{sample_basename}.filtered.snp.vcf \
      --filter-expression "DP < 5" --filter-name "LowCoverage" \
      --filter-expression "QUAL < 30.0" --filter-name "VeryLowQual" \
      --filter-expression "FS > 100.0" --filter-name "StrandBias" \
      --filter-expression "SOR>10.0" --filter-name "HighStrandsOR"
  }

  output {
    File output_vcf = "~{sample_basename}.filtered.snp.vcf"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
    requested_memory_mb_per_core: 6000
    cpu: 1
    runtime_minutes: 60
  }
}

task VariantFiltrationINDEL {
  input {
  File input_vcf
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  # Runtime params
  String gatk_path
  Int threads
  String docker
  }
  
  command {
  set -e
    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}" \
      VariantFiltration \
      -R ~{reference_fa} \
      -V ~{input_vcf} \
      -O ~{sample_basename}.filtered.indel.vcf \
      --filter-expression "DP < 5" --filter-name "LowCoverage" \
      --filter-expression  "QUAL < 30.0" --filter-name "VeryLowQual" \
      --filter-expression  "FS > 200.0" --filter-name "StrandBias" \
      --filter-expression  "SOR>10.0" --filter-name "HighStrandsOR"
  }

  output {
    File output_vcf = "~{sample_basename}.filtered.indel.vcf"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
    requested_memory_mb_per_core: 6000
    cpu: 1
    runtime_minutes: 60
  }
}

task CombineSNPindel {
  input {
  File input_snp_vcf
  File input_indel_vcf
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  # Runtime params
  String picard_path
  Int threads
  String docker
  }
  
  command {
  set -e
     java -jar ~{picard_path} MergeVcfs I=~{input_snp_vcf} I=~{input_indel_vcf} O=~{sample_basename}.both.vcf
  }

  output {
    File output_vcf = "~{sample_basename}.both.vcf"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
    requested_memory_mb_per_core: 6000
    cpu: 1
    runtime_minutes: 60
  }
}

task SelectFinalVariants {
  input {
  File input_vcf
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  # Runtime params
  String gatk_path
  Int threads
  String docker
  }
  
  command {
  set -e
     ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}" \
      SelectVariants \
      -R ~{reference_fa} \
      -V ~{input_vcf} \
      -O ~{sample_basename}.vcf \
      --exclude-filtered
  }

  output {
    File output_vcf = "~{sample_basename}.vcf"
    File output_vcf_index = "~{sample_basename}.vcf.idx"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
    requested_memory_mb_per_core: 6000
    cpu: 1
    runtime_minutes: 60
  }
}

# Convert a pair of FASTQs to uBAM
task PairedFastQsToUnmappedBAM {
  input {
    # Command parameters
    String sample_name
    File fastq_1
    File fastq_2
    String readgroup_name
    String library_name
    String platform_unit
    String run_date
    String platform_name
    String sequencing_center
    String gatk_path

    # Runtime parameters
    Int addtional_disk_space_gb = 10
    Int machine_mem_gb = 7
    Int preemptible_attempts = 3
    String docker
  }
    Int command_mem_gb = machine_mem_gb - 1
    Int disk_space_gb = ceil((size(fastq_1, "GB") + size(fastq_2, "GB")) * 2 ) + addtional_disk_space_gb
  command {
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}g" \
    FastqToSam \
    --FASTQ ~{fastq_1} \
    --FASTQ2 ~{fastq_2} \
    --OUTPUT ~{readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME ~{readgroup_name} \
    --SAMPLE_NAME ~{sample_name} \
    --LIBRARY_NAME ~{library_name} \
    --PLATFORM_UNIT ~{platform_unit} \
    --RUN_DATE ~{run_date} \
    --PLATFORM ~{platform_name} \
    --SEQUENCING_CENTER ~{sequencing_center} 
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries: 3
    requested_memory_mb_per_core: 6000
    cpu: 1
    runtime_minutes: 700
  }
  output {
    File output_unmapped_bam = "~{readgroup_name}.unmapped.bam"
  }
}

task SamSplitter {
  input {
    File input_bam
    Int n_reads
    Int preemptible_tries
    Int compression_level
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  # Since the output bams are less compressed than the input bam we need a disk multiplier that's larger than 2.
  Float disk_multiplier = 2.5
  Int disk_size = ceil(disk_multiplier * unmapped_bam_size + 20)

  command {
    set -e
    mkdir output_dir

    total_reads=$(samtools view -c ~{input_bam})

    java -Dsamjdk.compression_level=~{compression_level} -Xms3000m -jar /usr/gitc/picard.jar SplitSamByNumberOfReads \
      INPUT=~{input_bam} \
      OUTPUT=output_dir \
      SPLIT_TO_N_READS=~{n_reads} \
      TOTAL_READS_IN_INPUT=$total_reads
  }
  output {
    Array[File] split_bams = glob("output_dir/*.bam")
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    preemptible: preemptible_tries
    maxRetries: 3
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
    requested_memory_mb_per_core: 6000
    cpu: 1
    runtime_minutes: 600
  }
}

# Combine multiple *unsorted* BAM files
# Note that if/when WDL supports optional outputs, we should merge this task with the sorted version
task GatherUnsortedBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename
    Int compression_level
    Int preemptible_tries
  }

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms2000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=false
    }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    preemptible: preemptible_tries
    maxRetries: 3
    requested_memory_mb_per_core: 6000
    cpu: 1
    runtime_minutes: 400
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}

# Sort BAM file by coordinate order
task SortSam {
  input {
    File input_bam
    String output_bam_basename
    Int preemptible_tries
    Int compression_level
  }
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Float sort_sam_disk_multiplier = 6.25
  Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms8000m -jar /usr/gitc/picard.jar \
      SortSam \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=1000000
    
    samtools index ~{output_bam_basename}.bam
    
    cp ~{output_bam_basename}.bam.bai ~{output_bam_basename}.bai
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    memory: "5000 MiB"
    preemptible: preemptible_tries
    maxRetries: 3
    requested_memory_mb_per_core: 1000
    cpu: 45
    runtime_minutes: 900
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}



task StringToArray {
  input {
    String input_string
    String separator
  }
  command <<<
    echo '~{input_string}' | tr '~{separator}' \\n | tr -d "[:blank:]" > intervals.list
    echo '~{input_string}' | tr '~{separator}' \\n | tr -d "[:blank:]"
  >>>
  runtime {
    docker:"biocontainers/bcftools:v1.9-1-deb_cv1"
    requested_memory_mb_per_core: 500
    cpu: 1
    runtime_minutes: 5
  }
  output {
    Array[String] values = read_lines(stdout())
    File intervals_list = "intervals.list"
  }
}
