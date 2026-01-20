version 1.0

### import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/AnnotationPipeline.wdl" as Annotation
import "../manta/manta_workflow.wdl" as Manta

# WORKFLOW DEFINITION 
workflow SoftSearchWF {
  input {
    File input_bam
    File input_bam_index
    File reference_fa
    File reference_fai
    File reference_dict
    Array[String] scatter_regions
    String sample_basename
  }
  
  scatter (chromosome in scatter_regions ) {
    call SoftSearch {
      input:
        input_bam = input_bam,
        input_bam_index = input_bam_index,
        ref_fasta=reference_fa,
        ref_fasta_index=reference_fai,
        chromosome=chromosome,
        sample_basename=sample_basename
    }

    call CompressAndIndexVCF {
      input:
        input_vcf = SoftSearch.output_vcf,
        sample_basename = sample_basename,
        docker = "dceoy/bcftools"
      }
  }

  call MergeVCFs as MergeVCFs {
    input:
      input_vcfs = CompressAndIndexVCF.output_vcfgz,
      input_vcfs_indexes = CompressAndIndexVCF.output_vcfgz_index,
      reference_dict = reference_dict,
      sample_basename = sample_basename,
      docker = "broadinstitute/gatk:4.2.0.0",
      gatk_path = "/gatk/gatk"
  }
  
  call SoftSearch_filter as SoftSearch_filter {
    input:
      input_vcf = MergeVCFs.output_vcf,
      sample_basename = sample_basename
  }

  call Manta.annotSV as SoftSearch_annotSV {
    input:
      genome_build = "GRCh37",
      input_vcf = SoftSearch_filter.output_vcf,
      output_tsv_name = sample_basename + ".SoftSearch.annotSV.tsv"
  }

  output {
    File output_complete_vcf = MergeVCFs.output_vcf
    File output_complete_vcf_index = MergeVCFs.output_vcf_index
    File output_vcf = SoftSearch_filter.output_vcf
    File? output_tsv_name = SoftSearch_annotSV.sv_variants_tsv
  }
}

# SoftSearch task
task SoftSearch {
  input {
    File input_bam
    File input_bam_index
    File ref_fasta
    File ref_fasta_index
    String chromosome
    String sample_basename
  }

  command <<<
    echo chromosome: ~{chromosome}
    
    # Generate the genome file
    awk -v OFS='\t' {'print $1,$2'} ~{ref_fasta_index} > hg19.genome

    # Get population masked regions from the softsearch repository
    ## wget --no-check-certificate https://raw.githubusercontent.com/AlesMaver/CMGpipeline/svcalling/softsearch/breakpoint_mask.bed
    wget --no-check-certificate https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/softsearch/breakpoint_mask.bed
    # Re-attempt the wget without the https_proxy set
    unset https_proxy 
    ## wget --no-check-certificate https://raw.githubusercontent.com/AlesMaver/CMGpipeline/svcalling/softsearch/breakpoint_mask.bed
    wget --no-check-certificate https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/softsearch/breakpoint_mask.bed

    # Merge population and softsearch masks
    cat breakpoint_mask.bed /softsearch/library/blacklist_fixed.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > blacklist.bed

    # Perform SoftSearch
    perl /softsearch/script/SoftSearch.pl -b ~{input_bam} -o ~{sample_basename}.softSearch.vcf -f ~{ref_fasta} -v -blacklist blacklist.bed -genome hg19.genome -c ~{chromosome}

  >>>

  runtime {
    docker: "alesmaver/softsearch"
    maxRetries: 3
    requested_memory_mb_per_core: 5000
    cpu: 2
    runtime_minutes: 600
  }
  output {
    File output_vcf = "~{sample_basename}.softSearch.vcf"
  }
}

# Merge GVCFs generated per-interval for the same sample
task MergeVCFs {
  input {
    # Command parameters
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String sample_basename
    File? reference_dict

    # Runtime parameters
    String gatk_path
    String docker
  }
  
  command {
    set -e
    ~{gatk_path} --java-options -Xmx4G  \
      SortVcf \
      ~{"--SEQUENCE_DICTIONARY " + reference_dict} \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{sample_basename}.softSearch.vcf.gz
  }
  runtime {
    docker: docker
    maxRetries: 3
    requested_memory_mb_per_core: 4000
    cpu: 1
    runtime_minutes: 60
  }
  output {
    File output_vcf = "~{sample_basename}.softSearch.vcf.gz"
    File output_vcf_index = "~{sample_basename}.softSearch.vcf.gz.tbi"
  }
}

task SoftSearch_filter {
  input {
    File input_vcf
    String sample_basename
    String docker = "biocontainers/bcftools:v1.9-1-deb_cv1"
  }
  
  command {
  set -e
    echo Filtering soft search vcf file
    bcftools view -i'FORMAT/nSC>5 && FORMAT/lSC>10' -Oz -o ~{sample_basename}.softSearch.filtered.vcf.gz ~{input_vcf}
    
  }
  runtime {
    docker: docker
    maxRetries: 3
    requested_memory_mb_per_core: 5000
    cpu: 1
    runtime_minutes: 30
  }
  output {
    File output_vcf = "~{sample_basename}.softSearch.filtered.vcf.gz"
  }
}


task CompressAndIndexVCF {
  input {
    File input_vcf
    String sample_basename
    String docker
  }

  command {
    set -e
    # sorting the variants part of the input file
    cat ~{input_vcf} | bcftools view -h  > ~{sample_basename}.header.vcf 
    cat ~{input_vcf} | bcftools view -H | sort -k1,1V -k2,2n > ~{sample_basename}.variants.vcf
    cat ~{sample_basename}.header.vcf ~{sample_basename}.variants.vcf > ~{sample_basename}.sorted.vcf
  
    bcftools view -Oz ~{sample_basename}.sorted.vcf > ~{sample_basename}.vcf.gz
    bcftools index -t ~{sample_basename}.vcf.gz
  }
  runtime {
    docker: docker
    maxRetries: 3
    requested_memory_mb_per_core: 5000
    cpu: 1
    runtime_minutes: 90
  }
  output {
    File output_vcfgz = "~{sample_basename}.vcf.gz"
    File output_vcfgz_index = "~{sample_basename}.vcf.gz.tbi"
  }
}

