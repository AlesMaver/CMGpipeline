version 1.0

import "../CRAM_conversions.wdl" as CramConversions
import "../manta/manta_workflow.wdl" as Manta

# WORKFLOW DEFINITION 
workflow SoftSearchWF {
  input {
    File? input_bam
    File? input_bam_index
    File? input_cram
    File? input_cram_index
    File reference_fa
    File reference_fai
    File reference_dict
    String sample_basename
  }

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

  call FileToArray
  
  call LocaliseInput as LocaliseInput {
    input:
      input_bam = select_first([input_bam, CramToBam.output_bam]),
      input_bam_index = select_first([input_bam_index, CramToBam.output_bai])
  }

  scatter (chromosome in FileToArray.scatter_regions ) {
    call SoftSearch {
      input:
        input_bam = LocaliseInput.output_bam,
        input_bam_index = LocaliseInput.output_bam_index,
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

  # call MergeVCFs as MergeVCFs {
  #   input:
  #     input_vcfs = CompressAndIndexVCF.output_vcfgz,
  #     input_vcfs_indexes = CompressAndIndexVCF.output_vcfgz_index,
  #     reference_dict = reference_dict,
  #     sample_basename = sample_basename,
  #     docker = "broadinstitute/gatk:4.2.0.0",
  #     gatk_path = "/gatk/gatk"
  # }
  
  call MergeVCFsBcftools as MergeVCFsBcftools {
    input:
      input_vcfs = CompressAndIndexVCF.output_vcfgz,
      input_vcfs_indexes = CompressAndIndexVCF.output_vcfgz_index,
      reference_dict = reference_dict,
      sample_basename = sample_basename
  }

  call SoftSearch_filter as SoftSearch_filter {
    input:
      input_vcf = MergeVCFsBcftools.output_vcf,
      sample_basename = sample_basename
  }

  call Manta.annotSV as SoftSearch_annotSV {
    input:
      genome_build = "GRCh37",
      input_vcf = SoftSearch_filter.output_vcf,
      output_tsv_name = sample_basename + ".SoftSearch.annotSV.tsv"
  }

  output {
    File output_complete_vcf = MergeVCFsBcftools.output_vcf
    File output_complete_vcf_index = MergeVCFsBcftools.output_vcf_index
    File output_vcf = SoftSearch_filter.output_vcf
    File? output_tsv_name = SoftSearch_annotSV.sv_variants_tsv
  }
}


# Localise input BAM file to prevent multiple copy events for several shards (when inputs are on a different block device)
task LocaliseInput {
  input {
    File input_bam
    File input_bam_index
  }

  String output_bam = basename(input_bam)
  String output_bam_index = basename(input_bam_index)

  command <<<
    cp ~{input_bam} ./
    cp ~{input_bam_index} ./
  >>>

  runtime {
    docker: "alesmaver/softsearch"
    maxRetries: 3
    requested_memory_mb_per_core: 500
    cpu: 1
    runtime_minutes: 20
  }
  output {
    File output_bam = "~{output_bam}"
    File output_bam_index = "~{output_bam_index}"
  }
}

task FileToArray {

    String github_fname = "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/softsearch/wgs_calling_regions.v1_mod.list.txt"
    String fname = "wgs_calling_regions.v1_mod.list.txt"

    command <<<
        # Get chromosome interval list file
        wget  ~{github_fname}
        # Re-attempt the wget without the https_proxy set
        unset https_proxy 
        wget  ~{github_fname}     
        
        cat ~{fname}
    >>>
      
    runtime {
        docker: "alesmaver/softsearch"
        runtime_minutes: 5
    }
    output {
        Array[String] scatter_regions = read_lines(stdout())
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

    # Create an empty vcf file, so that the workflow continues in case a single shard fails (usually recurring problematic shards)
    touch ~{sample_basename}.softSearch.vcf

    # Get population masked regions from the softsearch repository
    ## wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/svcalling/softsearch/breakpoint_mask.bed
    wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/softsearch/breakpoint_mask.bed
    wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/softsearch/highly_variable_regions.bed
    # Re-attempt the wget without the https_proxy set
    unset https_proxy 
    ## wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/svcalling/softsearch/breakpoint_mask.bed
    wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/softsearch/breakpoint_mask.bed
    wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/softsearch/highly_variable_regions.bed

    # Merge population and softsearch masks
    cat breakpoint_mask.bed /softsearch/library/blacklist_fixed.bed highly_variable_regions.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > blacklist.bed

    # Perform SoftSearch
    perl /softsearch/script/SoftSearch.pl -b ~{input_bam} -o ~{sample_basename}.softSearch.vcf -f ~{ref_fasta} -v -blacklist blacklist.bed -genome hg19.genome -c ~{chromosome}

  >>>

  runtime {
    docker: "alesmaver/softsearch"
    maxRetries: 3
    continueOnReturnCode: true
    requested_memory_mb_per_core: 1000
    cpu: 5
    runtime_minutes: 600
  }
  output {
    File output_vcf = "~{sample_basename}.softSearch.vcf"
  }
}

task CompressAndIndexVCF {
  input {
    File input_vcf
    String sample_basename
    String docker
  }

  command <<<
    set -e

    if [ -s "~{input_vcf}" ]; then
        echo "Input file is not empty. Proceeding normally..."
    else
        echo "Input file is empty. Creating an output file with default content."
        INPUT_VCF="~{input_vcf}"

        # Redirect printf output to the file
        printf "##fileformat=VCFv4.1\n" > "$INPUT_VCF"
        printf "##FILTER=<ID=PASS,Description=\"All filters passed\">\n" >> "$INPUT_VCF"
        printf "##fileDate=$(date +"%Y-%m-%d %H:%M:%S")\n" >> "$INPUT_VCF"
        printf "##source=SoftSearch.pl\n" >> "$INPUT_VCF"
        printf "##bcftools_viewVersion=1.14+htslib-1.14\n" >> "$INPUT_VCF"
        printf "##bcftools_viewCommand=view test.vcf; Date=$(date)\n" >> "$INPUT_VCF"
        printf "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t4\n" >> "$INPUT_VCF"
    fi

    # sorting the variants part of the input file
    cat ~{input_vcf} | bcftools view -h  > ~{sample_basename}.header.vcf 
    cat ~{input_vcf} | bcftools view -H | sort -k1,1V -k2,2n > ~{sample_basename}.variants.vcf
    cat ~{sample_basename}.header.vcf ~{sample_basename}.variants.vcf > ~{sample_basename}.sorted.vcf
  
    bcftools view -Oz ~{sample_basename}.sorted.vcf > ~{sample_basename}.vcf.gz
    bcftools index -t ~{sample_basename}.vcf.gz
  >>>
  
  runtime {
    docker: docker
    maxRetries: 3
    requested_memory_mb_per_core: 1000
    cpu: 5
    runtime_minutes: 90
  }
  output {
    File output_vcfgz = "~{sample_basename}.vcf.gz"
    File output_vcfgz_index = "~{sample_basename}.vcf.gz.tbi"
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
    requested_memory_mb_per_core: 1000
    cpu: 4
    runtime_minutes: 1800
  }
  output {
    File output_vcf = "~{sample_basename}.softSearch.vcf.gz"
    File output_vcf_index = "~{sample_basename}.softSearch.vcf.gz.tbi"
  }
}

# Merge GVCFs generated per-interval for the same sample using bcftools
task MergeVCFsBcftools {
  input {
    # Command parameters
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String sample_basename
    File reference_dict

    # Runtime parameters
    String bcftools_path = "bcftools"
    String docker = "biocontainers/bcftools:v1.9-1-deb_cv1"
  }
  
  command <<<
    set -e
    ~{bcftools_path} concat --naive -Oz -o ~{sample_basename}.softSearch.merged.vcf.gz ~{sep=' ' input_vcfs}

    # SoftSearch produces an incomplete header, so we need to add the contigs from the reference dictionary
    awk '/^@SQ/ {print "##contig=<ID="$2",length="$3">"}' ~{reference_dict} | sed 's/SN://; s/LN://' > contigs_header.txt
    bcftools view -h  ~{sample_basename}.softSearch.merged.vcf.gz | grep -v "^#CHROM" > original_header.txt
    cat original_header.txt contigs_header.txt > new_header.txt
    bcftools view -h ~{sample_basename}.softSearch.merged.vcf.gz | grep "^#CHROM" >> new_header.txt

    bcftools view ~{sample_basename}.softSearch.merged.vcf.gz -Ov  | bcftools reheader -h new_header.txt | sed 's/MATEMATEID/MATEID/g' | sed 's/ch\([0-9]\+\)/chr\1/g' | bcftools view -Oz -o reheadered.vcf.gz
    bcftools annotate -x INFO/END -o cleaned_reheadered.vcf.gz -Oz reheadered.vcf.gz
    bcftools sort cleaned_reheadered.vcf.gz -o sorted_cleaned_reheadered.vcf.gz -Oz

    mv sorted_cleaned_reheadered.vcf.gz ~{sample_basename}.softSearch.vcf.gz
    ~{bcftools_path} index -t ~{sample_basename}.softSearch.vcf.gz
  >>>

  runtime {
    docker: docker
    maxRetries: 3
    requested_memory_mb_per_core: 1000
    cpu: 4
    runtime_minutes: 180
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
    requested_memory_mb_per_core: 1000
    cpu: 5
    runtime_minutes: 30
  }
  output {
    File output_vcf = "~{sample_basename}.softSearch.filtered.vcf.gz"
  }
}

