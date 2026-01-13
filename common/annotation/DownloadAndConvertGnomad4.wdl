version 1.0

workflow DownloadAndConvertGnomad4 {
    input {
        # Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX", "chrY"]
        Array[String] chromosomes = ["chr20","chr21", "chrY"]

        String baseUrl = "test"

        String chain_url
        File target_fasta
        String out_prefix
        Int? mem_gb
    }

    scatter (chromosome in chromosomes ) {
        call DownloadFile {
            input:
                chromosome = chromosome,
                baseUrl = baseUrl
        }
    }

    call FileToArray

    scatter (chromosome in FileToArray.scatter_regions ) {
        call FilterVariants {
            input:
                input_vcfs = DownloadFile.output_vcf,
                input_vcfs_indices = DownloadFile.output_vcf_index,
                contig = chromosome, 
                output_vcf = "scattered.vcf.gz"

        }

        call picard {
            input:
                vcf_file = FilterVariants.output_vcf,
                chain_url = chain_url,
                target_fasta = target_fasta,
                out_prefix = out_prefix,
                mem_gb = mem_gb
        }

        call SortAndIndexVcf {
            input:
                input_vcf = picard.out_file
        }
    }

    scatter (chromosome in chromosomes) {
        call FilterVariantsChromosome {
            input:
                input_vcfs = SortAndIndexVcf.output_vcf,
                input_vcfs_indices = SortAndIndexVcf.output_vcf_index,
                contig = chromosome,
                output_vcf = "~{chromosome}_filtered.vcf.gz"
        }

        call ConcatenateAndSort {
            input:
                input_vcf = FilterVariantsChromosome.output_vcf,
                output_vcf = "~{chromosome}_sorted.vcf.gz"
        }
    }

    call ConcatenateAndIndex {
        input:
            input_vcfs = ConcatenateAndSort.output_vcf,
            input_vcfs_indices = ConcatenateAndSort.output_vcf_index,
            output_vcf = "~{out_prefix}.vcf.gz"
    }

  # # Concatenate VCFs 
  # call concatSortVcf as concatSortVcf {
  #   input:
  #     input_vcfs = SortAndIndexVcf.output_vcf,
  #     input_vcfs_indices = SortAndIndexVcf.output_vcf_index, 
  #     output_name = out_prefix
  # }

  output {
    File output_vcf = ConcatenateAndIndex.output_vcf
    File output_vcf_index = ConcatenateAndIndex.output_vcf_index
  }

}

task DownloadFile {
    input {
        String chromosome
        String baseUrl
    }

    String chromosome_nochr = sub(chromosome, "chr", "")

    command {
        echo Downloading ~{chromosome};
        curl 'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.~{chromosome}.vcf.bgz' -H 'Accept-Encoding: gzip, deflate, br' --compressed > ~{chromosome}.vcf.bgz
        curl 'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/exomes/gnomad.exomes.v4.0.sites.~{chromosome}.vcf.bgz.tbi' -H 'Accept-Encoding: gzip, deflate, br' --compressed > ~{chromosome}.vcf.bgz.tbi
        
    }

    output {
        File output_vcf = "~{chromosome}.vcf.bgz"
        File output_vcf_index = "~{chromosome}.vcf.bgz.tbi"
    }

    runtime {
        docker: "dceoy/bedops"
    }
}


task SortAndIndexVcf {
    input {
        String input_vcf
    }

    command <<<
        bcftools annotate -x ^INFO/AC_joint,^INFO/AF_joint,^INFO/AN_joint,^INFO/nhomalt_joint ~{input_vcf} -Oz -o thin.vcf
        bcftools sort thin.vcf -Oz -o sorted.vcf.gz
        bcftools index -t sorted.vcf.gz
    >>>

    output {
        File output_vcf = "sorted.vcf.gz"
        File output_vcf_index = "sorted.vcf.gz.tbi"
    }

    runtime {
        docker: "alesmaver/bcftools"
    }
}

##############################
task SortVcf {
    input {
      File input_vcf
      File input_vcf_indices
      String output_name ="output"
    }
  
  command <<<
  set -e
    mkdir $PWD/sort_tmp
    bcftools sort ~{input_vcf} -Oz -o ~{output_name}.vcf.gz --temp-dir $PWD/sort_tmp -m 9000000000
    bcftools index --threads 10 -t ~{output_name}.vcf.gz
  >>>

  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    requested_memory_mb_per_core: 1000
    cpu: 10
    #runtime_minutes: 90
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
    File output_vcf_index = "~{output_name}.vcf.gz.tbi"
  }
}

##############################
task concatSortVcf {
    input {
      Array[File] input_vcfs
      Array[File] input_vcfs_indices
      String output_name = "output"
    }
  
  command <<<
  set -e
    mkdir $PWD/sort_tmp
    bcftools concat -a --threads 10 -f ~{write_lines(input_vcfs)} -Oz -o ~{output_name}_unsorted.vcf.gz
    bcftools sort ~{output_name}_unsorted.vcf.gz -Oz -o ~{output_name}.vcf.gz --temp-dir $PWD/sort_tmp -m 9000000000
    bcftools index --threads 10 -t ~{output_name}.vcf.gz
  >>>

  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    requested_memory_mb_per_core: 1000
    cpu: 10
    #runtime_minutes: 90
  }
  output {
    File output_vcf = "~{output_name}.vcf.gz"
    File output_vcf_index = "~{output_name}.vcf.gz.tbi"
  }
}


task picard {
    input {
        File vcf_file
        String chain_url
        File target_fasta
        String out_prefix
        Int mem_gb = 16
    }

    String chain_file = basename(chain_url)

    command <<<
        curl ~{chain_url} --output ~{chain_file}
        java -Xmx~{mem_gb}g -jar /usr/picard/picard.jar CreateSequenceDictionary \
            --REFERENCE ~{target_fasta}
        java -Xmx~{mem_gb}g -jar /usr/picard/picard.jar LiftoverVcf \
            --CHAIN ~{chain_file} \
            --INPUT ~{vcf_file} \
            --OUTPUT ~{out_prefix}.vcf.gz \
            --REJECT rejected_variants.vcf.gz \
            --REFERENCE_SEQUENCE ~{target_fasta} \
            --RECOVER_SWAPPED_REF_ALT true \
            --ALLOW_MISSING_FIELDS_IN_HEADER true \
            --MAX_RECORDS_IN_RAM 10000
        zcat rejected_variants.vcf.gz | grep -v "^#" | wc -l > num_rejects.txt
    >>>

    output {
        File out_file = "~{out_prefix}.vcf.gz"
        File rejects_file = "rejected_variants.vcf.gz"
        Int num_rejects = read_int("num_rejects.txt")
    }

    runtime {
        docker: "broadinstitute/picard:2.27.5"
        memory: "~{mem_gb}GB"
    }
}

task FilterVariants {
    input {
        Array[String] input_vcfs
        Array[String] input_vcfs_indices
        String contig
        String output_vcf
    }

    command <<<
        bcftools concat -a -r ~{contig} ~{sep=" " input_vcfs} | bcftools annotate -x ^INFO/AC_joint,^INFO/AF_joint,^INFO/AN_joint,^INFO/nhomalt_joint -Oz -o concatenated_~{output_vcf}
        bcftools sort concatenated_~{output_vcf} -Oz -o sorted_~{output_vcf}
        bcftools index -t sorted_~{output_vcf}
    >>>

    runtime {
      docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
      requested_memory_mb_per_core: 1000
      cpu: 10
  }
  output {
    File output_vcf = "sorted_~{output_vcf}"
    File output_vcf_index = "sorted_~{output_vcf}.tbi"
  }
}

task FilterVariantsChromosome {
    input {
        Array[String] input_vcfs
        Array[String] input_vcfs_indices
        String contig
        String output_vcf
    }

    command <<<
        bcftools concat -a -r ~{contig} ~{sep=" " input_vcfs} -Oz -o concatenated_~{output_vcf}
        bcftools sort concatenated_~{output_vcf} -Oz -o sorted_~{output_vcf}
        bcftools index -t sorted_~{output_vcf}
    >>>

    runtime {
      docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
      requested_memory_mb_per_core: 1000
      cpu: 10
  }
  output {
    File output_vcf = "sorted_~{output_vcf}"
    File output_vcf_index = "sorted_~{output_vcf}.tbi"
  }
}

task ConcatenateAndSort {
    input {
        String input_vcf
        String output_vcf
    }

    command <<<
        bcftools concat ~{input_vcf} -O z -o ~{output_vcf} && \
        bcftools sort ~{output_vcf} -O z -o ~{output_vcf} && \
        bcftools index -t ~{output_vcf}
    >>>

    runtime {
      docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
      requested_memory_mb_per_core: 1000
      cpu: 10
    }

    output {
      File output_vcf = "~{output_vcf}"
      File output_vcf_index = "~{output_vcf}.tbi"
    }
}

task ConcatenateAndIndex {
    input {
        Array[File] input_vcfs
        Array[File] input_vcfs_indices
        String output_vcf
    }

    command <<<
        bcftools concat --naive --threads 10 ~{sep=" " input_vcfs} -O z -o ~{output_vcf}
        bcftools index -t ~{output_vcf}
    >>>

    runtime {
      docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
      requested_memory_mb_per_core: 1000
      cpu: 10
    }

    output {
      File output_vcf = "~{output_vcf}"
      File output_vcf_index = "~{output_vcf}.tbi"
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
