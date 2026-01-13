version 1.0

import "https://raw.githubusercontent.com/UW-GAC/primed-file-conversion/main/liftover_vcf_picard.wdl" as liftover_vcf_picard
import "http://wdl_server/bravo-pipeline/BravoDataPreparation.wdl" as BravoDataPreparation
import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/AnnotationPipeline.wdl" as AnnotationPipeline
workflow DownloadAndConvertRegeneron {
    input {
        Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX", "chrY"]
        String baseUrl = "test"

        String chain_url
        File target_fasta
        String out_prefix
        Int? mem_gb
    }

    scatter (chromosome in chromosomes) {
        call DownloadFile {
            input:
                chromosome = chromosome,
                baseUrl = baseUrl
        }

        call RenameChrs {
            input:
                input_vcf = DownloadFile.downloadedFile
        }

      call picard {
        input: vcf_file = RenameChrs.output_vcf,
               chain_url = chain_url,
               target_fasta = target_fasta,
               out_prefix = out_prefix,
               mem_gb = mem_gb
      }

        # call liftover_vcf_picard.liftover_vcf as liftover {
        #     input:
        #         vcf_file = RenameChrs.output_vcf,
        #         chain_url = chain_url,
        #         target_fasta = target_fasta,
        #         out_prefix = out_prefix,
        #         mem_gb = mem_gb
        # }

        call SortAndIndexVcf {
            input:
                input_vcf = picard.out_file
        }

    }

  # Concatenate VCFs 
  call concatSortVcf as concatSortVcf {
    input:
      input_vcfs = SortAndIndexVcf.output_vcf,
      input_vcfs_indices = SortAndIndexVcf.output_vcf_index
  }

  output {
    File output_vcf = concatSortVcf.output_vcf
    File topmed_output_vcf_index = concatSortVcf.output_vcf_index
  }

}

task DownloadFile {
    input {
        String chromosome
        String baseUrl
    }

    String chromosome_nochr = sub(chromosome, "chr", "")

    command {
        echo Downloading $chromosome;
        curl 'https://rgc-research.regeneron.com/me/downloads/20231004/rgc_me_variant_frequencies_~{chromosome}_20231004.vcf.gz' -H 'Accept-Encoding: gzip, deflate, br' --compressed > ~{chromosome}.vcf.gz
        curl 'https://rgc-research.regeneron.com/me/downloads/20231004/rgc_me_variant_frequencies_~{chromosome}_20231004.vcf.gz.tbi' -H 'Accept-Encoding: gzip, deflate, br' --compressed > ~{chromosome}.vcf.gz.tbi
    }

    output {
        File downloadedFile = "~{chromosome}.vcf.gz"
    }

    runtime {
        docker: "dceoy/bedops"
    }
}

task RenameChrs {
    input {
        File input_vcf
    }

    String file_basename = basename(input_vcf)

    command {
        wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/0490e00417c09d1d988fe4ca3fd2e3659b24e067/references/rename_chrs
        bcftools annotate --rename-chrs rename_chrs ~{input_vcf} -Oz -o renamed.~{file_basename}
    }

    output {
        File output_vcf = "renamed.~{file_basename}"
    }

    runtime {
        docker: "alesmaver/bcftools"
    }
}

task SortAndIndexVcf {
    input {
        String input_vcf
    }

    command {
        bcftools sort ~{input_vcf} -Oz -o sorted.vcf.gz
        bcftools index -t sorted.vcf.gz
    }

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
      String output_name ="output"
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
