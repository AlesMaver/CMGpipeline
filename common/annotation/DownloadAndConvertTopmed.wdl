version 1.0

import "https://raw.githubusercontent.com/UW-GAC/primed-file-conversion/main/liftover_vcf_picard.wdl" as liftover_vcf_picard
import "http://wdl_server/bravo-pipeline/BravoDataPreparation.wdl" as BravoDataPreparation
import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/AnnotationPipeline.wdl" as AnnotationPipeline
workflow DownloadAndConvertTopmed {
    input {
        Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"]
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

        call liftover_vcf_picard.liftover_vcf as liftover {
            input:
                vcf_file = DownloadFile.downloadedFile,
                chain_url = chain_url,
                target_fasta = target_fasta,
                out_prefix = out_prefix,
                mem_gb = mem_gb
        }

        call SortAndIndexVcf {
            input:
                input_vcf = liftover.out_file
        }

    }

  # Concatenate VCFs 
  call concatSortVcf as concatSortVcf {
    input:
      input_vcfs = SortAndIndexVcf.output_vcf,
      input_vcfs_indices = SortAndIndexVcf.output_vcf_index
  }

  output {
    File topmed_output_vcf = concatSortVcf.output_vcf
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
        echo Downloading TOPMED freeze 8 $chromosome;
        curl 'https://bravo.sph.umich.edu/freeze8/hg38/downloads/vcf/~{chromosome_nochr}' -H 'Accept-Encoding: gzip, deflate, br' -H 'Cookie: _gid=GA1.2.851051206.1693377892;um_cookie_consent=na;_gcl_au=1.1.324996009.1693413116;_ga_4L4323B4K3=GS1.1.1693413115.1.0.1693413115.60.0.0;_ga_NVLEYCF35H=GS1.1.1693412983.1.1.1693413305.0.0.0;_gat=1;_ga_29C4RBNCNZ=GS1.1.1693493731.5.0.1693493737.0.0.0;_gat_gtag_UA_73910830_2=1;remember_token=ales.maver@gmail.com|cdfe3abffa9cda418cb03b54953435866e80628856377ac2a1821b52ec5c5e3fd3c59cb829f750fbfd443f26a7ddaeb777e7a769a9bcbece9c704289ba8497b0;session=.eJwdj7tOxTAQBf_FdRTZ69c6FSJKR0XBq4nW9vrmooSL7AQKxL8TUY9mdM6PmEvltohhrwd3Yr5mMQjOUXMo2Slw1mKMnMDm4lB5pyNFyCZGDBzYSGex-KSChWJTORFGExwkcpSk9kBEKoFhp_E0XMJkKErKqBCytxJC0hqAyTFkFcB7i6ITlTfeIte5cbp95CYGdEbKXnai7bTzufK1jfqJvtvDNE4vi9-fx_Xzkd6WaX2v98fZONrp_x-ilVu_0RfXu8tG17VPt038_gEQCE0L.ZPCp8Q.6HBM6ErfzXN36Aag406Tp2JusNU;_ga=GA1.1.1277966882.1650988017;_ga_HD76LS6C66=GS1.1.1693493737.140.1.1693493748.0.0.0;_ga_5B596RM26L=GS1.1.1693493751.5.0.1693493756.0.0.0;' --compressed > ~{chromosome_nochr}.BRAVO_TOPMed_Freeze_8.vcf.gz
    }

    output {
        File downloadedFile = "~{chromosome_nochr}.BRAVO_TOPMed_Freeze_8.vcf.gz"
    }

    runtime {
        docker: "dceoy/bedops"
    }
}

task SortAndIndexVcf {
    input {
        String input_vcf
    }

    command {
        bcftools sort ~{input_vcf} -Oz -o TOPMED.sorted.vcf.gz
        bcftools index -t TOPMED.sorted.vcf.gz
    }

    output {
        File output_vcf = "TOPMED.sorted.vcf.gz"
        File output_vcf_index = "TOPMED.sorted.vcf.gz.tbi"
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