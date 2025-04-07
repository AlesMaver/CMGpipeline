version 1.0

# VEP wdl: scattered version

#--------------------------
# v delu... in progress...
#--------------------------

workflow VEP {
  input {
    String sample_basename
    File input_vcf
    # not needed File input_vcf_index
    String filename_infix = ""
    ### tole uporabi pri klicu call VEP v anotacijah: output_vcf = sample_basename + ".VEP.hg19.annotated.vcf.gz" | sample_basename + ".DeepVariant.VEP.hg19.annotated.vcf.gz" 
    File chromosome_list
    String? targetRegions
  }

  String filename_suffix = ".VEP.hg19.annotated.vcf.gz"

  Array[String] chromosomes = read_lines(chromosome_list)

  call VcfZippingAndIndexing {
        input:
	    input_vcf = input_vcf
  }

  if( defined(targetRegions) ) {
      call StringToArray {
        input:
            input_string = select_first([targetRegions, ""]),
            separator = ";"
      }
  }

  scatter (chromosome in select_first([StringToArray.values, chromosomes]) ) {

    call VcfPartitioning {
        input:
            sample_basename = sample_basename,
	    input_vcf = VcfZippingAndIndexing.output_vcf,
	    input_vcf_index = VcfZippingAndIndexing.output_vcf_index,
	    chromosome = chromosome
    }

    call RunVEP {
        input:
            sample_basename = sample_basename,
            input_vcf = VcfPartitioning.output_vcf,
            annotated_vcf = sample_basename + filename_infix + filename_suffix,
	    chromosome = chromosome
    }

  }  # end-of-scatter

  call MergeVCFs {
    input:
      input_vcfs = RunVEP.output_vcf,
      input_vcfs_indexes = RunVEP.output_vcf_index,
      sample_basename = sample_basename,
      output_filename_gz = sample_basename + filename_infix + filename_suffix
  }

  # workflow VEP output files:
  output {
      #File output_vcf = MergeVCFs.output_vcfgz
      #File output_vcf_index = MergeVCFs.output_vcfgz_index
  }

}

##################
# TASK DEFINITIONS
##################


task VcfZippingAndIndexing {
  input {
    File input_vcf
  }

  String file_name = basename(input_vcf)
  command {
    set -e
    bcftools sort ~{input_vcf} -Oz -o ~{file_name}.gz
    #bgzip -c ~{input_vcf} > ~{file_name}.gz
    bcftools index -t ~{file_name}.gz
  }
  runtime {
    docker: "dceoy/bcftools"
    maxRetries: 1
    requested_memory_mb_per_core: 1000
    cpu: 4
    runtime_minutes: 60
  }
  output {
    File output_vcf = "~{file_name}.gz"
    File output_vcf_index = "~{file_name}.gz.tbi"
  }
}

task VcfPartitioning {
  input {
    String sample_basename
    File input_vcf
    File input_vcf_index
    String chromosome
  }
  
  command {
    set -e
    echo ~{chromosome}
    bcftools view -r chromosome ~{input_vcf}  > ~{sample_basename}.~{chromosome}.vcf.gz
  }

  runtime {
    docker: "dceoy/bcftools"
    maxRetries: 1
    requested_memory_mb_per_core: 1000
    cpu: 4
    runtime_minutes: 60
  }
  output {
    File output_vcf = "~{sample_basename}.part.vcf.gz"
  }
}


task RunVEP {
    input {
      String sample_basename
      File input_vcf
      String annotated_vcf
      String chromosome
    }

    command {
        set -e
        echo ~{chromosome}
        echo ~{annotated_vcf}
        echo "Starting VEP analysis..."
        vep -i ~{input_vcf} \
            -o ~{annotated_vcf} \
            --fork 48 --offline --format vcf --vcf --force_overwrite --compress_output bgzip -v \
            --merged \
            --cache --dir_cache /opt/vep/.vep \
            --plugin AlphaMissense,file=/opt/vep/.vep/Plugins/AlphaMissense/AlphaMissense_hg19.tsv.gz \
            --nearest symbol \
            --shift_hgvs 0 \
            --allele_number \
            --assembly GRCh37 \
            --no_stats 
        tabix -p vcf ~{annotated_vcf}
        echo "Finishing VEP analysis."
        ls -ls ~{annotated_vcf}*
    }

    runtime {
        docker: "alesmaver/vep_grch37"
        requested_memory_mb_per_core: 1000
        cpu: 8
        runtime_minutes: 60
    }

    output {
        File output_vcf = annotated_vcf
        File output_vcf_index = annotated_vcf + ".tbi"
    }
}



task MergeVCFs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String sample_basename
    String output_filename_gz
    String gatk_docker = "broadinstitute/gatk:latest"
    String gatk_path = "/gatk/gatk"
  }
  
  command {
    set -e
    ~{gatk_path} --java-options -Xmx4G  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{output_filename_gz}
  }
  runtime {
    docker: gatk_docker
    maxRetries: 3
    requested_memory_mb_per_core: 1000
    cpu: 8
    runtime_minutes: 120
  }
  output {
    File output_vcfgz = "~{output_filename_gz}"
    File output_vcfgz_index = "~{output_filename_gz}.tbi"
  }
}


task StringToArray {
  input {
    String input_string
    String separator
  }
  command {
    echo '~{input_string}' | tr '~{separator}' \\n | tr -d "[:blank:]" > intervals.list
    echo '~{input_string}' | tr '~{separator}' \\n | tr -d "[:blank:]"
  }
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



