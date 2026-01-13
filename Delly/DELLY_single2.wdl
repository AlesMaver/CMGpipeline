version 1.0

import "../manta/manta_workflow.wdl" as manta

workflow DELLY {
  input {
    File input_bam
    File input_bai
    File population_bcf
    File population_bcf_index
    File reference_fasta
  }

  String sample_basename = sub(basename(input_bam), "[\_,\.].*", "" )

  call DELLY_call {
      input:
          input_bam = input_bam,
          input_bai = input_bai,
          reference_fasta = reference_fasta
  }

  Array[File] input_bcfs = [DELLY_call.bcf_file, population_bcf]

  call DELLY_merge {
    input:
        input_bcfs = input_bcfs
  }

  call DELLY_genotype {
      input:
          input_bam = input_bam,
          input_bai = input_bai,
          reference_fasta = reference_fasta,
          sites_bcf_file = DELLY_merge.sites_bcf_file
  }

  call DELLY_merge_genotype {
    input:
        input_geno_bcf = DELLY_genotype.geno_bcf_file,
        input_geno_bcf_index = DELLY_genotype.geno_bcf_index_file,
        input_population_bcf = population_bcf,
        input_population_bcf_index = population_bcf_index,
        sample_basename = sample_basename
  }

    call DELLY_filter {
        input:
            input_bcf = DELLY_merge_genotype.merged_geno_file,
            sample_basename = sample_basename
    }

    call manta.annotSV as annotSV {
        input:
            genome_build = "GRCh37",
            input_vcf = DELLY_filter.filtered_vcf,
            output_tsv_name = sample_basename + ".DELLY.AnnotSV.tsv"
    }

  output {
    File sample_bcfs = DELLY_genotype.geno_bcf_file
    File sample_bcf_indices = DELLY_genotype.geno_bcf_index_file
    File unfiltered_population_bcf = DELLY_merge_genotype.merged_geno_file
    File filtered_population_bcf = DELLY_filter.filtered_bcf
  }
}



task DELLY_call {
  input {
    File input_bam
    File input_bai
    File reference_fasta
  }
  
  String sample_basename = sub(basename(input_bam), "[\_,\.].*", "" )

  command {
    wget https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg19.excl.tsv
    delly call -g ~{reference_fasta} -o ~{sample_basename}.bcf -x human.hg19.excl.tsv ~{input_bam}
  }

  runtime {
    docker: "alesmaver/delly2"
    requested_memory_mb_per_core: 4000
    cpu: 2
    runtime_minutes: 200
  }
  output {
    File bcf_file = "~{sample_basename}.bcf" 
  }
}

task DELLY_merge {
  input {
    Array[File] input_bcfs
  }

  command {
    delly merge -o sites.bcf ~{sep=' ' input_bcfs}
  }

  runtime {
    docker: "alesmaver/delly2"
    requested_memory_mb_per_core: 4000
    cpu: 2
    runtime_minutes: 200
  }
  output {
    File sites_bcf_file = "sites.bcf" 
  }
}

task DELLY_genotype {
  input {
    File input_bam
    File input_bai
    File sites_bcf_file
    File reference_fasta
  }
  
  String sample_basename = sub(basename(input_bam), "[\_,\.].*", "" )

  command {
    wget https://raw.githubusercontent.com/dellytools/delly/main/excludeTemplates/human.hg19.excl.tsv
    delly call -g ~{reference_fasta} -v ~{sites_bcf_file} -o ~{sample_basename}.geno.bcf -x human.hg19.excl.tsv ~{input_bam}
  }

  runtime {
    docker: "alesmaver/delly2"
    requested_memory_mb_per_core: 4000
    cpu: 2
    runtime_minutes: 200
  }
  output {
    File geno_bcf_file = "~{sample_basename}.geno.bcf" 
    File geno_bcf_index_file = "~{sample_basename}.geno.bcf.csi" 
  }
}

task DELLY_merge_genotype {
  input {
    File input_geno_bcf
    File input_geno_bcf_index
    File input_population_bcf
    File input_population_bcf_index
    String sample_basename
  }

  command {
    bcftools query -l ~{input_population_bcf} > samples.txt
    if grep -q '~{sample_basename}' samples.txt
    then
        cp ~{input_population_bcf} merged.bcf
    else
        bcftools merge ~{input_population_bcf} ~{input_geno_bcf} | bcftools +fill-tags | bcftools view -Ob -o merged.bcf
    fi
    
    #bcftools merge -m id -O b -o merged.bcf ~{sep=' ' input_geno_bcf}
  }

  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    requested_memory_mb_per_core: 4000
    cpu: 2
    runtime_minutes: 200
  }
  output {
    File merged_geno_file = "merged.bcf" 
  }
}

task DELLY_filter {
  input {
    File input_bcf
    String sample_basename
  }

  command {
    bcftools index ~{input_bcf}
    bcftools view -i 'INFO/AC<4 || INFO/AC="."' ~{input_bcf} | bcftools view --no-update -s ~{sample_basename} | bcftools view --no-update -e 'GT="0/0" || GT="./."' | bcftools view -Ob -o ~{sample_basename}.bcf
    #delly filter -f germline -o population.bcf ~{input_bcf}
    bcftools view -Ov -o ~{sample_basename}.vcf ~{sample_basename}.bcf
  }

  runtime {
    docker: "alesmaver/delly2"
    requested_memory_mb_per_core: 4000
    cpu: 2
    runtime_minutes: 200
  }
  output {
    File filtered_bcf = "~{sample_basename}.bcf" 
    File filtered_vcf = "~{sample_basename}.vcf" 
  }
}   