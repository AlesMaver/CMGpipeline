version 1.0

workflow DELLY {
  input {
    Array[File] input_bams
    Array[File] input_bais

    File reference_fasta
  }

  Array[Pair[File,File]] bam_bai_pairs = zip(input_bams, input_bais)

  scatter (bam_bai_pair in bam_bai_pairs) {
        File input_bam = bam_bai_pair.left
        File input_bai = bam_bai_pair.right

        call DELLY_call {
            input:
                input_bam = input_bam,
                input_bai = input_bai,
                reference_fasta = reference_fasta
        }
  }

  call DELLY_merge {
    input:
        input_bcfs = DELLY_call.bcf_file
  }

  scatter (bam_bai_pair in bam_bai_pairs) {
        File input_bam2 = bam_bai_pair.left
        File input_bai2 = bam_bai_pair.right

        call DELLY_genotype {
            input:
                input_bam = input_bam2,
                input_bai = input_bai2,
                reference_fasta = reference_fasta,
                sites_bcf_file = DELLY_merge.sites_bcf_file
        }
  }

  call DELLY_merge_genotype {
    input:
        input_geno_bcfs = DELLY_genotype.geno_bcf_file,
        input_geno_bcf_indices = DELLY_genotype.geno_bcf_index_file
  }

call DELLY_filter {
    input:
        input_bcf = DELLY_merge_genotype.merged_geno_file
}

  output {
    File? delly_bcf = DELLY_filter.filtered_bcf
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
    docker: "alesmaver/delly"
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
    docker: "alesmaver/delly"
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
    docker: "alesmaver/delly"
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
    Array[File] input_geno_bcfs
    Array[File] input_geno_bcf_indices
  }

  command {
    bcftools merge -m id -O b -o merged.bcf ~{sep=' ' input_geno_bcfs}
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
  }

  command {
    bcftools index ~{input_bcf}
    delly filter -f germline -o germline.bcf ~{input_bcf}
    bcftools view -Ov -o germline.vcf germline.bcf
  }

  runtime {
    docker: "alesmaver/delly"
    requested_memory_mb_per_core: 4000
    cpu: 2
    runtime_minutes: 200
  }
  output {
    File filtered_bcf = "merged.bcf" 
  }
}   