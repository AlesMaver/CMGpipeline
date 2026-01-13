version 1.0

workflow star_salmon {
  input {
    String sample_id
    File input_fq1
    File input_fq2

    String hg19_STAR_genome_index_dir = '/home/cmg/dataSeq/RNAseq_DBs/hg19_STAR'
    String grch38_salmon_index_dir = '/home/cmg/dataSeq/RNAseq_DBs/grch38_salmon_index'

    String star_docker = 'mgibio/star'
    String salmon_docker = 'combinelab/salmon'
  }

  parameter_meta {
    sample_id: "sample name"
    input_fq1: "fq file containing R1"
    input_fq2: "fq file containing R2"
    hg19_STAR_genome_index_dir: "hg19 STAR aligner genome index"
    grch38_salmon_index_dir: "grch38 Salmon transcriptome index"
    star_docker: "star aligner docker"
    salmon_docker: "salmon docker"
  }

  meta {
      author: "Gaber Bergant"
      email: "gaber.bergant@kclj.si"
  }
  
  call RunSTAR {
      input:
        sample_id = sample_id,
        input_fq1 = input_fq1,
        input_fq2 = input_fq2,
        hg19_STAR_genome_index_dir = hg19_STAR_genome_index_dir,
        star_docker = star_docker
    }

  call RunSalmon {
      input:
        sample_id = sample_id,
        input_fq1 = input_fq1,
        input_fq2 = input_fq2,
        grch38_salmon_index_dir = grch38_salmon_index_dir,
        salmon_docker = salmon_docker
    }

  output {
    File? star_bam = RunSTAR.star_bam
  }

}

task RunSTAR {
  input {
    String sample_id
    File input_fq1
    File input_fq2
    File hg19_STAR_genome_index_dir
    String star_docker
  }

  output {
    File star_bam = "~{sample_id}.bam"
  }

  command <<<

    echo "[ RUNNING ] STAR aligner on sample ~{sample_id}"
    STAR \
      --genomeDir ~{hg19_STAR_genome_index_dir} \
      --readFilesIn ~{input_fq1} ~{input_fq2} \
      --readFilesCommand zcat \
      --outFileNamePrefix ~{sample_id} \
      --chimSegmentMin 20 \
      --twopassMode Basic

  >>>
  
  runtime {
    docker: star_docker
    maxRetries: 1
    requested_memory_mb_per_core: 4000
    cpu: 16
    runtime_minutes: 240
  }

}

task RunSalmon {
  input {
    String sample_id
    File input_fq1
    File input_fq2
    File grch38_salmon_index_dir
    String salmon_docker
  }

  output {
    File salmon_quant = "~{sample_id}/quants.sf"
  }

  command <<<

    echo "[ RUNNING ] Salmon quant on sample ~{sample_id}"
    salmon quant \
      -i ~{grch38_salmon_index_dir} \
      -l A \
      -1 ~{input_fq1} \
      -2 ~{input_fq2} \
      --gcBias \
      --validateMappings \
      --no-version-check \
      -o ~{sample_id}

  >>>
  
  runtime {
    docker: star_docker
    maxRetries: 1
    requested_memory_mb_per_core: 4000
    cpu: 16
    runtime_minutes: 240
  }

}
