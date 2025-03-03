version 1.0


# Convert CRAM to BAM file
# Obtained from Broad workflows, here: https://github.com/gatk-workflows/gatk4-germline-snps-indels/blob/master/haplotypecaller-gvcf-gatk4.wdl
task CramToBam {
  input {
    # Command parameters
    File? ref_fasta
    File? ref_fasta_index
    File? ref_dict
    File? input_cram # Declared this as an optional input because the input of workflow is also optional
    String? sample_name

    # Runtime parameters
    String docker
    String samtools_path
  }
  
  command {
    set -e
    set -o pipefail

    ~{samtools_path} view -h -T ~{ref_fasta} ~{input_cram} |
    ~{samtools_path} view -b -o ~{sample_name}.bam -
    ~{samtools_path} index -b ~{sample_name}.bam
    mv ~{sample_name}.bam.bai ~{sample_name}.bai
  }
  runtime {
    docker: docker
    maxRetries: 3
    requested_memory_mb_per_core: 1000
    cpu: 16
    runtime_minutes: 480
 }
  output {
    File output_bam = "~{sample_name}.bam"
    File output_bai = "~{sample_name}.bai"
  }
}

# Convert BAM to CRAM
# Obtained from Broad workflows, here: https://github.com/gatk-workflows/gatk4-genome-processing-pipeline/blob/master/tasks/Utilities.wdl
task ConvertToCram {
  input {
    File input_bam
    File ref_fasta
    File ref_fasta_index
    String sample_basename
  }

  command <<<
    set -e
    set -o pipefail

    samtools view -C -T ~{ref_fasta} ~{input_bam} | \
    tee ~{sample_basename}.cram | \
    md5sum | awk '{print $1}' > ~{sample_basename}.cram.md5

    # Create REF_CACHE. Used when indexing a CRAM
    seq_cache_populate.pl -root ./ref/cache ~{ref_fasta}
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s

    samtools index ~{sample_basename}.cram
    
    mv ~{sample_basename}.cram ~{sample_basename}.hg19.cram
    mv ~{sample_basename}.cram.crai ~{sample_basename}.hg19.cram.crai
    mv ~{sample_basename}.cram.md5 ~{sample_basename}.hg19.cram.md5
  >>>
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    maxRetries: 3
    requested_memory_mb_per_core: 1000
    cpu: 6
    runtime_minutes: 360
  }
  output {
    File output_cram = "~{sample_basename}.hg19.cram"
    File output_cram_index = "~{sample_basename}.hg19.cram.crai"
    File output_cram_md5 = "~{sample_basename}.hg19.cram.md5"
  }
}
