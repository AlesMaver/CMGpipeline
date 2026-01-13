version 1.0

task PrepareMaskedGenomeFasta {
  input {
    # Command parameters
    File reference_fixed_fa
    File reference_fixed_fai

    String? targetRegions # FORMAT "chr1:123033-130000;chrX:1-1000"

    # Runtime parameters
    String docker
  }
  
  command <<<
    set -e

    echo "~{targetRegions}"  | tr ';' '\n' | tr ':' '\t' | tr '-' '\t' > targetRegions.bed

    awk -v OFS='\t' {'print $1,$2'} ~{reference_fixed_fai} > ~{reference_fixed_fa}.genome

    bedtools complement  -i targetRegions.bed -g ~{reference_fixed_fa}.genome > targetRegions.complement.bed

    bedtools maskfasta -fi ~{reference_fixed_fa} -bed targetRegions.complement.bed -fo reference.masked.fa
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 4000
    cpu: 4
    runtime_minutes: 800
  }
  output {
    File reference_masked_fa = "reference.masked.fa"
    File targetRegions_bed = "targetRegions.bed"
  }
}


task PrepareMaskedBWAIndex {
  input {
    # Command parameters
    File reference_masked_fa
    String reference_masked_fa_filename = basename(reference_masked_fa)
    # Runtime parameters
    String docker
  }
  
  command {
    set -e

    cp ~{reference_masked_fa} ./

    bwa index ~{reference_masked_fa_filename}
  }
  runtime {
    docker: docker
    requested_memory_mb_per_core: 5000
    cpu: 4
    runtime_minutes: 800
  }
  output {
    File reference_masked_fasta = "~{reference_masked_fa_filename}"
    File reference_masked_sa = "~{reference_masked_fa_filename}.sa"
    File reference_masked_amb = "~{reference_masked_fa_filename}.amb"
    File reference_masked_ann = "~{reference_masked_fa_filename}.ann"
    File reference_masked_bwt = "~{reference_masked_fa_filename}.bwt"
    File reference_masked_pac = "~{reference_masked_fa_filename}.pac"
  }
}






