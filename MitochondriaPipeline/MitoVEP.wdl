version 1.0

workflow MitoVEP {
  input {
    String sample_basename
    File input_vcf
    File reference_fasta
    File gnomad_mito_sites_vcf
    File gnomad_mito_sites_vcf_index
    String assembly = "GRCh38"
    Int fork = 48
    String output_filename = sample_basename + ".hg38.mitorsaw.VEP.vcf.gz"
    String output_filename_most_severe = sample_basename + ".mitorsaw.VEP.mostSevere.vcf.gz"
    String docker = "alesmaver/vep"
  }

  call RunMitoVEP {
    input:
      input_vcf = input_vcf,
      output_vcf = output_filename,
      output_vcf_most_severe = output_filename_most_severe,
      reference_fasta = reference_fasta,
      gnomad_mito_sites_vcf = gnomad_mito_sites_vcf,
      gnomad_mito_sites_vcf_index = gnomad_mito_sites_vcf_index,
      assembly = assembly,
      fork = fork,
      docker = docker
  }

  output {
    File output_vcf = RunMitoVEP.output_vcf
    File output_vcf_index = RunMitoVEP.output_vcf_index
    File output_vcf_most_severe = RunMitoVEP.output_vcf_most_severe
    File output_vcf_most_severe_index = RunMitoVEP.output_vcf_most_severe_index
  }
}

task RunMitoVEP {
  input {
    File input_vcf
    String output_vcf
    String output_vcf_most_severe
    File reference_fasta
    File gnomad_mito_sites_vcf
    File gnomad_mito_sites_vcf_index
    String assembly
    Int fork
    String docker
  }

  command <<<
    set -e

    vep -i ~{input_vcf} \
      -o ~{output_vcf} \
      --fork ~{fork} \
      --offline \
      --format vcf \
      --vcf \
      --force_overwrite \
      --compress_output bgzip \
      --cache \
      --dir_cache /opt/vep/.vep \
      --assembly ~{assembly} \
      --nearest symbol \
      --shift_hgvs 0 \
      --allele_number \
      --no_stats \
      --symbol \
      --biotype \
      --canonical \
      --mane \
      --tsl \
      --appris \
      --ccds \
      --uniprot \
      --domains \
      --protein \
      --numbers \
      --variant_class \
      --gene_phenotype \
      --regulatory \
      --pubmed \
      --af \
      --af_gnomad \
      --max_af \
      --clin_sig_allele 1 \
      --hgvs \
      --hgvsg \
      --fasta ~{reference_fasta} \
      --custom ~{gnomad_mito_sites_vcf},gnomAD_MT,vcf,exact,0,AN,AC_hom,AC_het

    tabix -p vcf ~{output_vcf}


    vep -i ~{input_vcf} \
      -o ~{output_vcf_most_severe} \
      --fork ~{fork} \
      --offline \
      --format vcf \
      --vcf \
      --force_overwrite \
      --compress_output bgzip \
      --cache \
      --dir_cache /opt/vep/.vep \
      --assembly ~{assembly} \
      --nearest symbol \
      --shift_hgvs 0 \
      --allele_number \
      --no_stats \
      --symbol \
      --biotype \
      --canonical \
      --mane \
      --tsl \
      --appris \
      --ccds \
      --uniprot \
      --domains \
      --protein \
      --numbers \
      --variant_class \
      --gene_phenotype \
      --regulatory \
      --pubmed \
      --af \
      --af_gnomad \
      --max_af \
      --clin_sig_allele 1 \
      --pick \
      --hgvs \
      --hgvsg \
      --fasta ~{reference_fasta} \
      --custom ~{gnomad_mito_sites_vcf},gnomAD_MT,vcf,exact,0,AN,AC_hom,AC_het

    tabix -p vcf ~{output_vcf_most_severe}
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 1000
    cpu: fork
    runtime_minutes: 240
    bootDiskSizeGb: "150"
  }

  output {
    File output_vcf = output_vcf
    File output_vcf_index = output_vcf + ".tbi"
    File output_vcf_most_severe = output_vcf_most_severe
    File output_vcf_most_severe_index = output_vcf_most_severe + ".tbi"
  }
}
