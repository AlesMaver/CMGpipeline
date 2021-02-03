version 1.0

## Copyright CMG@KIGM, Ales Maver

task calculateBAF {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    String sample_basename

    File reference_fa

    File dbSNPcommon_bed
    File dbSNPcommon_bed_index

    # Runtime parameters
    String docker
  }
  
  command <<<
  set -e
  
  # The perl script converts mpileup to BAF wiggle file for IGV display
  echo '''#!/usr/bin/perl
  use strict;
  use warnings;
  use Getopt::Long qw(:config pass_through no_ignore_case);

  my ($min_reads) = (10);
  GetOptions (
          "min-reads:s" => \$min_reads,
  );

  while (<>) {
          my $line = $_;
          my @columns = split("\t",$line);
          my $chr = $columns[0];
          my $start = $columns[1];
          my $end = $start + 1;
          my $num_reads = $columns[3];
          my $calls = $columns[4];
          my $id = "mpileup_number_" . $.;
          if($num_reads < $min_reads){ # not enough coverage to have good confidence in the call
                  next;
          }
          my $num_ref = 0;
          while ($calls =~ /[,.]/g) { $num_ref++ }
          my $num_var = $num_reads - $num_ref;
          my $varAlleleFreq = ($num_var/$num_reads)*100;
          print("$chr\t$start\t$end\t$varAlleleFreq\n");
          }''' > mpileupToWig.pl

  samtools mpileup -q 15 -Q20 -f ~{reference_fa} -l ~{dbSNPcommon_bed} ~{input_bam} | perl mpileupToWig.pl --min-reads=20 > ~{sample_basename}.BAF.wig  
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 4000
    cpu: 2
    runtime_minutes: 1200
  }
  output {
    File output_BAF = "~{sample_basename}.BAF.wig"
  }
}

task CallROH {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    String sample_basename

    File reference_fa

    File dbSNPcommon_bed
    File dbSNPcommon_bed_index

    # This file is used for determining the locations of common gnomAD SNPs (>1%) for variant calling
    File gnomAD_maf01_vcf
    File gnomAD_maf01_vcf_index

    # This gnomad maf01 tab file contains the frequencies of common gnomAD SNPs (>1%) required for the bcftools roh function
    # The file was generated from the gnomAD vcf file, filtered for common (>1%) SNPs
    # This is performed using the following two commands
    # bcftools query -f'%CHROM\t%POS\t%REF,%ALT\t%INFO/AF\n' gnomad.genomes.r2.1.1.sites.hg19.maf01.vcf.bgz | bgzip -c > gnomad.genomes.r2.1.1.sites.hg19.maf01.tab.gz
    # tabix -s1 -b2 -e2 -f gnomad.genomes.r2.1.1.sites.hg19.maf01.tab.gz # Index as recommended, otherwise bcftools roh won't work
    File gnomAD_maf01_tab
    File gnomAD_maf01_tab_index

    # Runtime parameters
    String docker
  }
  
  command <<<
  set -e
  bcftools mpileup -q 15 -Q20 -f ~{reference_fa} -T ~{dbSNPcommon_bed} ~{input_bam} | bcftools call -m | bcftools view -i 'DP>10 && QUAL>100' -V indels -Oz -o ~{sample_basename}.dbSNP.vcf.gz
  bcftools index -t ~{sample_basename}.dbSNP.vcf.gz
  bcftools roh --AF-file ~{gnomAD_maf01_tab} -G30 -I ~{sample_basename}.dbSNP.vcf.gz ~{sample_basename}.dbSNP.vcf.gz > ~{sample_basename}.bcftoolsROH.output
  cat ~{sample_basename}.bcftoolsROH.output | grep "^[^#]" | grep "^RG" | awk -F'\t' '{if($7>20 && $8>30 && $6>1000000)print $3,$4,$5,$8}' OFS='\t' > ~{sample_basename}.ROHcalls.qual.wig
  cat ~{sample_basename}.bcftoolsROH.output | grep "^[^#]" | grep "^RG" | awk -F'\t' '{if($7>20 && $8>30 && $6>1000000)print $3,$4,$5,$7}' OFS='\t' > ~{sample_basename}.ROHcalls.size.wig
  cat ~{sample_basename}.bcftoolsROH.output | grep "^[^#]" | grep "^ST" | awk -F'\t' '{print $3,$4,$4,$5}' OFS='\t' > ~{sample_basename}.ROHintervals.state.wig
  cat ~{sample_basename}.bcftoolsROH.output | grep "^[^#]" | grep "^ST" | awk -F'\t' '{print $3,$4,$4,$6}' OFS='\t' > ~{sample_basename}.ROHintervals.qual.wig
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 3000
    cpu: 2
    runtime_minutes: 1200
  }
  output {
    File ROH_calls_qual = "~{sample_basename}.ROHcalls.qual.wig"
    File ROH_calls_size = "~{sample_basename}.ROHcalls.size.wig"
    File ROH_intervals_state = "~{sample_basename}.ROHintervals.state.wig"
    File ROH_intervals_qual = "~{sample_basename}.ROHintervals.qual.wig"
    File BAF_vcf = "~{sample_basename}.dbSNP.vcf.gz"
  }
}

task CallPlink {
  input {
    # Command parameters
    File input_vcf
    String sample_basename

    # Runtime parameters
    String docker
  }
  
  command <<<
  set -e
  /usr/local/bin/plink1.9 --vcf ~{input_vcf} --make-bed --out ~{sample_basename} --no-sex --no-parents --no-fid --no-pheno --allow-extra-chr
  /usr/local/bin/plink1.9 -bfile ~{sample_basename} --homozyg --allow-extra-chr --out ~{sample_basename}
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 3000
    cpu: 2
    runtime_minutes: 400
  }
  output {
    File ROHplink_calls = "~{sample_basename}.hom"
    File ROHplink_calls_windows = "~{sample_basename}.hom.summary"
    File ROHplink_calls_summary = "~{sample_basename}.hom.indiv"
  }
}
