version 1.0

workflow Tryptase {
  input {
    ##File input_bam
    ##File input_bai
    File input_cram
    File input_crai
    File ref_fa
  }

  call Tryptase_call {
      input:
          ##input_bam = input_bam,
          ##input_bai = input_bai,
          input_cram = input_cram,
          input_crai = input_crai,
          ref_fa = ref_fa        
  }
 
  output {
    File tryptase_output_file = Tryptase_call.tryptase_output_file
    File tryptase_output1_file = Tryptase_call.tryptase_output1_file
  }
}

task Tryptase_call {
  input {
    ##File input_bam
    ##File input_bai
    File input_cram
    File input_crai
    File ref_fa
  }
  
  String sample_basename = sub(basename(input_cram), "[\_,\.].*", "" )
  
  command {
    
    # Extract reads mapping to the general tryptase locus from original BAM file
    echo Extract reads mapping to the general tryptase locus from original CRAM file
    samtools view -b -F0xF0C ~{input_cram} -T ~{ref_fa} chr16:1250000-1350000 > ~{sample_basename}.tryptase.bam
    ls -ls ~{sample_basename}.tryptase.bam

    # Convert to interleaved FASTQ
    echo Convert to interleaved FASTQ
    samtools bamshuf -Ou ~{sample_basename}.tryptase.bam temp | samtools bam2fq - > ~{sample_basename}.tryptase.fq
    ls -ls ~{sample_basename}.tryptase.fq

    # Map reads to consensus sequence, retaininly only mapped reads
    echo Map reads to consensus sequence, retaininly only mapped reads 
    bwa mem -pM -R '@RG\tID:1\tSM:~{sample_basename}' /usr/working/Tryptase/CONS.fa ~{sample_basename}.tryptase.fq | samtools view -F0xF0C -S -h - > ~{sample_basename}.cons.sam
    ls -ls ~{sample_basename}.cons.sam
    
    # Cluster mapped reads into distinct haplotypes
    echo Cluster mapped reads into distinct haplotypes - parseHaplotypes - stay patient...
    cat ~{sample_basename}.cons.sam \
      | awk '$10!~/N/' \
      | python /usr/working/Tryptase/parseHaplotypes.py \
      > ~{sample_basename}.fa
     
    ls -ls ~{sample_basename}.fa

    cat ~{sample_basename}.fa \
      | sed 's/ .*//' \
      | sed 's/^X*//' \
      | sed 's/X*$//' \
      | tr X N \
      | bwa mem -M /usr/working/Tryptase/TPS.fa - \
      | awk '$1!~/_[0-9]_/' \
      > ~{sample_basename}.consx.sam
      
    ls -ls ~{sample_basename}.consx.sam


    cat ~{sample_basename}.consx.sam |grep -v ^@ \
      | cut -f1,3,4,6,12-13|tr _ "\t" \
      | awk '$3>5' \
      | sort -k4,4 -k3nr \
      > ~{sample_basename}.tryptase.output.txt
      
    ls -ls ~{sample_basename}.tryptase.output.txt
    
    # --------------------------------------------
    # vse skupaj ponovimo še s CONS1.fa in TPS1.fa
    # --------------------------------------------
    echo Ponovimo vse skupaj še s CONS1.fa in TPS1.fa
    
    # Map reads to consensus sequence, retaininly only mapped reads 
    echo Map reads to consensus sequence, retaininly only mapped reads 
    bwa mem -pM -R '@RG\tID:1\tSM:~{sample_basename}' /usr/working/Tryptase/CONS1.fa ~{sample_basename}.tryptase.fq | samtools view -F0xF0C -S -h - > ~{sample_basename}.cons1.sam
    ls -ls ~{sample_basename}.cons1.sam
    
    # Cluster mapped reads into distinct haplotypes
    echo Cluster mapped reads into distinct haplotypes  - parseHaplotypes - stay patient... 
    cat ~{sample_basename}.cons1.sam \
      | awk '$10!~/N/' \
      | python /usr/working/Tryptase/parseHaplotypes.py \
      > ~{sample_basename}.1.fa
     
    ls -ls ~{sample_basename}.1.fa

    cat ~{sample_basename}.1.fa \
      | sed 's/ .*//' \
      | sed 's/^X*//' \
      | sed 's/X*$//' \
      | tr X N \
      | bwa mem -M /usr/working/Tryptase/TPS1.fa - \
      | awk '$1!~/_[0-9]_/' \
      > ~{sample_basename}.consx1.sam
      
    ls -ls ~{sample_basename}.consx1.sam


    cat ~{sample_basename}.consx1.sam |grep -v ^@ \
      | cut -f1,3,4,6,12-13|tr _ "\t" \
      | awk '$3>5' \
      | sort -k4,4 -k3nr \
      > ~{sample_basename}.tryptase.output.1.txt
      
    ls -ls ~{sample_basename}.tryptase.output.1.txt
    
    echo END.
        
  }

  runtime {
    docker: "alesmaver/bwa_samtools_picard_tryptase"
    requested_memory_mb_per_core: 4000
    cpu: 2
    runtime_minutes: 600
    continueOnReturnCode: true
  }
  output {
    File tryptase_output_file = "~{sample_basename}.tryptase.output.txt"
    File tryptase_output1_file = "~{sample_basename}.tryptase.output.1.txt"
  }
}
