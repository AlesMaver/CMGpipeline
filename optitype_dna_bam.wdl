version 1.0

task OptitypeDnafromBam {
  input {
    String optitype_name = "optitype"
    File input_bam
  }

  Int space_needed_gb = 10 + round(5*size([input_bam], "GB"))
  runtime {
    memory: "64GB"
    docker: "mgibio/immuno_tools-cwl:1.0.1"
    disks: "local-disk ~{space_needed_gb} SSD"
    bootDiskSizeGb: 3*space_needed_gb
    runtime_minutes: 2400
  }

  command <<<
    set -e -o pipefail
    dnaref="/ref_data/optitype_ref/hla_reference_dna.fasta"    # Optitype DNA reference file
    name=~{optitype_name}
    echo $name

    echo Sorting bam
    sambamba sort  -n -t 4 -m 8G -o $name.qsorted.bam  ~{input_bam}                ## 4-threaded replacement sorting with sambamba:
    
    echo Running bedtools bamtofastq
    /usr/bin/bedtools bamtofastq -fq  $name.q.fwd.fastq -fq2  $name.q.rev.fastq -i $name.qsorted.bam 2>/dev/null;
    rm $name.qsorted.bam
    
    echo Aligning forward reads to reference HLA locus sequence
    /usr/local/bin/bwa mem -t 4 $dnaref $name.q.fwd.fastq > $name.aln.fwd.sam      # use bwa mem, store output IN TEMP, and skip samse step
    rm -f $name.q.fwd.fastq
    
    echo Aligning reverse reads to reference HLA locus sequence
    /usr/local/bin/bwa mem -t 4 $dnaref $name.q.rev.fastq > $name.aln.rev.sam      # use bwa mem, store output IN TEMP, and skip samse step
    rm -f $name.q.rev.fastq

    echo Select only the mapped reads from the sam files:
    /opt/samtools/bin/samtools view -S -F 4 $name.aln.fwd.sam > $name.aln.map.fwd.sam
    /opt/samtools/bin/samtools view -S -F 4 $name.aln.rev.sam > $name.aln.map.rev.sam
    rm -f  $name.aln.fwd.sam
    rm -f  $name.aln.rev.sam

    echo Convert sam files to fastq files, also stored in temp dir
    cat $name.aln.map.fwd.sam | grep -v ^@ | /usr/bin/awk '{print "@"$1"\n"$10"\n+\n"$11}' > $name.hla.fwd.fastq
    cat $name.aln.map.rev.sam | grep -v ^@ | /usr/bin/awk '{print "@"$1"\n"$10"\n+\n"$11}' > $name.hla.rev.fastq
    rm -f $name.aln.map.fwd.sam
    rm -f $name.aln.map.rev.sam

    echo step 5: run Optitype
    # run optitype 
    /usr/bin/python /usr/local/bin/OptiType/OptiTypePipeline.py -i $name.hla.fwd.fastq $name.hla.rev.fastq --dna -v -p $name -o .

    mv ./${name}_result.tsv ./${name}.optitype_result.tsv
    mv ./${name}_coverage_plot.pdf ./${name}.optitype_coverage_plot.pdf
  >>>

  output {
    File optitype_tsv = optitype_name + ".optitype_result.tsv"
    File optitype_plot = optitype_name + ".optitype_coverage_plot.pdf"
  }
}

workflow OptitypeDna_fromBam {
  input {
    String? optitype_name
    File input_bam
  }
  call OptitypeDnafromBam {
    input:
    optitype_name=optitype_name,
    input_bam = input_bam,
  }
}
