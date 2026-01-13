version 1.0

task OptitypeDnafromFastq {
  input {
    String optitype_name = "optitype"
    File input_fq1
    File input_fq2
  }

  Int space_needed_gb = 10 + round(5*size([input_fq1, input_fq2], "GB"))
  runtime {
    memory: "64GB"
    docker: "mgibio/immuno_tools-cwl:1.0.1"
    disks: "local-disk ~{space_needed_gb} SSD"
    bootDiskSizeGb: 3*space_needed_gb
  }

  command <<<
    set -e -o pipefail
    dnaref="/ref_data/optitype_ref/hla_reference_dna.fasta"    # Optitype DNA reference file
    name=~{optitype_name}
    echo $name

    echo Aligning forward reads to reference HLA locus sequence
    /usr/local/bin/bwa mem -t 4 $dnaref ~{input_fq1} > $name.aln.fwd.sam      # use bwa mem, store output IN TEMP, and skip samse step
    echo Aligning reverse reads to reference HLA locus sequence
    /usr/local/bin/bwa mem -t 4 $dnaref ~{input_fq2} > $name.aln.rev.sam      # use bwa mem, store output IN TEMP, and skip samse step

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

workflow OptitypeDnafromFastq {
  input {
    String? optitype_name
    File input_fq1
    File input_fq2
  }
  call OptitypeDnafromFastq {
    input:
    optitype_name=optitype_name,
    input_fq1 = input_fq1,
    input_fq2 = input_fq2,
  }
}
