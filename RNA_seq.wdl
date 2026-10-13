version 1.0

## orig: import "https://raw.githubusercontent.com/AlesMaver/bioWDL_RNAseq/develop/align-hisat2.wdl" as AlignHisat2
import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/RNA_seq/align-hisat2.wdl" as AlignHisat2
import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/GenerateCRAM.wdl" as GenerateCRAM

workflow RNA_seq {
  input {
        String sample_basename

        # data for AlignHisat2:
        ### Array[FastqPair]+ inputReads 
        # Provide either a FASTQ pair or a CRAM; when input_cram is set the
        # reads are extracted from it and the FASTQ inputs are ignored.
        File? input_fq1
        File? input_fq2
        File? input_cram

        String outputDir = "."
        String library
        Array[String] readgroups
        Array[File]+ indexFiles

        # references for CRAM generation:
        File reference_fa
        File reference_fai
  }

  #FastqPair inputReads =  { 
  #            "R1": input_fq1,
  #            "R2": input_fq2 
  #          }

  if (defined(input_cram)) {
    call CramToFastq {
      input:
        input_cram = select_first([input_cram]),
        reference_fa = reference_fa,
        reference_fai = reference_fai,
        sample_basename = sample_basename
    }
  }

  File fq1 = select_first([CramToFastq.fastq_r1, input_fq1])
  File fq2 = select_first([CramToFastq.fastq_r2, input_fq2])

  Array[FastqPair] inputReads =  [
    { "R1": fq1, "R2": fq2 }
  ]


  call AlignHisat2.AlignHisat2 as AlignHisat2 {
      input:
        inputReads = inputReads, 
        outputDir = outputDir,
        sample = sample_basename,
        library = library,
        readgroups = readgroups,
        ### platform = "illumina",
        indexFiles = indexFiles
  }

  call GenerateCRAM.GenerateCRAM as GenerateCRAM {
      input:       
        input_bam = AlignHisat2.bamFile["file"],
        reference_fa = reference_fa,
        reference_fai = reference_fai,
        sample_basename = sample_basename
  }

  output {
    #File output_bam = AlignHisat2.bamFile["file"] 
    #File output_bam_index = AlignHisat2.bamFile["index"]
    
    File output_cram = GenerateCRAM.output_cram 
    File output_cram_index = GenerateCRAM.output_cram_index 
  }
}

# Extract paired reads (primary alignments only, including unmapped) from a CRAM
task CramToFastq {
  input {
    File input_cram
    File reference_fa
    File reference_fai
    String sample_basename
    Int threads = 4
    String docker = "quay.io/biocontainers/samtools:1.21--h50ea8bc_0"
  }

  command <<<
    set -euo pipefail
    samtools collate -@ ~{threads} -u -O --reference ~{reference_fa} ~{input_cram} collate_tmp \
      | samtools fastq -@ ~{threads} -F 0x900 -n -0 /dev/null -s /dev/null \
          -1 ~{sample_basename}_R1.fastq.gz -2 ~{sample_basename}_R2.fastq.gz -
  >>>

  output {
    File fastq_r1 = "~{sample_basename}_R1.fastq.gz"
    File fastq_r2 = "~{sample_basename}_R2.fastq.gz"
  }

  runtime {
    docker: docker
    maxRetries: 3
    requested_memory_mb_per_core: 1000
    cpu: threads
    runtime_minutes: 300
  }
}

struct FastqPair {
    File R1
    String? R1_md5
    File? R2
    String? R2_md5
}
