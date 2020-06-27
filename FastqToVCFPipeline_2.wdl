version 1.0
## Copyright CMG@KIGM, Ales Maver

import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/AnnotationPipeline.wdl" as Annotation
import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/Conifer.wdl" as Conifer

# WORKFLOW DEFINITION 
workflow FastqToVCF {
  input {
    File input_fq1
    File input_fq2
    
    File illuminaAdapters

    File chromosome_list

    Int bwa_threads
    Int threads
    Int split_reads_num

    File reference_fixed_fa
    File reference_fixed_fai
    File reference_fixed_amb
    File reference_fixed_ann
    File reference_fixed_bwt
    File reference_fixed_pac
    File reference_fixed_sa

    File reference_fa
    File reference_fai
    File reference_dict
    
    File gnomAD_vcf
    File gnomAD_vcf_index

    File gnomADexomes_vcf
    File gnomADexomes_vcf_index

    File SLOpopulation_vcf
    File SLOpopulation_vcf_index

    File ClinVar_vcf
    File ClinVar_vcf_index

    File SpliceAI
    File SpliceAI_index

    File dbscSNV
    File dbscSNV_index

    File HPO
    File HPO_index
    File OMIM
    File OMIM_index
    File gnomadConstraints
    File gnomadConstraints_index
    File CGD
    File CGD_index
    File bcftools_annotation_header

    File dbNSFP
    File dbNSFP_index

    File dbsnp_vcf
    File dbsnp_vcf_index
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices

    String enrichment
    File enrichment_bed

    Array[File] input_reference_rpkms 
    Int CONIFER_svd
    Float CONIFER_threshold

    String cutadapt_docker = "kfdrc/cutadapt:latest"
    String gatk_docker = "broadinstitute/gatk:latest"
    String gatk_path = "/gatk/gatk"
    String vcfanno_docker = "clinicalgenomics/vcfanno:0.3.2"
    String bcftools_docker = "biocontainers/bcftools:v1.9-1-deb_cv1"
    String SnpEff_docker = "alesmaver/snpeff_v43:latest"
  }  

  String sample_basename = sub(basename(input_fq1), "[\_,\.].*", "" )

  Array[String] chromosomes = read_lines(chromosome_list)


  call CutAdapters as CutAdapters_fq1 {
    input:
      input_fq=input_fq1,
      sample_basename=sample_basename,
      illuminaAdapters=illuminaAdapters,

      # Runtime 
      docker = cutadapt_docker
  }
  
  call CutAdapters as CutAdapters_fq2 {
    input:
      input_fq=input_fq2,
      sample_basename=sample_basename,
      illuminaAdapters=illuminaAdapters,

      # Runtime 
      docker = cutadapt_docker
  }

  call PairedFastQsToUnmappedBAM {
      input:
        sample_name = sample_basename,
        fastq_1 = CutAdapters_fq1.output_fq_trimmed,
        fastq_2 = CutAdapters_fq2.output_fq_trimmed,
        readgroup_name = sample_basename,
        library_name = sample_basename,
        platform_unit = "NextSeq550",
        run_date = "2020-01-01",
        platform_name = "illumina",
        sequencing_center = "KIGM",
        gatk_path = gatk_path,
        docker = gatk_docker
    }

  call SamSplitter {
    input :
      input_bam = PairedFastQsToUnmappedBAM.output_unmapped_bam,
      n_reads = split_reads_num,
      preemptible_tries = 3,
      compression_level = 2
  }

  scatter(unmapped_bam in SamSplitter.split_bams) {
    call Align {
      input:
       input_bam=unmapped_bam,
       sample_basename=sample_basename,

       reference_fixed_fa=reference_fixed_fa,
       reference_fixed_fai=reference_fixed_fai,
       reference_fixed_amb=reference_fixed_amb,
       reference_fixed_ann=reference_fixed_ann,
       reference_fixed_bwt=reference_fixed_bwt,
       reference_fixed_pac=reference_fixed_pac,
       reference_fixed_sa=reference_fixed_sa,

       reference_fa=reference_fa,
       reference_dict=reference_dict,

       picard_path="/usr/picard/picard.jar",

       threads=bwa_threads,

       docker="alesmaver/bwa_samtools_picard"
    }
  }

  call GatherUnsortedBamFiles {
    input:
      input_bams = Align.output_bam,
      output_bam_basename = sample_basename,
      preemptible_tries = 3,
      compression_level = 2
  }

  # Sort aggregated+deduped BAM file and fix tags
  call SortSam {
    input:
      input_bam = GatherUnsortedBamFiles.output_bam,
      output_bam_basename = sample_basename + ".marked",
      compression_level = 2,
      preemptible_tries = 3
  }

  Float agg_bam_size = size(SortSam.output_bam, "GiB")

  # Create list of sequences for scatter-gather parallelization
  call CreateSequenceGroupingTSV {
    input:
      ref_dict = reference_dict,
      preemptible_tries = 3
  }

  # We need disk to localize the sharded input and output due to the scatter for BQSR.
  # If we take the number we are scattering by and reduce by 3 we will have enough disk space
  # to account for the fact that the data is not split evenly.
  Int num_of_bqsr_scatters = length(CreateSequenceGroupingTSV.sequence_grouping)
  Int potential_bqsr_divisor = num_of_bqsr_scatters - 10
  Int bqsr_divisor = if potential_bqsr_divisor > 1 then potential_bqsr_divisor else 1

  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    # Generate the recalibration model by interval
    call BaseRecalibrator {
      input:
        input_bam = SortSam.output_bam,
        input_bam_index = SortSam.output_bam_index,
        recalibration_report_filename = sample_basename + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbsnp_vcf = dbsnp_vcf,
        dbsnp_vcf_index = dbsnp_vcf_index,
        known_indels_sites_vcfs = known_indels_sites_vcfs,
        known_indels_sites_indices = known_indels_sites_indices,
        ref_dict = reference_dict,
        ref_fasta = reference_fa,
        ref_fasta_index = reference_fai,
        bqsr_scatter = bqsr_divisor,
        preemptible_tries = 3
    }
  }

  # Merge the recalibration reports resulting from by-interval recalibration
  # The reports are always the same size
  call GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = sample_basename + ".recal_data.csv",
      preemptible_tries = 3
  }

  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped) {
    # Apply the recalibration model by interval
    call ApplyBQSR {
      input:
        input_bam = SortSam.output_bam,
        input_bam_index = SortSam.output_bam_index,
        output_bam_basename = sample_basename,
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = reference_dict,
        ref_fasta = reference_fa,
        ref_fasta_index = reference_fai,
        bqsr_scatter = bqsr_divisor,
        compression_level = 2,
        preemptible_tries = 3,
        bin_base_qualities = true,
        somatic = false
    }
  }

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherSortedBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = sample_basename,
      total_input_size = agg_bam_size,
      compression_level = 2,
      preemptible_tries = 3
  }

  scatter (chromosome in chromosomes) {
    call HaplotypeCaller {
      input:
        input_bam = GatherSortedBamFiles.output_bam,
        input_bai = GatherSortedBamFiles.output_bam_index,
        sample_basename=sample_basename,

        reference_fa=reference_fa,
        reference_fai=reference_fai,
        reference_dict=reference_dict,

        chromosome=chromosome,

        threads=threads,

        gatk_path=gatk_path,
        docker=gatk_docker
    }
  }

  # Merge per-interval GVCFs
  call MergeVCFs {
    input:
      input_vcfs = HaplotypeCaller.output_vcf,
      input_vcfs_indexes = HaplotypeCaller.output_vcf_index,
      sample_basename = sample_basename,
      docker = gatk_docker,
      gatk_path = gatk_path
  }

  call SplitSNPindel {
  input:
    input_vcf = MergeVCFs.output_vcf,
    input_vcf_index = MergeVCFs.output_vcf_index,
    sample_basename=sample_basename,

    reference_fa=reference_fa,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    threads=threads,

    gatk_path=gatk_path,
    docker=gatk_docker
  }

  call VariantFiltrationSNP {
  input:
    input_vcf = SplitSNPindel.output_vcf_SNP,
    sample_basename=sample_basename,

    reference_fa=reference_fa,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    threads=threads,

    gatk_path=gatk_path,
    docker=gatk_docker
  }

  call VariantFiltrationINDEL {
  input:
    input_vcf = SplitSNPindel.output_vcf_INDEL,
    sample_basename=sample_basename,

    reference_fa=reference_fa,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    threads=threads,

    gatk_path=gatk_path,
    docker=gatk_docker
  }

  call CombineSNPindel {
  input:
    input_snp_vcf = VariantFiltrationSNP.output_vcf,
    input_indel_vcf = VariantFiltrationINDEL.output_vcf,
    sample_basename=sample_basename,

    reference_fa=reference_fa,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    threads=threads,

    picard_path="/usr/picard/picard.jar",
    docker="alesmaver/bwa_samtools_picard"
  }

  call SelectFinalVariants {
  input:
    input_vcf = CombineSNPindel.output_vcf,
    sample_basename=sample_basename,

    reference_fa=reference_fa,
    reference_fai=reference_fai,
    reference_dict=reference_dict,

    threads=threads,

    gatk_path=gatk_path,
    docker=gatk_docker
  }

  call Annotation.AnnotateVCF as AnnotateVCF{
    input:
      input_vcf = SelectFinalVariants.output_vcf,
      chromosome_list = chromosome_list,
      
      gnomAD_vcf = gnomAD_vcf,
      gnomAD_vcf_index = gnomAD_vcf_index,

      gnomADexomes_vcf = gnomADexomes_vcf,
      gnomADexomes_vcf_index = gnomADexomes_vcf_index,

      SLOpopulation_vcf = SLOpopulation_vcf,
      SLOpopulation_vcf_index = SLOpopulation_vcf_index,

      ClinVar_vcf = ClinVar_vcf,
      ClinVar_vcf_index = ClinVar_vcf_index,

      SpliceAI = SpliceAI,
      SpliceAI_index = SpliceAI_index,

      dbscSNV = dbscSNV,
      dbscSNV_index = dbscSNV_index,

      HPO = HPO,
      HPO_index = HPO_index,
      OMIM = OMIM,
      OMIM_index = OMIM_index,
      gnomadConstraints = gnomadConstraints,
      gnomadConstraints_index = gnomadConstraints_index,
      CGD = CGD,
      CGD_index = CGD_index,
      bcftools_annotation_header = bcftools_annotation_header,

      fasta_reference = reference_fa,
      fasta_reference_index = reference_fai,
      fasta_reference_dict = reference_dict,

      dbNSFP = dbNSFP,
      dbNSFP_index = dbNSFP_index,

      bcftools_docker = bcftools_docker,
      SnpEff_docker = SnpEff_docker,
      gatk_docker = gatk_docker,
      gatk_path = gatk_path,
      vcfanno_docker = vcfanno_docker
  }

  call Conifer.Conifer as Conifer{
  input:
    input_bam = SortSam.output_bam,
    input_bam_index = SortSam.output_bam_index

    Array[File] input_reference_rpkms 
    Int CONIFER_svd
    Float CONIFER_threshold

    String enrichment
    File enrichment_bed
  }


  output {
    File output_bam = SortSam.output_bam
    File output_bam_index = SortSam.output_bam_index

    File output_vcf_raw = MergeVCFs.output_vcf
    File output_vcf_raw_index = MergeVCFs.output_vcf_index

    File output_vcf = SelectFinalVariants.output_vcf
    File output_vcf_index = SelectFinalVariants.output_vcf_index

    File output_annotated_vcf = AnnotateVCF.output_vcf
    File output_annotated_vcf_index = AnnotateVCF.output_vcf_index
    File output_table = AnnotateVCF.output_table

    File output_conifer_calls = Conifer.output_conifer_calls
    Array[File] output_plotcalls = Conifer.output_plotcalls
  }
}





##################
# TASK DEFINITIONS
##################

task CutAdapters {
  input {
    # Command parameters
    File input_fq
    String sample_basename
    
    File illuminaAdapters

    # Runtime parameters
    String docker
  }
  
  command {
  set -e
     cutadapt -j 20 -a file:~{illuminaAdapters} --mask-adapter -o ~{sample_basename}.trimmed.fq.gz ~{input_fq}
  }
  runtime {
    docker: docker
  }
  output {
    File output_fq_trimmed = "~{sample_basename}.trimmed.fq.gz"
  }
}

task Align {
  input {
  File input_bam

  File reference_fixed_fa
  File reference_fixed_fai
  File reference_fixed_amb
  File reference_fixed_ann
  File reference_fixed_bwt
  File reference_fixed_pac
  File reference_fixed_sa

  File reference_fa
  File reference_dict

  String sample_basename

  String picard_path="/usr/picard/picard.jar"

  # Runtime params
  Int threads
  String docker
  }
  
  command {
    java -Xms1000m -Xmx1000m -jar ~{picard_path} SamToFastq INPUT=~{input_bam} FASTQ=/dev/stdout INTERLEAVE=true | bwa mem -p -t ~{threads} ~{reference_fixed_fa} /dev/stdin | samtools sort -@ ~{threads} - | java -jar ~{picard_path} AddOrReplaceReadGroups I=/dev/stdin O=/dev/stdout RGID=4 RGLB=~{sample_basename} RGPL=illumina RGPU=unit1 RGSM=~{sample_basename} | java -jar ~{picard_path} ReorderSam I=/dev/stdin O=~{sample_basename}.sorted.bam REFERENCE_SEQUENCE=~{reference_fa} SEQUENCE_DICTIONARY=~{reference_dict}

    java -jar ~{picard_path} MarkDuplicates  I=~{sample_basename}.sorted.bam O=~{sample_basename}.marked.bam M=~{sample_basename}.metrics.txt CREATE_INDEX=true
  }

  output {
    File output_bam = "~{sample_basename}.marked.bam"
    File output_bam_index = "~{sample_basename}.marked.bai"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
  }
}

# Generate sets of intervals for scatter-gathering over chromosomes
task CreateSequenceGroupingTSV {
  input {
    File ref_dict
    Int preemptible_tries
  }
  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  command <<<
    python <<CODE
    with open("~{ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
    CODE
  >>>
  runtime {
    preemptible: preemptible_tries
    docker: "us.gcr.io/broad-gotc-prod/python:2.7"
    memory: "2 GiB"
  }
  output {
    Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
    Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")
  }
}


task BaseRecalibrator {
  input {
    File input_bam
    File input_bam_index
    String recalibration_report_filename
    Array[String] sequence_group_interval
    File dbsnp_vcf
    File dbsnp_vcf_index
    Array[File] known_indels_sites_vcfs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int bqsr_scatter
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.10.1"
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Float dbsnp_size = size(dbsnp_vcf, "GiB")
  Int disk_size = ceil((size(input_bam, "GiB") / bqsr_scatter) + ref_size + dbsnp_size) + 20

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command {
    gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
      -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
      -Xloggc:gc_log.log -Xms5g" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{recalibration_report_filename} \
      --known-sites ~{dbsnp_vcf} \
      --known-sites ~{sep=" -known-sites " known_indels_sites_vcfs} \
      -L ~{sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "6 GiB"
    disks: "local-disk " + disk_size + " HDD"
    maxRetries: 3
  }
  output {
    File recalibration_report = "~{recalibration_report_filename}"
  }
}

# Apply Base Quality Score Recalibration (BQSR) model
task ApplyBQSR {
  input {
    File input_bam
    File input_bam_index
    String output_bam_basename
    File recalibration_report
    Array[String] sequence_group_interval
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    Int compression_level
    Int bqsr_scatter
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.10.1"
    Int memory_multiplier = 1
    Int additional_disk = 20
    Boolean bin_base_qualities = true
    Boolean somatic = false
  }

  Float ref_size = size(ref_fasta, "GiB") + size(ref_fasta_index, "GiB") + size(ref_dict, "GiB")
  Int disk_size = ceil((size(input_bam, "GiB") * 3 / bqsr_scatter) + ref_size) + additional_disk

  Int memory_size = ceil(3500 * memory_multiplier)

  Boolean bin_somatic_base_qualities = bin_base_qualities && somatic

  parameter_meta {
    input_bam: {
      localization_optional: true
    }
  }

  command {
    gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
      -XX:+PrintGCDetails -Xloggc:gc_log.log \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=~{compression_level} -Xms3000m" \
      ApplyBQSR \
      --create-output-bam-md5 \
      --add-output-sam-program-record \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{output_bam_basename}.bam \
      -bqsr ~{recalibration_report} \
      ~{true='--static-quantized-quals 10' false='' bin_base_qualities} \
      ~{true='--static-quantized-quals 20' false='' bin_base_qualities} \
      ~{true='--static-quantized-quals 30' false='' bin_base_qualities} \
      ~{true='--static-quantized-quals 40' false='' bin_somatic_base_qualities} \
      ~{true='--static-quantized-quals 50' false='' bin_somatic_base_qualities} \
      -L ~{sep=" -L " sequence_group_interval}
  }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    maxRetries: 3
    memory: "~{memory_size} MiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File recalibrated_bam = "~{output_bam_basename}.bam"
    File recalibrated_bam_checksum = "~{output_bam_basename}.bam.md5"
  }
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
task GatherBqsrReports {
  input {
    Array[File] input_bqsr_reports
    String output_report_filename
    Int preemptible_tries
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.0.10.1"
  }

  command {
    gatk --java-options "-Xms3000m" \
      GatherBQSRReports \
      -I ~{sep=' -I ' input_bqsr_reports} \
      -O ~{output_report_filename}
    }
  runtime {
    docker: gatk_docker
    preemptible: preemptible_tries
    memory: "3500 MiB"
    disks: "local-disk 20 HDD"
  }
  output {
    File output_bqsr_report = "~{output_report_filename}"
  }
}

# Combine multiple *sorted* BAM files
task GatherSortedBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename
    Float total_input_size
    Int compression_level
    Int preemptible_tries
  }

  # Multiply the input bam size by two to account for the input and output
  Int disk_size = ceil(2 * total_input_size) + 20

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms2000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true
    }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    preemptible: preemptible_tries
    maxRetries: 3
    memory: "3 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}


task RecalibrateBAM {
  input {
  File input_bam
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  File dbsnp_vcf
  File dbsnp_vcf_index
  Array[File] known_indels_sites_vcfs
  Array[File] known_indels_sites_indices

  # Runtime params
  String gatk_path
  Int threads
  String docker
  }
  
  command {
  set -e
    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      BaseRecalibrator \
      -R ~{reference_fa} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{sample_basename}.sample.recal_data.grp \
      --known-sites ~{dbsnp_vcf} \
      --known-sites ~{sep=" -known-sites " known_indels_sites_vcfs}

      ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      ApplyBQSR \
      --create-output-bam-md5 \
      --add-output-sam-program-record \
      -R ~{reference_fa} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{sample_basename}.recalibrated.bam \
      -bqsr ~{sample_basename}.sample.recal_data.grp 
  }

  output {
    File output_bam = "~{sample_basename}.recalibrated.bam"
    File output_bai = "~{sample_basename}.recalibrated.bai"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
  }
}

task HaplotypeCaller {
  input {
  File input_bam
  File input_bai
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  String chromosome

  # Runtime params
  String gatk_path
  Int threads
  String docker
  }

  String chromosome_intervals = sub(chromosome, ",", " -L ")
  
  command {
  set -e
    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}" \
      HaplotypeCaller \
      -R ~{reference_fa} \
      -stand-call-conf 50.0 \
      -L ~{chromosome_intervals} \
      -I ~{input_bam} \
      -O ~{sample_basename}.raw.vcf.gz
  }

  output {
    File output_vcf = "~{sample_basename}.raw.vcf.gz"
    File output_vcf_index = "~{sample_basename}.raw.vcf.gz.tbi"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
  }
}

# Merge GVCFs generated per-interval for the same sample
task MergeVCFs {
  input {
    # Command parameters
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String sample_basename

    # Runtime parameters
    String gatk_path
    String docker
  }
  
  command {
  set -e

    ~{gatk_path} --java-options -Xmx4G  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{sample_basename}.raw.vcf.gz
  }
  runtime {
    docker: docker
    maxRetries: 3
  }
  output {
    File output_vcf = "~{sample_basename}.raw.vcf.gz"
    File output_vcf_index = "~{sample_basename}.raw.vcf.gz.tbi"
  }
}

task SplitSNPindel {
  input {
  File input_vcf
  File input_vcf_index
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  # Runtime params
  String gatk_path
  Int threads
  String docker
  }
  
  command {
  set -e
    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      SelectVariants \
      -R ~{reference_fa} \
      -V ~{input_vcf} \
      -O ~{sample_basename}.raw.snp.vcf \
      -select-type SNP

      ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      SelectVariants \
      -R ~{reference_fa} \
      -V ~{input_vcf} \
      -O ~{sample_basename}.raw.indel.vcf \
      -select-type INDEL

  }

  output {
    File output_vcf_SNP = "~{sample_basename}.raw.snp.vcf"
    File output_vcf_INDEL = "~{sample_basename}.raw.indel.vcf"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
  }
}

task VariantFiltrationSNP {
  input {
  File input_vcf
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  # Runtime params
  String gatk_path
  Int threads
  String docker
  }
  
  command {
  set -e
    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      VariantFiltration \
      -R ~{reference_fa} \
      -V ~{input_vcf} \
      -O ~{sample_basename}.filtered.snp.vcf \
      --filter-expression "DP < 5" --filter-name "LowCoverage" \
      --filter-expression "QUAL < 30.0" --filter-name "VeryLowQual" \
      --filter-expression "FS > 60.0" --filter-name "StrandBias" \
      --filter-expression "SOR>5.0" --filter-name "HighStrandsOR"
  }

  output {
    File output_vcf = "~{sample_basename}.filtered.snp.vcf"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
  }
}

task VariantFiltrationINDEL {
  input {
  File input_vcf
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  # Runtime params
  String gatk_path
  Int threads
  String docker
  }
  
  command {
  set -e
    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}" \
      VariantFiltration \
      -R ~{reference_fa} \
      -V ~{input_vcf} \
      -O ~{sample_basename}.filtered.indel.vcf \
      --filter-expression "DP < 5" --filter-name "LowCoverage" \
      --filter-expression  "QUAL < 30.0" --filter-name "VeryLowQual" \
      --filter-expression  "FS > 200.0" --filter-name "StrandBias" \
      --filter-expression  "SOR>10.0" --filter-name "HighStrandsOR"
  }

  output {
    File output_vcf = "~{sample_basename}.filtered.indel.vcf"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
  }
}

task CombineSNPindel {
  input {
  File input_snp_vcf
  File input_indel_vcf
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  # Runtime params
  String picard_path
  Int threads
  String docker
  }
  
  command {
  set -e
     java -jar ~{picard_path} MergeVcfs I=~{input_snp_vcf} I=~{input_indel_vcf} O=~{sample_basename}.both.vcf
  }

  output {
    File output_vcf = "~{sample_basename}.both.vcf"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
  }
}

task SelectFinalVariants {
  input {
  File input_vcf
  String sample_basename

  File reference_fa
  File reference_fai
  File reference_dict

  # Runtime params
  String gatk_path
  Int threads
  String docker
  }
  
  command {
  set -e
     ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}" \
      SelectVariants \
      -R ~{reference_fa} \
      -V ~{input_vcf} \
      -O ~{sample_basename}.vcf \
      --exclude-filtered
  }

  output {
    File output_vcf = "~{sample_basename}.vcf"
    File output_vcf_index = "~{sample_basename}.vcf.idx"
  }

  runtime {
    docker: "~{docker}"
    maxRetries: 3
  }
}

# Convert a pair of FASTQs to uBAM
task PairedFastQsToUnmappedBAM {
  input {
    # Command parameters
    String sample_name
    File fastq_1
    File fastq_2
    String readgroup_name
    String library_name
    String platform_unit
    String run_date
    String platform_name
    String sequencing_center
    String gatk_path

    # Runtime parameters
    Int addtional_disk_space_gb = 10
    Int machine_mem_gb = 7
    Int preemptible_attempts = 3
    String docker
  }
    Int command_mem_gb = machine_mem_gb - 1
    Int disk_space_gb = ceil((size(fastq_1, "GB") + size(fastq_2, "GB")) * 2 ) + addtional_disk_space_gb
  command {
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}g" \
    FastqToSam \
    --FASTQ ~{fastq_1} \
    --FASTQ2 ~{fastq_2} \
    --OUTPUT ~{readgroup_name}.unmapped.bam \
    --READ_GROUP_NAME ~{readgroup_name} \
    --SAMPLE_NAME ~{sample_name} \
    --LIBRARY_NAME ~{library_name} \
    --PLATFORM_UNIT ~{platform_unit} \
    --RUN_DATE ~{run_date} \
    --PLATFORM ~{platform_name} \
    --SEQUENCING_CENTER ~{sequencing_center} 
  }
  runtime {
    docker: docker
    memory: machine_mem_gb + " GB"
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
    maxRetries: 3
  }
  output {
    File output_unmapped_bam = "~{readgroup_name}.unmapped.bam"
  }
}

task SamSplitter {
  input {
    File input_bam
    Int n_reads
    Int preemptible_tries
    Int compression_level
  }

  Float unmapped_bam_size = size(input_bam, "GiB")
  # Since the output bams are less compressed than the input bam we need a disk multiplier that's larger than 2.
  Float disk_multiplier = 2.5
  Int disk_size = ceil(disk_multiplier * unmapped_bam_size + 20)

  command {
    set -e
    mkdir output_dir

    total_reads=$(samtools view -c ~{input_bam})

    java -Dsamjdk.compression_level=~{compression_level} -Xms3000m -jar /usr/gitc/picard.jar SplitSamByNumberOfReads \
      INPUT=~{input_bam} \
      OUTPUT=output_dir \
      SPLIT_TO_N_READS=~{n_reads} \
      TOTAL_READS_IN_INPUT=$total_reads
  }
  output {
    Array[File] split_bams = glob("output_dir/*.bam")
  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    preemptible: preemptible_tries
    maxRetries: 3
    memory: "3.75 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }
}

# Combine multiple *unsorted* BAM files
# Note that if/when WDL supports optional outputs, we should merge this task with the sorted version
task GatherUnsortedBamFiles {
  input {
    Array[File] input_bams
    String output_bam_basename
    Int compression_level
    Int preemptible_tries
  }

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms2000m -jar /usr/gitc/picard.jar \
      GatherBamFiles \
      INPUT=~{sep=' INPUT=' input_bams} \
      OUTPUT=~{output_bam_basename}.bam \
      CREATE_INDEX=false \
      CREATE_MD5_FILE=false
    }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    preemptible: preemptible_tries
    maxRetries: 3
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
  }
}

# Sort BAM file by coordinate order
task SortSam {
  input {
    File input_bam
    String output_bam_basename
    Int preemptible_tries
    Int compression_level
  }
  # SortSam spills to disk a lot more because we are only store 300000 records in RAM now because its faster for our data so it needs
  # more disk space.  Also it spills to disk in an uncompressed format so we need to account for that with a larger multiplier
  Float sort_sam_disk_multiplier = 6.25
  Int disk_size = ceil(sort_sam_disk_multiplier * size(input_bam, "GiB")) + 20

  command {
    java -Dsamjdk.compression_level=~{compression_level} -Xms4000m -jar /usr/gitc/picard.jar \
      SortSam \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_basename}.bam \
      SORT_ORDER="coordinate" \
      CREATE_INDEX=true \
      CREATE_MD5_FILE=true \
      MAX_RECORDS_IN_RAM=300000

  }
  runtime {
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.3-1564508330"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    memory: "5000 MiB"
    preemptible: preemptible_tries
    maxRetries: 3
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
  }
}




