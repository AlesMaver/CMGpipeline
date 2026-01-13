version 1.0

import "./CRAM_conversions.wdl" as CramConversions

workflow ExpansionHunter {
  input {
    String sample_id
    File? bam_file
    File? bai_file
    File? cram_file
    File? crai_file
    File reference_fasta
    File? reference_fasta_index
    File? reference_dict
    String expansion_hunter_docker
    String? patient_sex  # Optional input for patient sex
  }

  parameter_meta {
    sample_id: "sample name"
    bam_file: ".bam file to search for repeat expansions"
    reference_fasta: ".fasta file with reference used to align bam file"
    expansion_hunter_docker: "expansion hunter docker including annotation software"
  }

  meta {
      author: "Gaber Bergant and Aleš Maver"
      email: "cmg.kimg@kclj.si"
  }

  if (defined(cram_file)) {
    call CramConversions.CramToBam as CramToBam {
        input:
          sample_name = sample_id,
          input_cram = cram_file,
          ref_fasta = reference_fasta,
          ref_fasta_index = reference_fasta_index,
          ref_dict = reference_dict,
          docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817",
          samtools_path = "samtools"
    }
  }

  call RunExpansionHunter {
      input:
        sample_id = sample_id,
        bam_file = select_first([bam_file, CramToBam.output_bam]),
        bai_file = select_first([bai_file, CramToBam.output_bai]),
        reference_fasta = reference_fasta,
        expansion_hunter_docker = expansion_hunter_docker,
        patient_sex = patient_sex  # Pass the optional input to the task
    }

  call AnnotateExpansionHunter {
      input:
        sample_id = sample_id,
        expansion_hunter_docker = expansion_hunter_docker,
        expansion_hunter_vcf = RunExpansionHunter.expansion_hunter_vcf
    }

  output {
    File? expansion_hunter_vcf_annotated = AnnotateExpansionHunter.expansion_hunter_vcf_annotated
  }

}

task RunExpansionHunter {
  input {
    String sample_id
    File bam_file
    File bai_file
    File reference_fasta
    String expansion_hunter_docker
    String? patient_sex  # Optional input for patient sex
  }

  output {
    File expansion_hunter_vcf = "~{sample_id}.vcf"
  }

  command <<<
    echo ~{sample_id}
    # echo ~{patient_sex}
    echo ' '
    echo "[ PREPARATION ] Downloading variant catalog JSON"
    wget "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/ExpansionHunter_configuration/variant_catalog.json"
    unset https_proxy
    wget "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/ExpansionHunter_configuration/variant_catalog.json"


    echo "[ PREPARATION ] Fixing BAI file ending to BAM.BAI as required for the ExpansionHunter"
    BAIFILE=~{bai_file}
    BAMBAIFILE=${BAIFILE%.bai}.bam.bai
    cp "$BAIFILE" "$BAMBAIFILE"

    echo "[ RUNNING ] expansion hunter on sample ~{sample_id}"
    ExpansionHunter \
      --reads "~{bam_file}" \
      --reference "~{reference_fasta}" \
      --variant-catalog variant_catalog.json \
      --output-prefix "~{sample_id}" \
      ~{if defined(patient_sex) then ("--sex " + patient_sex) else ""}

  >>>
  
  runtime {
    docker: expansion_hunter_docker
    maxRetries: 1
    requested_memory_mb_per_core: 1000
    cpu: 1
    runtime_minutes: 180
  }

}

task AnnotateExpansionHunter {
  input {
    String sample_id
    String expansion_hunter_docker
    File expansion_hunter_vcf
  }
  
  output {
    File expansion_hunter_vcf_annotated = "~{sample_id}.ExpansionHunter.annotated.vcf"
  }
  
  command <<<

    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    # annotated_vcf = "${sample_id}.annotated.vcf"

    echo "[ PREPARATION ] Downloading repeats file JSON"
    wget "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/ExpansionHunter_configuration/variant_catalog_hg19.json"
    unset https_proxy
    wget "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/ExpansionHunter_configuration/variant_catalog_hg19.json"

    echo "[ RUNNING ] expansion hunter vcf annotation on sample ~{sample_id}"
    stranger \
      --repeats-file variant_catalog_hg19.json \
      ~{expansion_hunter_vcf} > ~{sample_id}.ExpansionHunter.annotated.vcf

  >>>

  runtime {
    docker: expansion_hunter_docker
    maxRetries: 1
    requested_memory_mb_per_core: 1000
    cpu: 1
    runtime_minutes: 30
  }

}

