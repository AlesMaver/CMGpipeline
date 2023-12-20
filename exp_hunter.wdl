version 1.0

workflow ExpansionHunter {
  input {
    String sample_id
    File bam_file
    File bai_file
    File reference_fasta
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
  
  call RunExpansionHunter {
      input:
        sample_id = sample_id,
        bam_file = bam_file,
        bai_file = bai_file,
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
      --reads ~{bam_file} \
      --reference ~{reference_fasta} \
      --variant-catalog variant_catalog.json \
      ~{"--sex " + patient_sex} \  # Conditionally add the -s option
      --output-prefix ~{sample_id}

  >>>
  
  runtime {
    docker: expansion_hunter_docker
    maxRetries: 1
    requested_memory_mb_per_core: 1000
    cpu: 1
    runtime_minutes: 30
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
