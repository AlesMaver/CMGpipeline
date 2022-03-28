version 1.0
## Copyright CMG@KIGM, Ales Maver

# The manta workflow currently holds the annotSV workflow, therefore importing it from there
import "./manta/manta_workflow.wdl" as manta
# Conifer is also used to analyse WGS data as it proved to be useful in detecting rare CNVs while removing the common CNVs as noise
# We use the BED file with 1Mb windows across the genome, prepared using these steps
# The chrY and chrM chromosomes were too short to make accurate calling using conifer
# awk {'print $1, $2'} OFS='\t' hg19.fa.fai |head -n 26 | grep -v '_' | grep -v 'chrM' | grep -v 'chrY' > hg19.genome
# bedtools makewindows -w 1000 -g hg19.genome  > WGS1Mb.bed

# WORKFLOW DEFINITION 
workflow Conifer {
  input {
    File input_bam
    File input_bam_index

    #Array[File]? input_reference_rpkms 
    Array[File] input_reference_rpkms = select_first([input_reference_rpkms, [""]])
    Int? CONIFER_svd
    Float? CONIFER_threshold

    String? enrichment
    File? enrichment_bed
  }  

  String sample_basename = sub(basename(input_bam), "[\_,\.].*", "" )
  
  call MakeRPKM {
      input:
        input_bam=input_bam,
        input_bam_index=input_bam_index,
        sample_basename=sample_basename,
        enrichment=enrichment,
        enrichment_bed=enrichment_bed
  }

  call CONIFER_Analyze {
      input:
        input_rpkm=MakeRPKM.output_rpkm,
        input_reference_rpkms=input_reference_rpkms,
        CONIFER_svd=CONIFER_svd,
        sample_basename=sample_basename,
        enrichment=enrichment,
        enrichment_bed=enrichment_bed

  }

  call CONIFER_Call {
      input:
        input_hdf5=CONIFER_Analyze.output_hdf5,
        CONIFER_threshold=CONIFER_threshold,
        sample_basename=sample_basename
  }

    call CONIFER_Plotcalls {
      input:
        input_hdf5=CONIFER_Analyze.output_hdf5,
        input_conifer_calls=CONIFER_Call.output_conifer_calls,
        sample_basename=sample_basename
  }

  call CONIFER_Export {
      input:
        input_hdf5=CONIFER_Analyze.output_hdf5,
        sample_basename=sample_basename,
        enrichment=enrichment
  }

  call manta.annotSV as annotSV {
      input:
        genome_build = "GRCh37",
        input_vcf = CONIFER_Call.output_conifer_annotSV_input_bed,
        output_tsv_name = sample_basename + "_CONIFER_annotSV.tsv"
  }

  output {
    File output_conifer_calls = CONIFER_Call.output_conifer_calls
    File output_conifer_calls_wig = CONIFER_Call.output_conifer_calls_wig
    Array[File] output_plotcalls = CONIFER_Plotcalls.output_plotcalls
    File CNV_bed = CONIFER_Export.CNV_bed
    File CNV_wig = CONIFER_Export.CNV_wig
    File output_rpkm = MakeRPKM.output_rpkm
    File? annotSV_tsv = annotSV.sv_variants_tsv
  }
}

##################
# TASK DEFINITIONS
##################

task MakeRPKM {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    String sample_basename

    String? enrichment
    File? enrichment_bed

    # Runtime parameters
    String docker = "molecular/conifer"
  }
  
  command {
  set -e
  python /home/bio/conifer_v0.2.2/conifer.py rpkm --probes ~{enrichment_bed} --input ~{input_bam} --output ~{enrichment}_~{sample_basename}.txt 
  }
  runtime {
    docker: docker
    requested_memory_mb_per_core: 3000
    cpu: 2
    runtime_minutes: 500
  }
  output {
    File output_rpkm = "~{enrichment}_~{sample_basename}.txt"
  }
}

task CONIFER_Analyze {
  input {
    # Command parameters
    File input_rpkm
    Array[File] input_reference_rpkms 
    String? CONIFER_svd
    String sample_basename

    String? enrichment
    File? enrichment_bed

    # Runtime parameters
    String docker = "molecular/conifer"
  }

  command <<<
  set -e
  RPKM_DIR=$(dirname ~{input_reference_rpkms[0]})
  cp ~{input_rpkm} $RPKM_DIR

  if [[ ~{enrichment_bed} == *"WGS"* ]]; then MIN_RPKM=" --min_rpkm=0.01 "; else MIN_RPKM=" "; fi

  # Removed this parameter from the analyze command as it was causing issues with python module loads --plot_scree ~{sample_basename}.screeplot.png
  # Need to prefix the python command with HOME= to make home writable in a rootless container
  HOME=$(dirname ~{input_rpkm}) python /home/bio/conifer_v0.2.2/conifer.py analyze $MIN_RPKM --probes ~{enrichment_bed} --rpkm_dir $RPKM_DIR --output ~{sample_basename}.analysis.hdf5 --svd ~{CONIFER_svd} --write_svals ~{sample_basename}.singular_values.txt --plot_scree ~{sample_basename}.screeplot.png --write_sd ~{sample_basename}.sd_values.txt
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 3000
    cpu: 2
    runtime_minutes: 120
    docker_user: "root"
  }
  output {
    File output_hdf5 = "~{sample_basename}.analysis.hdf5"
  }
}

task CONIFER_Call {
  input {
    # Command parameters
    File input_hdf5
    Float? CONIFER_threshold
    String sample_basename

    # Runtime parameters
    String docker = "molecular/conifer"
  }
  
  command <<<
  set -e
  python /home/bio/conifer_v0.2.2/conifer.py call --threshold ~{CONIFER_threshold} --input ~{input_hdf5} --output ~{sample_basename}.CONIFER_CALLS_POPULATION.txt

  echo "Filtering sample specific CNV calls..."
  head -n 1 ~{sample_basename}.CONIFER_CALLS_POPULATION.txt > ~{sample_basename}.CONIFER_CALLS.txt

  # Check if any calls have been generated - otherwise grep returns an error which interrupts the script and stops the call execution
  if grep -q "~{sample_basename}" ~{sample_basename}.CONIFER_CALLS_POPULATION.txt; then
    cat ~{sample_basename}.CONIFER_CALLS_POPULATION.txt | grep ~{sample_basename} >> ~{sample_basename}.CONIFER_CALLS.txt
    cat ~{sample_basename}.CONIFER_CALLS.txt | grep -v "start" | awk -F'\t' '{ if ($5 == "dup") $5="1"; if ($5 == "del") $5="-1";print $2,$3,$4,$5}' OFS='\t' > ~{sample_basename}.CNV.wig
    cat ~{sample_basename}.CONIFER_CALLS.txt | grep -v "start" | awk -F'\t' '{ if ($5 == "dup") $5="DUP"; if ($5 == "DEL") $5="-1";print $2,$3,$4,$5,"~{sample_basename}"}' OFS='\t' > ~{sample_basename}.CNV.annotSV.input.bed
    echo "Sample specific CNV calls generated."
  else
    echo "No CNV calls found for sample"
    echo "Creating an empty WIG file."
    touch ~{sample_basename}.CNV.wig
  fi
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 3000
    cpu: 2
    runtime_minutes: 120
    docker_user: "root"
  }
  output {
    File output_conifer_calls = "~{sample_basename}.CONIFER_CALLS.txt"
    File output_conifer_calls_wig = "~{sample_basename}.CNV.wig"
    File output_conifer_annotSV_input_bed ="~{sample_basename}.CNV.annotSV.input.bed"
  }
}

task CONIFER_Plotcalls {
  input {
    # Command parameters
    File input_hdf5
    File input_conifer_calls
    String sample_basename

    # Runtime parameters
    String docker = "molecular/conifer"
  }
  
  command {
  set -e
  export LC_ALL="C.UTF-8"
  export LC_CTYPE="C.UTF-8"
  # Need to prefix the python command to make home writable in a rootless container
  HOME=$(dirname ~{input_hdf5}) python /home/bio/conifer_v0.2.2/conifer.py plotcalls --input ~{input_hdf5} --calls ~{input_conifer_calls} --output ./
  }

  runtime {
    docker: docker
    requested_memory_mb_per_core: 9000
    cpu: 2
    runtime_minutes: 600
    docker_user: "root"    
  }
  output {
    Array[File] output_plotcalls = glob("*.png")
  }
}

task CONIFER_Export {
  input {
    # Command parameters
    File input_hdf5
    String sample_basename
    String? enrichment

    # Runtime parameters
    String docker = "molecular/conifer"
  }
  
  command <<<
  set -e
  export LC_ALL="C.UTF-8"
  export LC_CTYPE="C.UTF-8"
  python /home/bio/conifer_v0.2.2/conifer.py export --input ~{input_hdf5} --sample ~{enrichment}_~{sample_basename} --output ./
  cat ~{enrichment}_~{sample_basename}.bed | awk -F'\t' '{print $1,$2,$3,$5}' OFS='\t' > ~{sample_basename}.CNV.genome.wig
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 4000
    cpu: 1
    runtime_minutes: 60
    docker_user: "root"
  }
  output {
    File CNV_bed = "~{enrichment}_~{sample_basename}.bed" 
    File CNV_wig = "~{sample_basename}.CNV.genome.wig"
  }
}

task CONIFER_Annotate {
  input {
    # Command parameters
    File conifer_calls_wig

    String sample_basename

    File HPO
    File HPO_index
    File OMIM
    File OMIM_index
    File gnomadConstraints
    File gnomadConstraints_index
    File CGD
    File CGD_index

    # Runtime parameters
    String docker = "pegi3s/bedtools"
  }
  
  command <<<
  set -e
  bedtools intersect -a ~{conifer_calls_wig} -b ~{HPO} ~{OMIM} ~{gnomadConstraints} ~{CGD} -wa -wb -loj | awk '{print $1,$2,$3,$4,$9,$10}' OFS='\t' > ~{sample_basename}.CONIFER.ANNOTATED.txt
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 4000
    cpu: 1
    runtime_minutes: 60
  }
  output {
    File conifer_calls_annotated = "~{sample_basename}.CONIFER.ANNOTATED.txt"
  }
}