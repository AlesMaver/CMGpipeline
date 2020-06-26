version 1.0
## Copyright CMG@KIGM, Ales Maver

# WORKFLOW DEFINITION 
workflow Conifer {
  input {
    File input_bam
    File input_bam_index

    Array[File] input_reference_rpkms 
    Int CONIFER_svd

    String enrichment
    File enrichment_bed
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

    String enrichment
    File enrichment_bed

    # Runtime parameters
    String docker = "ielis/conifer:latest"
  }
  
  command {
  set -e
  python /root/conifer/conifer.py rpkm --probes ~{enrichment_bed} --input ~{input_bam} --output ~{enrichment}_~{sample_basename}.txt 
  }
  runtime {
    docker: docker
  }
  output {
    File output_rpkm = "sample_basename.~{enrichment}.rpkm"
  }
}

task CONIFER_Analyze {
  input {
    # Command parameters
    File input_rpkm
    Array[File] input_reference_rpkms 
    String CONIFER_svd
    String sample_basename

    String enrichment
    File enrichment_bed

    # Runtime parameters
    String docker = "ielis/conifer:latest"
  }
  
  command {
  set -e
  RPKMDIR=dirname $(readlink -f ~{input_reference_rpkms[0]})
  echo $RPKMDIR
  cp ~{input_rpkm} $RPKMDIR

  python /root/conifer/conifer.py analyze --probes ~{enrichment_bed} --rpkm_dir $RPKMDIR --output ~{sample_basename}.analysis.hdf5 --svd ~{CONIFER_svd} --write_svals --svd ~{sample_basename}.singular_values.txt --plot_scree ~{sample_basename}.screeplot.png") --write_sd ~{sample_basename}.sd_values.txt
  }
  runtime {
    docker: docker
  }
  output {
    File output_hdf5 = "~{sample_basename}.analysis.hdf5"
  }
}

