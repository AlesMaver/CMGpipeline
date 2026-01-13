version 1.0

workflow VEP {
  input {
    String sample_basename
    File input_vcf
    # not needed File input_vcf_index
    String filename_infix = ""
    ### tole uporabi pri klicu call VEP v anotacijah: output_vcf = sample_basename + ".VEP.hg19.annotated.vcf.gz" | sample_basename + ".DeepVariant.VEP.hg19.annotated.vcf.gz" 
  }

String filename_suffix = ".VEP.hg19.annotated.vcf.gz"
call RunVEP {
      input:
        sample_basename = sample_basename,
        input_vcf = input_vcf,
        annotated_vcf = sample_basename + filename_infix + filename_suffix
    }

  output {
      File output_vcf = RunVEP.output_vcf
      File output_vcf_index = RunVEP.output_vcf_index
  }

}

task RunVEP {
    input {
      String sample_basename
      File input_vcf
      String annotated_vcf
    }

    command {
        set -e
        echo ~{annotated_vcf}
        echo "Starting VEP analysis..."
        vep -i ~{input_vcf} \
            -o ~{annotated_vcf} \
            --offline --format vcf --vcf --force_overwrite --compress_output bgzip -v \
            --merged \
            --cache --dir_cache /opt/vep/.vep \
            --plugin AlphaMissense,file=/opt/vep/.vep/Plugins/AlphaMissense/AlphaMissense_hg19.tsv.gz \
            --nearest symbol \
            --shift_hgvs 0 \
            --allele_number \
            --assembly GRCh37 \
            --no_stats \
            --fork 48
        tabix -p vcf ~{annotated_vcf}
        echo "Finishing VEP analysis."
        ls -ls ~{annotated_vcf}*
    } 

    runtime {
        docker: "alesmaver/vep_grch37"
        requested_memory_mb_per_core: 1000
        cpu: 32
        # 60 minutes is enough for an exome. what is the min time for genome?
        runtime_minutes: 360
    }

    output {
        File output_vcf = annotated_vcf
        File output_vcf_index = annotated_vcf + ".tbi"
    }

}
