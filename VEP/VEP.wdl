version 1.0

workflow VEP {
  input {
    String sample_basename
    File input_vcf
    #File input_vcf_index
    String output_vcf
    ### tole uporabi pri klicu call VEP v anotacijah: output_vcf = sample_basename + ".VEP.hg19.annotated.vcf.gz" | sample_basename + ".DeepVariant.VEP.hg19.annotated.vcf.gz" 

  }

call RunVEP {
      input:
        sampleName = sample_basename,
        input_vcf = input_vcf,
        #input_vcf_index = input_vcf_index,
        output_vcf = annotated_vcf
    }

  output {
      #File outputVCF = RunVEP.outputVCF
      #File outputVCFIndex = RunVEP.outputVCFIndex
  }

}

task RunVEP {
    input {
      String sample_basename
      File input_vcf
      #File input_vcf_index
      String annotated_vcf


    }

    command {
        set -e
        vep -i ~{input_vcf} \
            -o ~{annotated_vcf} \
            --fork 48 --offline --format vcf --vcf --force_overwrite --compress_output bgzip -v \
            --merged \
            --cache --dir_cache /opt/vep/.vep \
            --plugin AlphaMissense,file=/opt/vep/.vep/Plugins/AlphaMissense/AlphaMissense_hg19.tsv.gz \
            --nearest symbol \
            --shift_hgvs 0 \
            --allele_number \
            --assembly GRCh37 \
            --no_stats 
        tabix -p vcf ~{annotated_vcf}
        ls -ls ~{annotated_vcf}*
    }

    runtime {
        docker: "alesmaver/vep"
        requested_memory_mb_per_core: 1000
        cpu: 8
        runtime_minutes: 60
    }

    output {
        #File outputVCF = outputVcf
        #File outputVCFIndex = outputVcf + ".tbi"
    }

}
