version 1.0

workflow VEP {
  input {
    String sample_basename
    File input_vcf
    String filename_infix = ""
  }

  String filename_suffix = ".VEP.hg19.annotated.vcf.gz"
  call RunVEP as RunVEP {
      input:
        sample_basename = sample_basename,
        input_vcf = input_vcf,
        annotated_vcf = sample_basename + filename_infix + filename_suffix
  }

  # Get the return code from RunVEP's metadata
  #Int return_code = select(first: select(RunVEP: workflow.RunVEP.metadata.calls))['returnCode']    
  # Conditional execution based on return code: 0 is success, 79 is OTL, if anything else, the runVEP should crash
  # when 79: the runVEP outputs are incorrect, so delete the content of them
  #if (return_code == 79) {
  #      call Cleanup as Cleanup {
  #          input: 
  #            input_vcf = RunVEP.output_vcf,
  #            output_filename = basename(RunVEP.output_vcf)
  #      }
  #}

  if ( !defined(RunVEP.output_vcf_index)) {
        call Cleanup as Cleanup {
            input: 
              input_vcf = RunVEP.output_vcf,
              output_filename = basename(RunVEP.output_vcf)
        }
  }
  output {
      File? output_vcf = select_first([Cleanup.output_vcf, RunVEP.output_vcf]) 
      File? output_vcf_index = select_first([Cleanup.output_vcf_index, RunVEP.output_vcf_index])
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
        echo "Finishing VEP analysis."
        ls -ls ~{annotated_vcf}*
    }

    runtime {
        docker: "alesmaver/vep_grch37"
        continueOnReturnCode: [0, 79]
        requested_memory_mb_per_core: 1000
        cpu: 32
        # 60 minutes is enough for an exome and genome. plus a little save margine (90).
        runtime_minutes: 40
    }

    output {
        File? output_vcf = annotated_vcf
        File? output_vcf_index = annotated_vcf + ".tbi"
    }
}

task Cleanup {
    input {
      File input_vcf
      String output_filename
    }

    command {
        set -e
        echo "Running cleanup..."
        rm ~{input_vcf}
        touch ~{output_filename}
        touch ~{output_filename}.tbi
        ls -ls ~{output_filename}*
    }

    runtime {
        docker: "alesmaver/vep_grch37"
        requested_memory_mb_per_core: 1000
        cpu: 1
        runtime_minutes: 5
    }

    output {
        File? output_vcf = output_filename
        File? output_vcf_index = output_filename + ".tbi"
    }

}
