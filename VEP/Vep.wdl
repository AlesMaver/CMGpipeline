version 1.0

# Comment by TM (2025-04-18)
# VEP annotation program has a bug: certain vcf files produced by Deep Variant program cannot be annotated. The VEP program 
# prints out the error message, but the program "runs" further. The whole container is then stopped by the system schedular
# when the time-limit is reached. 
# If this happens, the VEP annotated vcf index file is missing and the corresponding vcf.gz file is only partially built.
# In this case the Cleanup task is called to remove the partial file and to create two empty files (just for the record).


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

  if ( !defined(RunVEP.output_vcf_index)) {
        call Cleanup as Cleanup {
            input: 
              input_vcf = RunVEP.output_vcf,
              output_filename = select_first([basename(RunVEP.output_vcf), ""])
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
        File  output_vcf = annotated_vcf
        File? output_vcf_index = annotated_vcf + ".tbi"
    }
}

task Cleanup {
    input {
      File? input_vcf
      String? output_filename
    }

    command {
        ## don't use. zgrep on a corrupted file returns code 2 .... set -e
        echo "Running cleanup..."

        ## remove the file and create a new completely empty file and its index
        # rm ~{input_vcf}
        # touch ~{output_filename}
        # touch ~{output_filename}.tbi

        ## grep out all commented lines and put them into a temp file, gzip it and index it.
        zgrep "#" ~{input_vcf} > temp.vcf 
        bgzip -c temp.vcf > ~{output_filename}
        tabix -p vcf ~{output_filename}

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
