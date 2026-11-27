version 1.0

# Comment by TM (2025-04-18)
# VEP annotation program has a bug: certain vcf files produced by Deep Variant program cannot be annotated. The VEP program 
# prints out the error message, but the program "runs" further. The whole container is then stopped by the system schedular
# when the time-limit is reached. 
# If this happens, the VEP annotated vcf index file is missing and the corresponding vcf.gz file is only partially built.
# Previously a Cleanup task then replaced the partial file with a header-only VCF, so a timeout looked like a successful
# run with 0 variants (PX15640, 2026-11). Removed: RunVEP now fails on timeout or when output record count != input,
# is retried once, and otherwise fails the workflow.


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
            --fork 48 --offline --format vcf --vcf --force_overwrite --compress_output bgzip -v \
            --merged \
            --cache --dir_cache /opt/vep/.vep \
            --plugin AlphaMissense,file=/opt/vep/.vep/Plugins/AlphaMissense/AlphaMissense_hg19.tsv.gz \
            --nearest symbol \
            --shift_hgvs 1 \
            --hgvs \
            --allele_number \
            --assembly GRCh37 \
            --everything \
            --hgvsg \
            --no_stats 
        tabix -p vcf ~{annotated_vcf}

        # VEP --vcf writes every input record, so the counts must match. A mismatch means
        # VEP died or hung part-way (truncated output); fail so the run cannot pass silently.
        # An input with 0 records legitimately gives 0 output records.
        n_in=$(zcat -f ~{input_vcf} | grep -vc '^#' || true)
        n_out=$(zcat -f ~{annotated_vcf} | grep -vc '^#' || true)
        echo "Input records: $n_in, annotated records: $n_out"
        if [ "$n_in" -ne "$n_out" ]; then
            echo "ERROR: VEP output has $n_out records, input has $n_in" >&2
            exit 1
        fi
        echo "Finishing VEP analysis."
        ls -ls ~{annotated_vcf}*
    }

    runtime {
        docker: "alesmaver/vep_grch37"
        requested_memory_mb_per_core: 1000
        cpu: 64
        # Most runs take 0:37-1:42; PX15640 (2026-11) exceeded 2 h and silently fell through to Cleanup (header-only VCF).
        runtime_minutes: 240
        maxRetries: 1
    }

    output {
        File  output_vcf = annotated_vcf
        File  output_vcf_index = annotated_vcf + ".tbi"
    }
}
