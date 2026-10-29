version 1.0

workflow SVTopo {
  input {
    String sample_id
    File bam
    File? bam_index
    File? vcf
    File? supporting_reads_json
    File? reference_fasta
    File? exclude_regions
    String assembly = "hg19"
    Int threads = 4
    Boolean allow_unphased = false
  }

  call run_svtopo {
    input:
      sample_id = sample_id,
      bam = bam,
      vcf = vcf,
      supporting_reads_json = supporting_reads_json,
      reference_fasta = reference_fasta,
      exclude_regions = exclude_regions,
      assembly = assembly,
      threads = threads,
      allow_unphased = allow_unphased
  }

  output {
    File svtopo_archive = run_svtopo.svtopo_archive
  }
}

task run_svtopo {
  input {
    String sample_id
    File bam
    File? vcf
    File? supporting_reads_json
    File? reference_fasta
    File? exclude_regions
    String assembly = "hg19"
    Int threads = 4
    Boolean allow_unphased = false
  }

  String results_dir = sample_id + ".svtopo"
  String tarball_name = sample_id + ".svtopo.tar.gz"

  command <<<
    set -euo pipefail

    mkdir -p "~{results_dir}"
    workdir="~{results_dir}/work"
    mkdir -p "$workdir"

    bam="~{bam}"
    reference="~{reference_fasta}"

    if [[ "$bam" == *.cram ]]; then
      if [[ -z "$reference" ]]; then
        echo "[run_svtopo] CRAM input requires a reference fasta" >&2
        exit 1
      fi

      if [[ ! -f "${reference}.fai" ]]; then
        echo "[run_svtopo] Reference index not found, creating it in workdir" >&2
        ref_local="$workdir/$(basename "$reference")"
        if [[ ! -f "$ref_local" ]]; then
          ln -s "$reference" "$ref_local"
        fi
        samtools faidx "$ref_local"
        reference="$ref_local"
      fi

      bam="$workdir/~{sample_id}.~{assembly}.bam"
      if [[ ! -f "$bam" ]]; then
        echo "[run_svtopo] Converting CRAM to BAM: ~{bam} -> $bam" >&2
        samtools view -@ "~{threads}" -b -T "$reference" -o "$bam" "~{bam}"
        samtools index -@ "~{threads}" "$bam"
      fi
    elif [[ ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
      echo "[run_svtopo] Indexing BAM: $bam" >&2
      samtools index -@ "~{threads}" "$bam"
    fi

    exclude_regions="~{exclude_regions}"
    if [[ -z "$exclude_regions" ]]; then
      exclude_regions="$workdir/cnv.excluded_regions.~{assembly}.bed.gz"
      if [[ ! -f "$exclude_regions" ]]; then
        echo "[run_svtopo] Downloading default exclude-regions BED for ~{assembly}" >&2
        wget -q -O "$exclude_regions" \
          "https://github.com/PacificBiosciences/HiFiCNV/raw/refs/heads/main/data/excluded_regions/cnv.excluded_regions.~{assembly}.bed.gz"
      fi
    fi

    svtopo_args=(
      --bam "$bam"
      --svtopo-dir "~{results_dir}"
      --prefix "~{sample_id}"
      --exclude-regions "$exclude_regions"
    )

    if [[ -n "~{vcf}" ]]; then
      svtopo_args+=(--vcf "~{vcf}")
    fi
    if [[ -n "~{supporting_reads_json}" ]]; then
      svtopo_args+=(--variant-readnames "~{supporting_reads_json}")
    fi
    if [[ "~{allow_unphased}" == "true" ]]; then
      svtopo_args+=(--allow-unphased)
    fi

    echo "[run_svtopo] Running: svtopo ${svtopo_args[*]}" >&2
    svtopo "${svtopo_args[@]}"

    echo "[run_svtopo] Running: svtopovz --svtopo-dir ~{results_dir}" >&2
    svtopovz --svtopo-dir "~{results_dir}"

    tar -czf "~{tarball_name}" "~{results_dir}"
  >>>

  runtime {
    docker: "alesmaver/svtopo"
    requested_memory_mb_per_core: 2000
    cpu: 4
    runtime_minutes: 720
  }

  output {
    File svtopo_archive = "~{tarball_name}"
  }
}
