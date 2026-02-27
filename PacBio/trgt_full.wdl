version 1.0

workflow trgt {
  meta {
     title: "TRGT workflow"
     summary: "Tandem Repeat Genotyping Tool (TRGT) workflow"
     description: "Workflow to genotype tandem repeats using TRGT. The workflow extracts reads overlapping repeat regions, performs genotyping, and generates both allele and waterfall plots including methylation visualization."
  }

  parameter_meta {
    sample_name: {
      help: "Name of the sample being processed",
      label: "Sample name"
    }
    sex: {
      help: "Biological sex of the sample",
      label: "Sample sex",
      choices: ["M", "F"]
    }
    mapped_bam: {
      help: "Input BAM file containing aligned reads",
      label: "Mapped BAM"
    }
    mapped_bam_bai: {
      help: "Index file for the input BAM",
      label: "Mapped BAM index"
    }
    ref_fasta: {
      help: "Reference genome in FASTA format",
      label: "Reference FASTA"
    }
    ref_index: {
      help: "Index file for the reference FASTA",
      label: "Reference FASTA index"
    }
    trgt_bed: {
      help: "BED file containing tandem repeat regions to analyze",
      label: "TRGT tandem repeat catalog BED"
    }
    reads_overlapping_repeats: {
      help: "BAM file containing reads that overlap with repeat regions",
      label: "Repeat region reads BAM"
    }
    reads_overlapping_repeats_bai: {
      help: "Index for the BAM containing reads overlapping repeat regions",
      label: "Repeat region reads BAM index"
    }
    spanning_bam: {
      help: "Output BAM containing reads spanning repeat regions",
      label: "TRGT spanning reads BAM"
    }
    spanning_bam_bai: {
      help: "Index for the spanning reads BAM",
      label: "TRGT spanning reads BAM index"
    }
    vcf: {
      help: "Output VCF containing repeat genotypes",
      label: "TRGT repeats VCF"
    }
    images_motifs_allele: {
      help: "Archive containing SVG plots showing repeat motifs in allele view",
      label: "Motif allele plots"
    }
    images_meth_allele: {
      help: "Archive containing SVG plots showing methylation in allele view",
      label: "Methylation allele plots"
    }
    images_motifs_waterfall: {
      help: "Archive containing SVG plots showing repeat motifs in waterfall view",
      label: "Motif waterfall plots"
    }
    images_meth_waterfall: {
      help: "Archive containing SVG plots showing methylation in waterfall view",
      label: "Methylation waterfall plots"
    }
  }

  input {
    String sample_name
    String sex
    File mapped_bam
    File mapped_bam_bai

    File ref_fasta
    File ref_index

    File trgt_bed

    String docker_smrttools
  }

  call get_reads_overlapping_repeats {
    input:
      sample_name      = sample_name,
      mapped_bam       = mapped_bam,
      mapped_bam_bai   = mapped_bam_bai,
      trgt_bed         = trgt_bed,
      docker_smrttools = docker_smrttools
  }

  call trgt_genotype {
    input:
      sample_name      = sample_name,
      sex              = sex,
      mapped_bam       = mapped_bam,
      mapped_bam_bai   = mapped_bam_bai,
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      trgt_bed         = trgt_bed,
      docker_smrttools = docker_smrttools
  }

  call trgt_plot as trgt_plot_motifs_allele {
    input:
      sample_name      = sample_name,
      vcf              = trgt_genotype.vcf,
      spanning_bam     = trgt_genotype.spanning_bam,
      spanning_bam_bai = trgt_genotype.spanning_bam_bai,
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      trgt_bed         = trgt_bed,
      methylation      = false,
      waterfall        = false,
      docker_smrttools = docker_smrttools
  }

  call trgt_plot as trgt_plot_meth_allele {
    input:
      sample_name      = sample_name,
      vcf              = trgt_genotype.vcf,
      spanning_bam     = trgt_genotype.spanning_bam,
      spanning_bam_bai = trgt_genotype.spanning_bam_bai,
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      trgt_bed         = trgt_bed,
      methylation      = true,
      waterfall        = false,
      docker_smrttools = docker_smrttools
  }

  call trgt_plot as trgt_plot_motifs_waterfall {
    input:
      sample_name      = sample_name,
      vcf              = trgt_genotype.vcf,
      spanning_bam     = trgt_genotype.spanning_bam,
      spanning_bam_bai = trgt_genotype.spanning_bam_bai,
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      trgt_bed         = trgt_bed,
      methylation      = false,
      waterfall        = true,
      docker_smrttools = docker_smrttools
  }

  call trgt_plot as trgt_plot_meth_waterfall {
    input:
      sample_name      = sample_name,
      vcf              = trgt_genotype.vcf,
      spanning_bam     = trgt_genotype.spanning_bam,
      spanning_bam_bai = trgt_genotype.spanning_bam_bai,
      ref_fasta        = ref_fasta,
      ref_index        = ref_index,
      trgt_bed         = trgt_bed,
      methylation      = true,
      waterfall        = true,
      docker_smrttools = docker_smrttools
  }

  output {
    File vcf                             = trgt_genotype.vcf
    File spanning_bam                    = trgt_genotype.spanning_bam
    File spanning_bam_bai                = trgt_genotype.spanning_bam_bai
    File images_motifs_allele            = trgt_plot_motifs_allele.images
    File images_meth_allele              = trgt_plot_meth_allele.images
    File images_motifs_waterfall         = trgt_plot_motifs_waterfall.images
    File images_meth_waterfall           = trgt_plot_meth_waterfall.images
    File reads_overlapping_repeats       = get_reads_overlapping_repeats.reads_overlapping_repeats
    File reads_overlapping_repeats_bai   = get_reads_overlapping_repeats.reads_overlapping_repeats_bai
  }
}

task get_reads_overlapping_repeats {
  meta {
    title: "Extract reads overlapping repeat regions"
    summary: "Extracts reads from BAM file that overlap with specified repeat regions"
    description: "Uses samtools to extract reads that overlap with expanded repeat regions, creating a subset BAM for downstream analysis"
  }

  parameter_meta {
    sample_name: {
      help: "Name of the sample being processed",
      label: "Sample name"
    }
    mapped_bam: {
      help: "Input BAM file containing aligned reads",
      label: "Mapped BAM"
    }
    mapped_bam_bai: {
      help: "Index file for the input BAM",
      label: "Mapped BAM index"
    }
    trgt_bed: {
      help: "BED file containing tandem repeat regions to analyze",
      label: "TRGT tandem repeat catalog BED"
    }
    expand_size: {
      help: "Padding size in bases to expand repeat regions by (default 3000)",
      label: "Region expansion padding size"
    }
    reads_overlapping_repeats: {
      help: "BAM file containing reads that overlap with repeat regions",
      label: "Repeat region reads BAM"
    }
    reads_overlapping_repeats_bai: {
      help: "Index for the BAM containing reads overlapping repeat regions",
      label: "Repeat region reads BAM index"
    }
  }

  input {
    String sample_name
    File mapped_bam
    File mapped_bam_bai

    File trgt_bed

    Int expand_size = 3000

    Int threads = 1
    Int mem_gb  = 4

    String docker_smrttools
  }

  String output_name = "~{sample_name}.repeats.bam"
  String index_name = "~{output_name}.bai"
  String output_name_with_index = "~{output_name}##idx##~{index_name}"
  String expanded_bed_name = "~{sample_name}.expanded.bed"

  Int disk_size = ceil(size(mapped_bam, 'GB') * 2 + 10)

  command <<<
    set -e

    awk -v sz=~{expand_size} \
      '{if ($2 < sz) { print $1"\t0\t"$3 + sz} else {print $1"\t"$2 - sz"\t"$3 + sz} }' \
      ~{trgt_bed} \
      > ~{expanded_bed_name}

    samtools view \
      --write-index \
      --use-index \
      --target-file ~{expanded_bed_name} \
      --output ~{output_name_with_index} \
      ~{mapped_bam}
  >>>

  output {
    File reads_overlapping_repeats = "~{output_name}"
    File reads_overlapping_repeats_bai = "~{index_name}"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
  }
}

task trgt_genotype {
  meta {
    title: "TRGT genotyper"
    summary: "Genotype tandem repeats using TRGT"
    description: "Performs genotyping of tandem repeat regions using TRGT."
  }

  parameter_meta {
    sample_name: {
      help: "Name of the sample being processed",
      label: "Sample name"
    }
    sex: {
      help: "Biological sex of the sample",
      label: "Sample sex",
      choices: ["M", "F"]
    }
    mapped_bam: {
      help: "Input BAM file containing aligned reads",
      label: "Mapped BAM"
    }
    mapped_bam_bai: {
      help: "Index file for the input BAM",
      label: "Mapped BAM index"
    }
    ref_fasta: {
      help: "Reference genome in FASTA format",
      label: "Reference FASTA"
    }
    ref_index: {
      help: "Index file for the reference FASTA",
      label: "Reference FASTA index"
    }
    trgt_bed: {
      help: "BED file containing tandem repeat regions to analyze",
      label: "TRGT tandem repeat catalog BED"
    }
    threads: {
      help: "Number of CPU threads to use (default: 8)",
      label: "CPU threads"
    }
    mem_gb: {
      help: "Memory allocation in gigabytes (default: 16)",
      label: "Memory (GB)"
    }
    vcf: {
      help: "Output VCF containing repeat genotypes",
      label: "TRGT repeats VCF"
    }
    spanning_bam: {
      help: "Output BAM containing reads spanning repeat regions",
      label: "TRGT spanning reads BAM"
    }
    spanning_bam_bai: {
      help: "Index for the spanning reads BAM",
      label: "TRGT spanning reads BAM index"
    }
  }

  input {
    String sample_name
    String? sex
    File mapped_bam
    File mapped_bam_bai

    File ref_fasta
    File ref_index

    File trgt_bed

    Int threads = 8
    Int mem_gb  = 16

    String docker_smrttools
  }

  String out_prefix = "~{sample_name}.trgt"
  String karyotype = if select_first([sex, "F"]) == "M" then "XY" else "XX"

  Int disk_size = ceil((size(mapped_bam, 'GB') + size(ref_fasta, 'GB')) * 2 + 20)

  command <<<
    set -e

    trgt --verbose genotype \
      --threads ~{threads} \
      --preset targeted \
      --genome ~{ref_fasta} \
      --karyotype ~{karyotype} \
      --reads ~{mapped_bam} \
      --repeats ~{trgt_bed} \
      --output-prefix ~{out_prefix}

    gunzip -c "~{out_prefix}.vcf.gz" > ~{out_prefix}.vcf
    samtools sort -o ~{out_prefix}.sorted.spanning.bam ~{out_prefix}.spanning.bam
    samtools index ~{out_prefix}.sorted.spanning.bam
  >>>

  output {
    File vcf              = "~{out_prefix}.vcf"
    File spanning_bam     = "~{out_prefix}.sorted.spanning.bam"
    File spanning_bam_bai = "~{out_prefix}.sorted.spanning.bam.bai"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
  }
}

task trgt_plot {
  meta {
    title: "TRGT plotter"
    summary: "Generates visualization plots for TRGT results"
    description: "Creates SVG plots for each repeat region showing either motif composition or methylation status in either allele or waterfall view styles. Outputs are bundled into a ZIP archive."
  }

  parameter_meta {
    sample_name: {
      help: "Name of the sample being processed",
      label: "Sample name"
    }
    vcf: {
      help: "Input VCF containing repeat genotypes from TRGT",
      label: "TRGT VCF"
    }
    spanning_bam: {
      help: "BAM file containing reads spanning repeat regions",
      label: "Spanning reads BAM"
    }
    spanning_bam_bai: {
      help: "Index for the spanning reads BAM",
      label: "Spanning reads BAM index"
    }
    ref_fasta: {
      help: "Reference genome in FASTA format",
      label: "Reference FASTA"
    }
    ref_index: {
      help: "Index file for the reference FASTA",
      label: "Reference FASTA index"
    }
    trgt_bed: {
      help: "BED file containing tandem repeat regions to analyze",
      label: "TRGT tandem repeat catalog BED"
    }
    methylation: {
      help: "Whether to show methylation status instead of motifs (default: false)",
      label: "Show methylation"
    }
    waterfall: {
      help: "Whether to use waterfall plot style instead of allele view (default: false)",
      label: "Use waterfall plot"
    }
    threads: {
      help: "Number of CPU threads to use (default: 1)",
      label: "CPU threads"
    }
    mem_gb: {
      help: "Memory allocation in gigabytes (default: 4)",
      label: "Memory (GB)"
    }
    images: {
      help: "Archive containing SVG plot files",
      label: "Plot archive"
    }
  }

  input {
    String sample_name
    File vcf
    File spanning_bam
    File spanning_bam_bai

    File ref_fasta
    File ref_index

    File trgt_bed

    Boolean methylation = false
    Boolean waterfall = false

    Int threads = 1
    Int mem_gb  = 4

    String docker_smrttools
  }

  String show          = if (methylation) then "meth" else "motifs"
  String plot_type     = if (waterfall) then "waterfall" else "allele"
  String output_prefix = "~{sample_name}.~{show}_~{plot_type}"

  Int disk_size = ceil((size(spanning_bam, 'GB') + size(ref_fasta, 'GB')) * 2 + 20)

  # Note: bcftools is provided by the SMRT Link toolchain in the container
  command <<<
    bcftools query -f '%INFO/TRID\n' ~{vcf} | while read -r i; do
      trgt plot \
        --repeat-id "$i" \
        --genome ~{ref_fasta} \
        --repeats ~{trgt_bed} \
        --spanning-reads ~{spanning_bam} \
        --vcf ~{vcf} \
        --show ~{show} \
        --plot-type ~{plot_type} \
        --image "~{output_prefix}.$(printf '%s' "$i").trgt_plot.svg"
    done;
    python3 -m zipfile -c ~{output_prefix}.trgt_plots.zip ./*.trgt_plot.svg
  >>>

  output {
    File images = "~{output_prefix}.trgt_plots.zip"
  }

  runtime {
    docker: docker_smrttools
    cpu: threads
    memory: "~{mem_gb} GB"
    disk: disk_size + " GB"
    disks: "local-disk " + disk_size + " SSD"
  }
}