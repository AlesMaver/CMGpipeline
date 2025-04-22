version 1.0

# Copyright (c) 2018 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import "./CRAM_conversions.wdl" as CramConversions
import "./VEP/Vep2.wdl" as VEP

workflow DeepVariant {
  input {
    String sample_basename
    File? inputBam
    File? inputBamIndex
    File? inputCram
    File? inputCramIndex
    File  referenceFasta
    File  referenceFastaIndex
    File? referenceDictionary
    String modelType
    # String outputVcf
    Int? numShards
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

  if (defined(inputCram)) {
    call CramConversions.CramToBam as CramToBam {
        input:
          sample_name = sample_basename,
          input_cram = inputCram,
          ref_fasta = referenceFasta,
          ref_fasta_index = referenceFastaIndex,
          ref_dict = referenceDictionary,
          docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817",
          samtools_path = "samtools"
    }
  }

  call RunDeepVariant {
      input:
        referenceFasta = referenceFasta,
        referenceFastaIndex = referenceFastaIndex,
        inputBam=select_first([inputBam, CramToBam.output_bam]),
        inputBamIndex=select_first([inputBamIndex, CramToBam.output_bai]),
        modelType = modelType,
        # outputVcf = outputVcf,
        outputVcf = sample_basename + ".DeepVariant.vcf.gz",
        numShards = numShards,
        sampleName = sample_basename
    }

  call VEP.VEP as VEPDeepVariant {
      input:
        sample_basename = sample_basename,
        input_vcf = RunDeepVariant.outputVCF,
        filename_infix = ".DeepVariant"
  }

  output {
      File outputVCF = RunDeepVariant.outputVCF
      File outputVCFIndex = RunDeepVariant.outputVCFIndex
      File? outputVCFStatsReport = RunDeepVariant.outputVCFStatsReport
      File? outputGVCF = RunDeepVariant.outputGVCF
      File? outputGVCFIndex = RunDeepVariant.outputGVCFIndex
      File? VEPdeepvariantannotatedVCF = VEPDeepVariant.output_vcf
      File? VEPdeepvariantannotatedVCFIndex = VEPDeepVariant.output_vcf_index      
  }
}

task RunDeepVariant {
    input {
        File referenceFasta
        File referenceFastaIndex
        File inputBam
        File inputBamIndex
        String modelType
        String outputVcf
        String? postprocessVariantsExtraArgs
        File? customizedModel
        Int? numShards = 30
        String? outputGVcf
        String? outputGVcfIndex
        File? regions
        String? sampleName
        Boolean? VCFStatsReport = true

        String memory = "3GiB"
        Int timeMinutes = 2800
        String dockerImage = "google/deepvariant:1.0.0"
    }

    command {
        set -e
        /opt/deepvariant/bin/run_deepvariant \
        --ref ~{referenceFasta} \
        --reads ~{inputBam} \
        --model_type ~{modelType} \
        --output_vcf ~{outputVcf} \
        ~{"--output_gvcf " + outputGVcf} \
        ~{"--customized_model " + customizedModel} \
        ~{"--num_shards " + numShards} \
        ~{"--regions "  + regions} \
        ~{"--sample_name " + sampleName} \
        ~{"--postprocess_variants_extra_args " + postprocessVariantsExtraArgs} \
        ~{true="--vcf_stats_report" false="--novcf_stats_report" VCFStatsReport}
    }

    runtime {
        # orig: memory: memory
        # orig: time_minutes: timeMinutes
        # orig: docker: dockerImage

        docker: "google/deepvariant:1.5.0"
        requested_memory_mb_per_core: 1000
        cpu: 30
        runtime_minutes: 720
    }

    output {
        File outputVCF = outputVcf
        File outputVCFIndex = outputVcf + ".tbi"
        ## Array[File] outputVCFStatsReport = glob("*.visual_report.html")
        File? outputVCFStatsReport = sampleName + ".DeepVariant.visual_report.html"
        File? outputGVCF = outputGVcf
        File? outputGVCFIndex = outputGVcfIndex
    }

    parameter_meta {
        # inputs
        referenceFasta: {description: "Genome reference to use.", category: "required"}
        referenceFastaIndex: {description: "Index for the genome reference file.", category: "required"}
        inputBam: {description: "Aligned, sorted, indexed BAM file containing the reads we want to call.", category: "required"}
        inputBamIndex: {description: "Index for the input bam file.", category: "required"}
        modelType: {description: "<WGS|WES|PACBIO>. Type of model to use for variant calling. Each model_type has an associated default model, which can be overridden by the --customized_model flag.", category: "required"}
        outputVcf: {description: "Path where we should write VCF file.", category: "required"}
        postprocessVariantsExtraArgs: {description: "A comma-separated list of flag_name=flag_value. 'flag_name' has to be valid flags for calpostprocess_variants.py.", category: "advanced"}
        customizedModel: {description: "A path to a model checkpoint to load for the `call_variants` step. If not set, the default for each --model_type will be used.", category: "advanced"}
        numShards: {description: "Number of shards for make_examples step.", category: "common"}
        outputGVcf: {description: "Path where we should write gVCF file.", category: "common"}
        outputGVcfIndex: {description: "Path to where the gVCF index file will be written. This is needed as a workaround, set it to `outputGVcf + '.tbi.'`", category: "common"}
        regions: {description: "List of regions we want to process, in BED/BEDPE format.", category: "advanced"}
        sampleName: {description: "Sample name to use instead of the sample name from the input reads BAM (SM tag in the header).", category: "common"}
        VCFStatsReport: {description: "Output a visual report (HTML) of statistics about the output VCF.", category: "common"}
        memory: {description: "The amount of memory this job will use.", category: "advanced"}
        timeMinutes: {description: "The maximum amount of time the job will run in minutes.", category: "advanced"}
        dockerImage: {description: "The docker image used for this task. Changing this may result in errors which the developers may choose not to address.", category: "advanced"}

        # outputs
        outputVCF: {description: "Output VCF file."}
        outputVCFIndex: {description: "Index of output VCF file."}
        outputVCFStatsReport: {description: "Statistics file."}
        outputGVCF: {description: "GVCF version of VCF file(s)."}
        outputGVCFIndex: {description: "Index of GVCF file(s)."}
    }
}
