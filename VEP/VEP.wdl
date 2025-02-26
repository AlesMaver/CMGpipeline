version 1.0

workflow VEP {
  input {
    String sample_basename
    File input_vcf
    File input_vcf_index
    String output_vcf
    ### tole uporabi pri klicu call VEP v anotacijah: output_vcf = sample_basename + ".VEP.annotated.vcf.gz" | sample_basename + ".DeepVariant.VEP.annotated.vcf.gz" 

  }

call RunVEP {
      input:
        sampleName = sample_basename,
        input_vcf = input_vcf,
        input_vcf_index = input_vcf_index,
        output_vcf = output_vcf,



    }

  output {
      File outputVCF = RunDeepVariant.outputVCF
      File outputVCFIndex = RunDeepVariant.outputVCFIndex
      File? outputVCFStatsReport = RunDeepVariant.outputVCFStatsReport
      File? outputGVCF = RunDeepVariant.outputGVCF
      File? outputGVCFIndex = RunDeepVariant.outputGVCFIndex
  }
}

task RunVEP {
    input {
      String sample_basename
      File input_vcf
      File input_vcf_index
      String output_vcf


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
        docker: "alesmaver/vep"
        requested_memory_mb_per_core: 1000
        cpu: 8
        runtime_minutes: 60
    }

    output {
        File outputVCF = outputVcf
        File outputVCFIndex = outputVcf + ".tbi"
    }

}
