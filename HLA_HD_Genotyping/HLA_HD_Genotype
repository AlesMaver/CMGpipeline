version 1.0

#Author: Aleksander Turk
#This is a Cromwell WDL to run HLA-HD

workflow HLAHD_genotyping_workflow {
    input {
        File input_fq1
        File input_fq2
        String sample_basename

    }

    call HLAHD_genotyping {
        input:
            input_fq1 = input_fq1,
            input_fq2 = input_fq2,
            sample_basename = sample_basename
    }

    output {
        #File hlahd_result = HLAHD_genotyping.result
    }
}

task HLAHD_genotyping {
    input {
        File input_fq1
        File input_fq2
        String sample_basename
    }

    command {
        bash hlahd.sh  -t 8  ~{input_fq1}  ~{input_fq2}  HLA_gene.split.txt  dictionary  ~{sample_basename} .
        cp -p ./~{sample_basename}/result/~{sample_basename}_final.result.txt ./~{sample_basename}.HLA_HD_genotype.final_result.txt
    }

    output {
        #File result = "./${sample_basename}/result/${sample_basename}_final.result.txt"
    }

    runtime {
        docker: "aleksanderturk/hlahd_image:latest"
        requested_memory_mb_per_core: 1000
        cpu: 8
        runtime_minutes: 240
        # memory: "16 GB"
    }
}
