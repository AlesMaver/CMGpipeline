version 1.0
#Placeholder to connect variables from different .wdl-s
import "x.wdl" as X

workflow AlphaMissense_anno {
    input {
        #input_dir - wherever the input data is stored on the server during processing, define in JSON file
        File input_vcf
        String input_dir
    }
    call AlphaMissense_docker {
        input:
            input_vcf = input_vcf
            input_dir = input_dir
    }

    output {
        File annotated_vcf = AlphaMissense_docker.annotated_vcf
        File annotated_vcf_index = AlphaMissense_docker.annotated_vcf_index
    }
}
task AlphaMissense_docker {
    input {
        File input_vcf
        String input_dir
    }

    command <<<
    /bin/bash -c "
    vep -i /input_data/~{input_vcf} \
        -o /input_data/~{annotated_vcf} \
        --fork 48 --offline --format vcf --vcf --force_overwrite --compress_output bgzip -v \
        --merged \
        --cache --dir_cache /opt/vep/.vep \
        --plugin AlphaMissense,file=/opt/vep/.vep/Plugins/AlphaMissense/AlphaMissense_hg19.tsv.gz \
        --nearest symbol \
        --shift_hgvs 0 \
        --allele_number \
        --assembly GRCh37 \
        --no_stats && \
    tabix -p vcf /input_data/~{annotated_vcf}"
    >>>
    

    output {
        #These still need to have their names defined better - currently they're placeholders
        File annotated_vcf = "~{input_dir}/~{basename(input_vcf)}_DockerVEP.vcf.gz"
        File annotated_vcf_index =  "~{input_dir}/~{basename(input_vcf)}_DockerVEP.vcf.gz.tbi"
    }

    runtime {
        #Note: never mount to /data, it's a reserved directory by Docker
        docker: "alesmaver/vep_grch37"
        volumes: "${input_dir}:/input_data"
        #Given that the VEP AM plugin is fairly resource-hungry, we could increase this
        memory: "4 GB"
        cpu: 2
    }
}