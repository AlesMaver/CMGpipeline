version 1.0

workflow Stripy {
    input {
        String sample_id
        File? bam_file
        File? bai_file
        File reference_fasta
        #File? reference_fasta_index
        # sex: male / female (case sensitive)
        String sex
        String reference_genome_name = "hg19"
    }


    call run_stripy {
        input:
            sample_id = sample_id,
            reference = reference_fasta,
            output = output_directory,
            genome = reference_genome_name,
            sex = if defined(sex) && (sex == "male" || sex == "female") then sex else "male",
            input_file = bam_file
    }
}


task run_stripy {
    input {
        String sample_id
        File reference
        String output
        #String loci
        String genome
        String sex
        File input_file
    }

    command <<<
        set -e

        echo ${sex}
        echo "[ PREPARATION ] Downloading variant catalog JSON"
        wget "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/ExpansionHunter_configuration/variant_catalog.json"
        unset https_proxy
        wget "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/ExpansionHunter_configuration/variant_catalog.json"

        echo "[ PREPARATION ] Preparing LOCI"
        loci=$(jq -r '[.[] | .LocusId] | join(",")' ./variant_catalog.json)

        echo "[ RUNNING ] Stri.py"
        # Constructing Docker run command (inside Docker already)
        ./batch.sh -o ./ -r ${ref_file} -l ${loci} -g ${genome} -s ${sex} -i ${input_file}

        #mv ./{input_file}.json ./sample_id.Stripy.json
        
    >>>

    runtime {
        docker: "gbergant/stripy_prod:2.5"
        requested_memory_mb_per_core: 1000
        cpu: 4
        runtime_minutes: 30
    }

    #output {
    #    File a = "~{sample_basename}.CONIFER_CALLS.txt"
    #}
}
