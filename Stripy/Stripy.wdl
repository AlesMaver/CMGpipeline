version 1.0

workflow stripy_workflow {
    input {
        String sample_id
        File? bam_file
        File? bai_file
        File reference_fasta
        File? reference_fasta_index

        String sex
        # male / female (case sensitive)

        String output_directory
        String reference_genome_name = "hg19"
    }

    call extract_loci {

    }

    call run_stripy {
        input:
            reference_fasta = reference_fasta,
            output = output_directory,
            loci = extract_loci.loci_string,
            genome = reference_fasta_name,
            sex = sex,
            bam_file = bam_file
    }
}

task extract_loci {

    command {
        echo "[ PREPARATION ] Downloading variant catalog JSON"
        wget "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/ExpansionHunter_configuration/variant_catalog.json"
        unset https_proxy
        wget "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/ExpansionHunter_configuration/variant_catalog.json"

        jq -r '[.[] | .LocusId] | join(",")' ./variant_catalog.json
    }

    output {
        String loci_string = read_string(stdout())
    }

    runtime {
        docker: "stedolan/jq"
    }
}

task run_stripy {
    input {
        File reference
        String output
        String loci
        String genome
        String sex
        File bam_file
    }

    command <<<
        set -e

        # Path to required files for docker volumes
        ref_dir=$(dirname "${reference_fasta}")
        input_dir=$(dirname "${input_path}")

        # Filenames
        ref_file=$(basename "${reference_path}")
        input_file=$(basename "${input_path}")

        # Constructing Docker run command (inside Docker already)
        ./batch.sh -o /mnt/results -r /mnt/ref/${ref_file} -l ${loci} -g ${genome} -s ${sex} -i /mnt/data/${input_file}
    >>>

    runtime {
        docker: "stripy:v2.2"
        volumes: ["${ref_dir}:/mnt/ref", "${output_path}:/mnt/results", "${input_dir}:/mnt/data"]
    }
}
