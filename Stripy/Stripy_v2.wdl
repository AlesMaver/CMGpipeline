version 1.0

workflow Stripy {
    input {
        String sample_basename
        File input_bam_or_cram
        File input_bam_or_cram_index
        File reference_fasta
        # sex: male / female (case sensitive)
        String sex
        String reference_genome_name = "hg19"
    }


    call run_stripy {
        input:
            sample_basename = sample_basename,
            reference = reference_fasta,
            genome = reference_genome_name,
            sex = if defined(sex) && (sex == "male" || sex == "female") then sex else "male",
            input_file = input_bam_or_cram,
            input_file_index = input_bam_or_cram_index
    }

    output {
        File? stripy_tsv  = run_stripy.stripy_tsv
        File? stripy_html = run_stripy.stripy_html
    }
}


task run_stripy {
    input {
        String sample_basename
        File reference
        String genome
        String sex
        File input_file
        File input_file_index
    }

    String base_name = basename(input_file)

    command <<<
        set -e

        echo ~{sample_basename}
        echo ~{sex}
        echo ' '
        echo "[ PREPARATION ] Downloading variant catalog JSON"
        #wget "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/ExpansionHunter_configuration/variant_catalog.json"
        unset https_proxy
        wget "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/ExpansionHunter_configuration/variant_catalog.json"

        echo "[ PREPARATION ] Preparing LOCI"
        loci=$(jq -r '[.[] | .LocusId] | join(",")' ./variant_catalog.json)

        echo "[ RUNNING ] Stri.py"
        # Constructing Docker run command (inside Docker already)
        echo ' '
        #echo /usr/local/bin/stripy-pipeline/batch.sh -o . -r ~{reference} -l "\"$loci\"" -g ~{genome} -s ~{sex} -i ~{input_file}
        #/usr/local/bin/stripy-pipeline/batch.sh -o . -r ~{reference} -l "\"$loci\"" -g ~{genome} -s ~{sex} -i ~{input_file}
        echo ' '
        echo python3 /usr/local/bin/stripy-pipeline/stri.py --input ~{input_file} --locus "\"$loci\"" --sex ~{sex} --genome ~{genome} --reference ~{reference} --output .
        python3 /usr/local/bin/stripy-pipeline/stri.py --input ~{input_file} --locus "\"$loci\"" --sex ~{sex} --genome ~{genome} --reference ~{reference} --output .
        mv ~{base_name}.html ~{sample_basename}.Stripy.html
        mv ~{base_name}.tsv  ~{sample_basename}.Stripy.tsv
        echo  ' '
        echo 'END.'
    >>>

    runtime {
        docker: "gbergant/stripy_prod:2.5"
        requested_memory_mb_per_core: 1000
        cpu: 4
        runtime_minutes: 30
    }

    output {
        File? stripy_tsv = "~{sample_basename}.Stripy.tsv"
        File? stripy_html = "~{sample_basename}.Stripy.html"
    }
}
