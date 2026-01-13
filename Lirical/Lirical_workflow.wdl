version 1.0

workflow Lirical_workflow {
    input {
        String sample_basename
        String hpo_ids
        File hg19_variants_mv_db
        File? input_vcf
        Boolean lirical_small = false
    }

    String output_prefix =  sample_basename + if ( defined(input_vcf) ) then ".vcf.lirical" else ".lirical"

    if(!lirical_small) {
        call Lirical {
            input:
                sample_basename = sample_basename,
                hpo_ids = hpo_ids,
                hg19_variants_mv_db = hg19_variants_mv_db,
                input_vcf = input_vcf,
                output_prefix = output_prefix
        }
    }

    if(lirical_small) {
        call LiricalSmall {
            input:
                sample_basename = sample_basename,
                hpo_ids = hpo_ids,
                hg19_variants_mv_db = hg19_variants_mv_db,
                input_vcf = input_vcf,
                output_prefix = output_prefix
        }
    }

    output {
        File? lirical_tsv = Lirical.lirical_tsv
        File? lirical_json = Lirical.lirical_json
        File? lirical_html = Lirical.lirical_html
        File? small_lirical_html = LiricalSmall.small_lirical_html
    }
}

task Lirical {
    input {
        String sample_basename
        String hpo_ids
        File hg19_variants_mv_db
        File? input_vcf
        String output_prefix
    }

    command {
        set -e

        java -jar /lirical/lirical-cli-2.0.0-RC2/lirical-cli-2.0.0-RC2.jar prioritize \
            -d /lirical/data/ \
            -e19 ~{hg19_variants_mv_db} \
            --assembly hg19 \
            -p "~{hpo_ids}" \
            ~{"--vcf " + input_vcf} \
            -g --display-all-variants --ddndv -t 0 \
            --sample-id ~{sample_basename} \
            -o . \
            -f html -f json -f tsv

        mv lirical.tsv ~{output_prefix}.tsv
        mv lirical.json ~{output_prefix}.json
        mv lirical.html ~{output_prefix}.html
    }

    runtime {
        docker: "alesmaver/lirical"
    }

    output {
        File? lirical_tsv = output_prefix + ".tsv"
        File? lirical_json = output_prefix + ".json"
        File? lirical_html = output_prefix + ".html"
    }
}


task LiricalSmall {
    input {
        String sample_basename
        String hpo_ids
        File hg19_variants_mv_db
        File? input_vcf
        String output_prefix
    }

    command {
        set -e

        java -jar /lirical/lirical-cli-2.0.0-RC2/lirical-cli-2.0.0-RC2.jar prioritize \
            -d /lirical/data/ \
            -e19 ~{hg19_variants_mv_db} \
            --assembly hg19 \
            -p "~{hpo_ids}" \
            ~{"--vcf " + input_vcf} \
            --sample-id ~{sample_basename} \
            -o . \
            -x small_lirical \
            -f html

        mv small_lirical.html small_~{output_prefix}.html
    }

    runtime {
        docker: "alesmaver/lirical"
    }

    output {
        File? small_lirical_html = "small_" + output_prefix + ".html"
    }
}
