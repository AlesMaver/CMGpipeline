version 1.0

workflow CreateGenesBed {
    call DownloadAndPrepareBed
}

task DownloadAndPrepareBed {
    input {
        File gencode_gff3
    }

    command <<<
        cat ~{gencode_gff3} \
            | gunzip --stdout - \
            | awk '$3 == "gene"' \
            | convert2bed -i gff --attribute-key="gene_name" - \
            > genes.bed
    >>>

    output {
        File genes_bed = "genes.bed"
    }

    runtime {
        docker: "dceoy/bedops"
        cpu: 1
        memory: "4G"
    }

}

task FilterGenesBED {
    input {
        File genes_bed
        Array[String] filter_genes
    }

    command <<<
        cat ~{genes_bed} | awk 'BEGIN{FS="\t"} FNR==NR{genes[$1]; next} $4 in genes' ~{write_lines(filter_genes)} - > filtered_genes.bed
    >>>

    output {
        File filtered_genes_bed = "filtered_genes.bed"
    }

    runtime {
        docker: "dceoy/bedops"
        cpu: 1
        memory: "4G"
    }

}
