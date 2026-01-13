version 1.0

task RunCommands {

  command <<<
    wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/gnomAD/missense/missenseByGene.bb
    wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
    chmod +x ./bigBedToBed
    ./bigBedToBed missenseByGene.bb missenseByGene.bed12
    cat missenseByGene.bed12 | awk 'BEGIN{OFS=FS="\t"} {temp=$4; $4=$14; $14=temp; print}' | bedtools bed12tobed6 | cut -f 1-4 > missenseByGene.bed
    bgzip missenseByGene.bed
    tabix -p bed missenseByGene.bed.gz
  >>>
  output {
    File output_file = "missenseConstrained_chisq.bed.gz"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow ExecuteCommands {
  call RunCommands {}
}

