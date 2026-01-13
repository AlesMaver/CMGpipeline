version 1.0

task RunCommands {

  command <<<
    wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/gnomAD/pLI/pliByGene.bb
    wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
    chmod +x ./bigBedToBed
    ./bigBedToBed pliByGene.bb pliByGene.bed12
    cat pliByGene.bed12 | awk 'BEGIN{OFS=FS="\t"} {temp=$4; $4=$15; $15=temp; print}' | bedtools bed12tobed6 | cut -f 1-4 > pliByGene.bed
    bgzip pliByGene.bed
    tabix -p bed pliByGene.bed.gz
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