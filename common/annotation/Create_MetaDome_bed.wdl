version 1.0

task RunCommands {

  command <<<
    wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/metaDome/metaDome.bb
    wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
    chmod +x ./bigBedToBed
    ./bigBedToBed metaDome.bb metaDome.bed12
    cat metaDome.bed12 | awk 'BEGIN{OFS=FS="\t"} {temp=$4; $4=$5; $5=temp; print}' | bedtools bed12tobed6 | cut -f 1-4 > metaDome.bed
    bgzip metaDome.bed
    tabix -p bed metaDome.bed.gz
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

