version 1.0

task RunCommands {

  command <<<
    wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/gnomAD/missense/missenseConstrained.bb
    wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
    chmod +x ./bigBedToBed
    ./bigBedToBed missenseConstrained.bb missenseConstrained.bed
    cat missenseConstrained.bed | awk 'BEGIN{OFS=FS="\t"} {temp=$4; $4=$16; $16=temp; print}' | bedtools bed12tobed6 | cut -f 1-4 > missenseConstrained_oe.bed
    cat missenseConstrained.bed | awk 'BEGIN{OFS=FS="\t"} {temp=$4; $4=$17; $17=temp; print}' | bedtools bed12tobed6 | cut -f 1-4 > missenseConstrained_chisq.bed
    bgzip missenseConstrained_oe.bed
    tabix -p bed missenseConstrained_oe.bed.gz
    bgzip missenseConstrained_chisq.bed
    tabix -p bed missenseConstrained_chisq.bed.gz
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