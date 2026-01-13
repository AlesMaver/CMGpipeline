version 1.0

task RunCommands {

  command <<<
    wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw
    wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig
    chmod +x ./bigWigToWig
    ./bigWigToWig All_hg19_RS.bw GERP_hg19.wig

    cat GERP_hg19.wig | wig2bed -d | bgzip > GERP_hg19.bed.gz
    zcat GERP_hg19.bed.gz | sort -k 1,1 -k2,2n | bgzip > GERP_hg19.sorted.bed.gz
    tabix -p bed GERP_hg19.sorted.bed.gz
  >>>
  output {
    File output_file = "GERP_hg19.bed.gz"
    File output_file_index = "GERP_hg19.bed.gz.tbi"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow ExecuteCommands {
  call RunCommands {}
}