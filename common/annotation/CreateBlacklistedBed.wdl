version 1.0

task RunCommands {

  command <<<
    wget https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz
    gunzip hg19-blacklist.v2.bed.gz
    bgzip hg19-blacklist.v2.bed
    tabix -p bed hg19-blacklist.v2.bed.gz
  >>>
  output {
    File output_file = "hg19-blacklist.v2.bed.gz"
    File output_file_index = "hg19-blacklist.v2.bed.gz.tbi"
  }
  runtime {
    docker: "ubuntu:latest"
  }
}

workflow ExecuteCommands {
  call RunCommands {}
}