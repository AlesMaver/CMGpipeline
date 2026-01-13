version 1.0

task RunCommands {

  command <<<
    # This will download and prepare references for hg19 and GRCh38 assemblies
    exomiser_version=2309
    project_version=13.3.0

    wget https://data.monarchinitiative.org/exomiser/latest/exomiser-cli-${project_version}-distribution.zip

    wget https://data.monarchinitiative.org/exomiser/latest/${exomiser_version}_hg19.zip 
    wget https://data.monarchinitiative.org/exomiser/latest/${exomiser_version}_phenotype.zip
    wget https://zenodo.org/records/6576087/files/ReMM.v0.4.hg19.tsv.gz?download=1 -O ReMM.v0.4.hg19.tsv.gz
    wget https://zenodo.org/records/6576087/files/ReMM.v0.4.hg19.tsv.gz.tbi?download=1 -O ReMM.v0.4.hg19.tsv.gz.tbi

    wget https://data.monarchinitiative.org/exomiser/latest/${exomiser_version}_hg38.zip
    wget https://data.monarchinitiative.org/exomiser/latest/${exomiser_version}_phenotype.zip
    wget https://zenodo.org/records/6576087/files/ReMM.v0.4.hg38.tsv.gz?download=1 -O ReMM.v0.4.hg38.tsv.gz
    wget https://zenodo.org/records/6576087/files/ReMM.v0.4.hg38.tsv.gz.tbi?download=1 -O ReMM.v0.4.hg38.tsv.gz.tbi

    unzip exomiser-cli-${project_version}-distribution.zip
    unzip ${exomiser_version}_hg19.zip -d exomiser-cli-${project_version}/data
    unzip ${exomiser_version}_phenotype.zip -d exomiser-cli-${project_version}/data

    # Test the installation 
    # cd /cmg1scratch/cromwell/VariantAnnotationPipeline_databases/exomiser/exomiser-cli-13.3.0
    # java -Xmx4g -jar exomiser-cli-13.3.0.jar  --analysis ./examples/test-analysis-exome.yml

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
