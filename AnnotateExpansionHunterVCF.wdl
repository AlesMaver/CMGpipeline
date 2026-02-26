version 1.0

import "./exp_hunter.wdl" as ExpHunter

workflow AnnotateExpansionHunterVCF {
  input {
    String sample_id
    File expansion_hunter_vcf
    String expansion_hunter_docker = "gbergant/expansionhunter:latest"
    Boolean trgt = false
    File? custom_catalog_file
  }

  parameter_meta {
    sample_id: "sample name"
    expansion_hunter_vcf: "ExpansionHunter VCF file to annotate"
    expansion_hunter_docker: "expansion hunter docker including annotation software (stranger)"
    trgt: "Set to true if file was produced with TRGT (adds -t flag to stranger)"
    custom_catalog_file: "Optional custom catalog file for TRGT or other tools (if not provided, uses default ExpansionHunter catalog)"
  }

  meta {
      author: "Aleš Maver"
      email: "cmg.kimg@kclj.si"
      description: "Simple workflow to annotate an existing ExpansionHunter VCF file using stranger"
  }

  call ExpHunter.AnnotateExpansionHunter {
      input:
        sample_id = sample_id,
        expansion_hunter_docker = expansion_hunter_docker,
        expansion_hunter_vcf = expansion_hunter_vcf,
        trgt = trgt,
        custom_catalog_file = custom_catalog_file
    }

  output {
    File expansion_hunter_vcf_annotated = AnnotateExpansionHunter.expansion_hunter_vcf_annotated
  }

}
