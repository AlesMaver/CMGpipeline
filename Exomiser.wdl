version development

# Define the input parameter for the directory
workflow Exomiser {
  input {
    Directory input_dir
    File input_vcf
    String hpo_ids
  }

  # Call the task to list files in the directory
  call RunTask { input: 
                    input_dir = input_dir, 
                    input_vcf = input_vcf,
                    hpo_ids = hpo_ids 
                    }

  output {
    File exomiser_output_vcf = RunTask.output_vcf
    File exomiser_output_vcf_index = RunTask.output_vcf_index
    File exomiser_output_json= RunTask.output_json
    File exomiser_output_html = RunTask.output_html
    File exomiser_output_small_html = RunTask.output_small_html
    File exomiser_output_genes_tsv = RunTask.output_genes_tsv
    File exomiser_output_variants_tsv = RunTask.output_variants_tsv
  }
}

# Define the task to run for each file in the directory
task RunTask {
  input {
    Directory input_dir
    File input_vcf
    String hpo_ids
  }

  String sample_basename = basename(input_vcf, ".vcf")

  command <<<
    ls ~{input_dir}  # You can replace this with any command you want to run on the file
    echo ~{input_dir}
    echo ~{hpo_ids}
    echo ~{sample_basename}

    # Update file locations
    sed -i 's|^exomiser\.data-directory=.*|exomiser.data-directory=~{input_dir}/data|' ~{input_dir}/application.properties
    
    # Download the yq directly if needed
    # curl -fsSL https://github.com/mikefarah/yq/releases/download/v4.40.3/yq_linux_amd64 -o yq && chmod +x yq
    ~{input_dir}/yq e '.analysis.vcf = "~{input_vcf}"' -i ~{input_dir}/examples/test-analysis-exome.yml
    ~{input_dir}/yq e '.analysis.vcf = "~{input_vcf}"' -i ~{input_dir}/examples/test-analysis-exome-small.yml
    ~{input_dir}/yq e '.outputOptions.outputFileName = "~{sample_basename}-exomiser-small"' -i ~{input_dir}/examples/test-analysis-exome-small.yml

    ~{input_dir}/yq e '.analysis.analysisMode = "PASS_ONLY"' -i ~{input_dir}/examples/test-analysis-exome.yml
    ~{input_dir}/yq e '.analysis.analysisMode = "PASS_ONLY"' -i ~{input_dir}/examples/test-analysis-exome-small.yml

    echo ~{input_dir}/yq e '.analysis.hpoIds = [~{hpo_ids}]' -i ~{input_dir}/examples/test-analysis-exome.yml
    ~{input_dir}/yq e '.analysis.hpoIds = [~{hpo_ids}]' -i ~{input_dir}/examples/test-analysis-exome.yml
    
    java -Xmx12g -XX:+UseG1GC -jar ~{input_dir}/exomiser-cli-13.3.0.jar --analysis ~{input_dir}/examples/test-analysis-exome.yml --spring.config.location=~{input_dir}/
    java -Xmx12g -XX:+UseG1GC -jar ~{input_dir}/exomiser-cli-13.3.0.jar --analysis ~{input_dir}/examples/test-analysis-exome-small.yml --spring.config.location=~{input_dir}/
    mv results/* ./
    
    # Clean up the large copied annotation dir for exomiser immediately as hard links are not available for directory inputs
    rm -rf ~{input_dir}
  >>>

  output {
    # You can define task outputs here if needed
    File output_vcf = "~{sample_basename}-exomiser.vcf.gz"
    File output_vcf_index = "~{sample_basename}-exomiser.vcf.gz.tbi"
    File output_json = "~{sample_basename}-exomiser.json"
    File output_html = "~{sample_basename}-exomiser.html"
    File output_small_html = "~{sample_basename}-exomiser-small.html"
    File output_genes_tsv = "~{sample_basename}-exomiser.genes.tsv"
    File output_variants_tsv = "~{sample_basename}-exomiser.variants.tsv"
  }

  runtime {
      docker: "openjdk:22-jdk"
      requested_memory_mb_per_core: 1000
      cpu: 30
      runtime_minutes: 240
  }

}
