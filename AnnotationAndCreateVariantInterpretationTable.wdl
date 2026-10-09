version 1.0

# Subworkflows
import "./AnnotationPipeline.wdl" as Annotation
import "./CreateInterpretationTable.wdl" as CreateInterpretationTable

workflow AnnotateAndTable {
  input {
    File input_variant_string

    String sample_basename

    File chromosome_list
    
    File gnomAD_vcf
    File gnomAD_vcf_index

    File gnomADexomes_vcf
    File gnomADexomes_vcf_index

    File SLOpopulation_vcf
    File SLOpopulation_vcf_index

    File TopMed_vcf
    File TopMed_vcf_index

    File Regeneron_vcf
    File Regeneron_vcf_index

    File GnomAD4_exomes_vcf
    File GnomAD4_exomes_vcf_index


    File GnomAD4_genomes_vcf
    File GnomAD4_genomes_vcf_index

    File blacklisted_regions_bed
    File blacklisted_regions_bed_index

    File GERP_bed
    File GERP_bed_index

    File ClinVar_vcf
    File ClinVar_vcf_index

    File SpliceAI
    File SpliceAI_index

    File dbscSNV
    File dbscSNV_index

    File HPO
    File HPO_index
    File OMIM
    File OMIM_index
    File gnomadConstraints
    File gnomadConstraints_index
    File CGD
    File CGD_index
    File bcftools_annotation_header

    File pext_bed
    File pext_bed_index

    File gnomAD_pLi_bed
    File gnomAD_pLi_bed_index

    File gnomAD_misz_bed
    File gnomAD_misz_bed_index

    File metaDome_bed
    File metaDome_bed_index

    File reference_fa
    File reference_fai
    File reference_dict

    File dbNSFP
    File dbNSFP_index

    String? targetRegions
    
    #String bgzip_docker = "dockerbiotools/bcftools:latest"
    String bcftools_docker = "dceoy/bcftools:latest"
    #String bcftools_docker = "biocontainers/bcftools:v1.9-1-deb_cv1"
    String SnpEff_docker = "alesmaver/snpeff_v50:latest"
    String gatk_docker = "broadinstitute/gatk:latest"
    String gatk_path = "/gatk/gatk"
    String vcfanno_docker = "clinicalgenomics/vcfanno:0.3.2"
    String R_docker = "alesmaver/r-base"
  }

  String final_output_xlsx = sample_basename + ".IMPORTVARIANT.FinalReportNew.xlsx"
  String final_output_vcf = sample_basename + ".IMPORTVARIANT.annotated.vcf.gz"

  call CreateVCFfromString {
      input:
        input_variant_string = input_variant_string,
        sample_basename = sample_basename,

        docker = R_docker
    }

  call GetVariantRegions {
      input:
        input_variant_string = input_variant_string,
        separator = ";"
    }

  call Annotation.AnnotateVCF as AnnotateVCF {
    input:
      sample_basename = sample_basename,

      input_vcf = CreateVCFfromString.output_vcf,
      chromosome_list = chromosome_list,
      
      gnomAD_vcf = gnomAD_vcf,
      gnomAD_vcf_index = gnomAD_vcf_index,

      gnomADexomes_vcf = gnomADexomes_vcf,
      gnomADexomes_vcf_index = gnomADexomes_vcf_index,

      SLOpopulation_vcf = SLOpopulation_vcf,
      SLOpopulation_vcf_index = SLOpopulation_vcf_index,

      TopMed_vcf = TopMed_vcf,
      TopMed_vcf_index = TopMed_vcf_index,

      GnomAD4_exomes_vcf = GnomAD4_exomes_vcf,
      GnomAD4_exomes_vcf_index = GnomAD4_exomes_vcf_index,

      GnomAD4_genomes_vcf = GnomAD4_genomes_vcf,
      GnomAD4_genomes_vcf_index = GnomAD4_genomes_vcf_index,

      blacklisted_regions_bed = blacklisted_regions_bed,
      blacklisted_regions_bed_index = blacklisted_regions_bed_index,

      GERP_bed = GERP_bed,
      GERP_bed_index = GERP_bed_index,

      Regeneron_vcf = Regeneron_vcf,
      Regeneron_vcf_index = Regeneron_vcf_index,

      SLOpopulation_vcf = SLOpopulation_vcf,
      SLOpopulation_vcf_index = SLOpopulation_vcf_index, 

      ClinVar_vcf = ClinVar_vcf,
      ClinVar_vcf_index = ClinVar_vcf_index,

      SpliceAI = SpliceAI,
      SpliceAI_index = SpliceAI_index,

      dbscSNV = dbscSNV,
      dbscSNV_index = dbscSNV_index,

      HPO = HPO,
      HPO_index = HPO_index,
      OMIM = OMIM,
      OMIM_index = OMIM_index,
      gnomadConstraints = gnomadConstraints,
      gnomadConstraints_index = gnomadConstraints_index,
      CGD = CGD,
      CGD_index = CGD_index,
      bcftools_annotation_header = bcftools_annotation_header,

      pext_bed = pext_bed,
      pext_bed_index = pext_bed_index,

      gnomAD_pLi_bed = gnomAD_pLi_bed,
      gnomAD_pLi_bed_index = gnomAD_pLi_bed_index,

      gnomAD_misz_bed = gnomAD_misz_bed,
      gnomAD_misz_bed_index = gnomAD_misz_bed_index,

      metaDome_bed = metaDome_bed,
      metaDome_bed_index = metaDome_bed_index,

      fasta_reference = reference_fa,
      fasta_reference_index = reference_fai,
      fasta_reference_dict = reference_dict,

      dbNSFP = dbNSFP,
      dbNSFP_index = dbNSFP_index,

      targetRegions = select_first([GetVariantRegions.regions, targetRegions]),

      #bcftools_docker = bcftools_docker,
      #SnpEff_docker = SnpEff_docker,
      gatk_docker = gatk_docker,
      gatk_path = gatk_path,
      vcfanno_docker = vcfanno_docker
  }

  # Get snpEff and dbNSFP annotations
  call GetGeneListFromVCF {
    input:
      input_vcf = AnnotateVCF.output_vcf,
      docker = SnpEff_docker
    }

  call CreateInterpretationTable.CreateInterpretationTable as CreateInterpretationTable {
    input:
      input_vcf = AnnotateVCF.output_vcf,
      input_vcf_index = AnnotateVCF.output_vcf_index,
      panel_gene_list = GetGeneListFromVCF.panel_gene_list
  }

  call RenameFile {
    input:
      input_file = CreateInterpretationTable.XLSX_OUTPUT,
      new_name = final_output_xlsx
  }

  # Distinct names so the annotated VCF does not clash with the main pipeline's {sample}.annotated.vcf.gz
  call RenameFile as RenameAnnotatedVCF {
    input:
      input_file = AnnotateVCF.output_vcf,
      new_name = final_output_vcf
  }

  call RenameFile as RenameAnnotatedVCFIndex {
    input:
      input_file = AnnotateVCF.output_vcf_index,
      new_name = final_output_vcf + ".tbi"
  }

  # Outputs that will be retained when execution is complete
  output {
    File XLSX_OUTPUT = RenameFile.renamed_file
    File ANNOTATED_VCF = RenameAnnotatedVCF.renamed_file
    File ANNOTATED_VCF_INDEX = RenameAnnotatedVCFIndex.renamed_file
  }
}

task GetGeneListFromVCF {
	input { 
		File input_vcf

		String docker
	}

	command {
		SNPSIFT_EXTRACTFIELDS='java -jar /home/biodocker/bin/snpEff/SnpSift.jar extractFields -s "," - "ANN[*].GENE" | uniq '
		zcat ~{input_vcf} | eval $SNPSIFT_EXTRACTFIELDS | tail -n +2 | paste -sd "," - 
	}

	runtime {
	    docker: docker
	    requested_memory_mb_per_core: 3000
	    cpu: 1
	    runtime_minutes: 180
	}
	output {
		String panel_gene_list = read_string(stdout())
	}
}

# Output VCF from string
task CreateVCFfromString {
  input {
    # Command parameters
    String input_variant_string
    String CreateVCFfromString_Rscript = "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/R_scripts/SCRIPTS_createVCFfromVariantString.R"
    String sample_basename

    # Runtime parameters
    String docker
  }

  command <<<
  wget -t 1 -T 20 ~{CreateVCFfromString_Rscript}
  # Repeat in case the proxy defined in the docker image would case problems accessing the GitHub repo
  unset https_proxy
  wget -t 1 -T 20 ~{CreateVCFfromString_Rscript}

  Rscript SCRIPTS_createVCFfromVariantString.R --variant="~{input_variant_string}" --sample="~{sample_basename}"

  # Sort the VCF file to prevent issues with downstream commands
  mv ~{sample_basename}.vcf ~{sample_basename}.unsorted.vcf
  cat ~{sample_basename}.unsorted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > ~{sample_basename}.vcf
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 1000
    cpu: 1
    runtime_minutes: 30
  }
  output {
    File output_vcf = "~{sample_basename}.vcf"
  }
}

task GetVariantRegions {
  input {
    String input_variant_string
    String separator
  }
  command <<<
    echo "~{input_variant_string}" | tr "~{separator}" \\n | awk -F' ' '{print $1":"$2-100"-"$2+100}' OFS='\t' | tr \\n ";" | sed 's/\;$//'
  >>>
  runtime {
    docker:"biocontainers/bcftools:v1.9-1-deb_cv1"
    requested_memory_mb_per_core: 500
    cpu: 1
    runtime_minutes: 5
  }
  output {
    Array[String] value = read_lines(stdout())
    String regions = value[0]
  }
}

task RenameFile {
  input {
    File input_file
    String new_name
  }

  command {
    mv ${input_file} ${new_name}
  }

  output {
    File renamed_file = "${new_name}"
  }
}
