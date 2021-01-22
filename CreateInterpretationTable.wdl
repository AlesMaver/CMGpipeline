version 1.0
## Copyright CMG@KIGM, Ales Maver

# WORKFLOW DEFINITION 
workflow CreateInterpretationTable {
  input {
    #File input_vcf
    #String? panel_gene_list

    File input_vcf
    File input_vcf_index

    Array[File]? relative_vcfs
    Array[File]? relative_vcf_indexes

    String? panel_gene_list
    File? mitoResults_txt

    String SnpEff_docker = "alesmaver/snpeff_v43:latest"
    String R_docker = "alesmaver/r-base"
  }  

  #File makeXSLSXoutputs_Rscript = "/home/ales/FastqToVCFPipeline/makeXLSXoutputs.R"
  String GenerateXLSXscriptUrl = "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/makeXLSXoutputs.R"
  String sample_basename = basename(input_vcf, ".annotated.vcf.gz")

  if( defined(relative_vcfs) ) {
    call MergeMultipleSampleVCFs {
    input:
      input_vcf = input_vcf,
      input_vcf_index = input_vcf_index,
      relative_vcfs = relative_vcfs,
      relative_vcf_indexes = relative_vcf_indexes,
      sample_basename = sample_basename
    }
  }

  # Get snpEff and dbNSFP annotations
    call GenerateVariantTable {
    input:
      single_vcf = input_vcf,
      merged_vcf = MergeMultipleSampleVCFs.output_vcf,
      merged_sample_list = MergeMultipleSampleVCFs.merged_sample_list,
      panel_gene_list = panel_gene_list,
      sample_basename = sample_basename,
      docker = SnpEff_docker
    }

  # Download R script for XLSX file generation
    call GetGenerateXLSXscript {
    input:
      GenerateXLSXscriptUrl = GenerateXLSXscriptUrl,
      timestamp = GenerateVariantTable.timestamp,
      docker = "davidsouthgate/alpine-bash-wget"
    }

  # Get snpEff and dbNSFP annotations
    call GenerateXLSX {
    input:
      RARE_FUNCTIONAL = GenerateVariantTable.RARE_FUNCTIONAL,
      HET_DOMINANT = GenerateVariantTable.HET_DOMINANT,
      COMPHET_RECESSIVE = GenerateVariantTable.COMPHET_RECESSIVE,
      HOM_RECESSIVE = GenerateVariantTable.HOM_RECESSIVE,
      CLINVAR_PATHOGENIC = GenerateVariantTable.CLINVAR_PATHOGENIC,
      CLINVAR_FILTERED = GenerateVariantTable.CLINVAR_FILTERED,
      CLINVAR_ALL = GenerateVariantTable.CLINVAR_ALL,
      PANEL_FILTERED = GenerateVariantTable.PANEL_FILTERED,
      PANEL_ALL = GenerateVariantTable.PANEL_ALL,
      mitoResults_txt = mitoResults_txt,

      sample_basename = sample_basename,

      makeXSLSXoutputs_Rscript = GetGenerateXLSXscript.makeXSLSXoutputs_Rscript,

      docker = R_docker
  }
 
  # Outputs that will be retained when execution is complete
  output {
    File XLSX_OUTPUT = GenerateXLSX.XLSX_OUTPUT
  }
}

##################
# TASK DEFINITIONS
##################
# Generate table of variants for interpretation
task MergeMultipleSampleVCFs {
  input {
    File input_vcf
    File input_vcf_index
    Array[File]? relative_vcfs
    Array[File]? relative_vcf_indexes
    String sample_basename
  }
  
  command <<<
  set -e
    bcftools merge -m none ~{input_vcf} ~{sep=' ' relative_vcfs} -Oz -o ~{sample_basename}.merged.vcf.gz
    bcftools query -l ~{sample_basename}.merged.vcf.gz > ~{sample_basename}.merged_sample_list.txt
  >>>

  runtime {
    docker: "biocontainers/bcftools:v1.9-1-deb_cv1"
    maxRetries: 3
    requested_memory_mb_per_core: 5000
    cpu: 1
    runtime_minutes: 90
  }
  output {
    File output_vcf = "~{sample_basename}.merged.vcf.gz"
    File merged_sample_list = "~{sample_basename}.merged_sample_list.txt"
  }
}

task GenerateVariantTable {
  input {
    # Command parameters
    File single_vcf
    File? merged_vcf
    File? merged_sample_list
    
    String sample_basename
    String? panel_gene_list
    # Runtime parameters
    String docker
  }

  File input_vcf = select_first([merged_vcf, single_vcf])

  command <<<
    set -e
    if [ -f ~{if defined(merged_sample_list) then merged_sample_list else "merged.sample_list"} ]; then
      GENOTYPE_EXTRACTFIELDS=$(cat ~{merged_sample_list} | awk '{print "GEN["$1"].GT GEN["$1"].AD GEN["$1"].DP GEN["$1"].GQ"}' | tr "\\n" " ")
    else
      GENOTYPE_EXTRACTFIELDS=$(echo ~{sample_basename} | awk '{print "GEN["$1"].GT GEN["$1"].AD GEN["$1"].DP GEN["$1"].GQ"}' | tr "\\n" " ")
    fi

    #SNPSIFT_EXTRACTFIELDS='/home/biodocker/bin/snpEff/scripts/vcfEffOnePerLine.pl | java -jar /home/biodocker/bin/snpEff/SnpSift.jar extractFields  - CHROM POS REF ALT QUAL '$GENOTYPE_EXTRACTFIELDS' "ANN[*].GENE" "Disease_name" "Categorization" Inheritance Age HPO "OMIM" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].RANK" "ANN[*].HGVS_P" "ANN[*].IMPACT" "ANN[*].EFFECT"  SLOpopulation.AC_Het SLOpopulation.AC_Hom SLOpopulation.AC_Hemi gnomAD.AC gnomAD.AF gnomAD.nhomalt gnomADexomes.AC gnomADexomes.AF gnomADexomes.nhomalt clinvar.CLNSIG clinvar.CLNDN clinvar.CLNHGVS clinvar.CLNSIGCONF clinvar.CLNSIGINCL dbNSFP_REVEL_rankscore  dbNSFP_MetaSVM_pred dbNSFP_CADD_phred  dbNSFP_DANN_rankscore dbNSFP_SIFT_pred  dbNSFP_SIFT4G_pred  dbNSFP_Polyphen2_HDIV_pred  dbNSFP_MutationTaster_pred dbNSFP_PrimateAI_pred dbNSFP_Polyphen2_HDIV_score SpliceAI.SpliceAI dbscSNV.ada_score dbscSNV.rf_score dbNSFP_GERP___NR  dbNSFP_GERP___RS dbNSFP_Interpro_domain pLI oe_mis pRec "LOF[*].GENE" "LOF[*].GENEID" "LOF[*].NUMTR" "LOF[*].PERC" "NMD[*].GENE" "NMD[*].GENEID" "NMD[*].NUMTR" "NMD[*].PERC" | uniq | grep -v "structural_interaction_variant"'

    SNPSIFT_EXTRACTFIELDS='java -jar /home/biodocker/bin/snpEff/SnpSift.jar extractFields -s "," - CHROM POS REF ALT QUAL '$GENOTYPE_EXTRACTFIELDS' "ANN[*].GENE" "Disease_name" "Categorization" Inheritance Age HPO "OMIM" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].RANK" "ANN[*].HGVS_P" "ANN[*].IMPACT" "ANN[*].EFFECT"  SLOpopulation.AC_Het SLOpopulation.AC_Hom SLOpopulation.AC_Hemi gnomAD.AC gnomAD.AF gnomAD.nhomalt gnomADexomes.AC gnomADexomes.AF gnomADexomes.nhomalt clinvar.ID clinvar.CLNSIG clinvar.CLNDN clinvar.CLNHGVS clinvar.CLNSIGCONF clinvar.CLNSIGINCL dbNSFP_REVEL_rankscore  dbNSFP_MetaSVM_pred dbNSFP_CADD_phred  dbNSFP_DANN_rankscore dbNSFP_SIFT_pred  dbNSFP_SIFT4G_pred  dbNSFP_Polyphen2_HDIV_pred  dbNSFP_MutationTaster_pred dbNSFP_PrimateAI_pred dbNSFP_Polyphen2_HDIV_score SpliceAI.SpliceAI dbscSNV.ada_score dbscSNV.rf_score dbNSFP_GERP___NR  dbNSFP_GERP___RS dbNSFP_Interpro_domain pLI oe_mis pRec "LOF[*].GENE" "LOF[*].GENEID" "LOF[*].NUMTR" "LOF[*].PERC" "NMD[*].GENE" "NMD[*].GENEID" "NMD[*].NUMTR" "NMD[*].PERC" | uniq '

    # Optional fields available in the VCF files, consider adding them later
    # dbNSFP_GTEx_V7_gene dbNSFP_GTEx_V7_tissue dbNSFP_Geuvadis_eQTL_target_gene dbNSFP_Aloft_Fraction_transcripts_affected  dbNSFP_Aloft_prob_Tolerant  dbNSFP_Aloft_prob_Recessive  dbNSFP_Aloft_prob_Dominant  dbNSFP_Aloft_pred dbNSFP_Aloft_Confidence 
    zcat ~{input_vcf} | java -jar /home/biodocker/bin/snpEff/SnpSift.jar filter "((( gnomAD.AF < 0.01 ) | !( exists gnomAD.AC )) & (( gnomAD.nhomalt < 5 ) | !( exists gnomAD.nhomalt )) & (( SLOpopulation.AC_Het < 20 ) | !( exists SLOpopulation.AC_Het )) & (!( exists SLOpopulation.AC_Hom ) | ( SLOpopulation.AC_Hom <= 3 ))) & (( ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'HIGH') | (dbscSNV.ada_score > 0.5) | (dbscSNV.rf_score > 0.5)) & isVariant( GEN[0].GT )" | eval $SNPSIFT_EXTRACTFIELDS > ~{sample_basename}.RARE_FUNCTIONAL.tab

    zcat ~{input_vcf} | java -jar /home/biodocker/bin/snpEff/SnpSift.jar filter "((( gnomADexomes.AC < 10 ) | !( exists gnomADexomes.AC )) & (( gnomAD.AC < 10 ) | !( exists gnomAD.AC )) & (( SLOpopulation.AC_Het < 6 ) | !( exists SLOpopulation.AC_Het ))) & (( ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'HIGH') | (dbscSNV.ada_score > 0.5) | (dbscSNV.rf_score > 0.5)) & (( Inheritance =~ '.*AD.*' ) | ( HPO =~ '.*Autosomal_dominant.*' ) | ( HPO  =~  '.*X-linked_dominant.*' )) & (!(clinvar.CLNSIG =~ '.*Likely_benign.*')) & (!(clinvar.CLNSIG =~ '.*Benign.*')) & isHet(GEN[0].GT)" | eval $SNPSIFT_EXTRACTFIELDS > ~{sample_basename}.HET_DOMINANT.tab

    zcat ~{input_vcf} | java -jar /home/biodocker/bin/snpEff/SnpSift.jar filter "((( gnomAD.AF < 0.05 ) | !( exists gnomAD.AC )) & (( gnomAD.nhomalt < 10 ) | !( exists gnomAD.nhomalt )) & (( SLOpopulation.AC_Het < 100 ) | !( exists SLOpopulation.AC_Het )) & (!( exists SLOpopulation.AC_Hom ) | ( SLOpopulation.AC_Hom <= 6 ))) & (( ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'HIGH') | (dbscSNV.ada_score > 0.5) | (dbscSNV.rf_score > 0.5)) & (( Inheritance =~ '.*AR.*' ) | ( Inheritance =~ '.*XL.*' ) | ( HPO =~ '.*Autosomal_recessive.*' ) | ( HPO  =~  '.*X-linked_recessive.*' )) & (!(clinvar.CLNSIG =~ '.*Likely_benign.*')) & (!(clinvar.CLNSIG =~ '.*Benign.*')) & !(isHom(GEN[0].GT)) & isVariant( GEN[0].GT )" | eval $SNPSIFT_EXTRACTFIELDS >  ~{sample_basename}.COMPHET_RECESSIVE.tab

    zcat ~{input_vcf} | java -jar /home/biodocker/bin/snpEff/SnpSift.jar filter "((( gnomADexomes.AF < 0.05 ) | !( exists gnomADexomes.AC )) & (( gnomADexomes.nhomalt < 10 ) | !( exists gnomADexomes.nhomalt )) & (( gnomAD.AF < 0.05 ) | !( exists gnomAD.AC )) & (( gnomAD.nhomalt < 10 ) | !( exists gnomAD.nhomalt )) & (( SLOpopulation.AC_Het < 100 ) | !( exists SLOpopulation.AC_Het )) & (!( exists SLOpopulation.AC_Hom ) | ( SLOpopulation.AC_Hom <= 6 ))) & (( ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'HIGH') | (dbscSNV.ada_score > 0.5) | (dbscSNV.rf_score > 0.5) | (( gnomADexomes.AF < 0.001 ) | !( exists gnomADexomes.AC ))  ) & (( Inheritance =~ '.*AR.*' ) | ( Inheritance =~ '.*XL.*' ) | ( HPO =~ '.*Autosomal_recessive.*' ) | ( HPO  =~  '.*X-linked_recessive.*' )) & (!(clinvar.CLNSIG =~ '.*Likely_benign.*')) & (!(clinvar.CLNSIG =~ '.*Benign.*')) & isHom(GEN[0].GT) & isVariant( GEN[0].GT )" | eval $SNPSIFT_EXTRACTFIELDS > ~{sample_basename}.HOM_RECESSIVE.tab

    zcat ~{input_vcf} | java -jar /home/biodocker/bin/snpEff/SnpSift.jar filter "((( gnomAD.AF < 0.05 ) | !( exists gnomAD.AC )) & (( gnomAD.nhomalt < 10 ) | !( exists gnomAD.nhomalt )) & (( SLOpopulation.AC_Het < 100 ) | !( exists SLOpopulation.AC_Het )) & (!( exists SLOpopulation.AC_Hom ) | ( SLOpopulation.AC_Hom <= 6 ))) & (( ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'HIGH') | (dbscSNV.ada_score > 0.5) | (dbscSNV.rf_score > 0.5)) & ((clinvar.CLNSIGCONF =~ '.*Pathogenic.*') | (clinvar.CLNSIGCONF =~ '.*Likely_pathogenic.*') | (clinvar.CLNSIG =~ '.*Likely_pathogenic.*') | (clinvar.CLNSIG =~ '.*Pathogenic.*')) & isVariant( GEN[0].GT )" | eval $SNPSIFT_EXTRACTFIELDS > ~{sample_basename}.CLINVAR_PATHOGENIC.tab

    zcat ~{input_vcf} | java -jar /home/biodocker/bin/snpEff/SnpSift.jar filter "((( gnomAD.AF < 0.05 ) | !( exists gnomAD.AC )) & (( gnomAD.nhomalt < 10 ) | !( exists gnomAD.nhomalt )) & (( SLOpopulation.AC_Het < 100 ) | !( exists SLOpopulation.AC_Het )) & (!( exists SLOpopulation.AC_Hom ) | ( SLOpopulation.AC_Hom <= 6 ))) & (( ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'HIGH') | (dbscSNV.ada_score > 0.5) | (dbscSNV.rf_score > 0.5)) & (exists clinvar.CLNSIG) & isVariant( GEN[0].GT )" | eval $SNPSIFT_EXTRACTFIELDS > ~{sample_basename}.CLINVAR_FILTERED.tab

    zcat ~{input_vcf} | java -jar /home/biodocker/bin/snpEff/SnpSift.jar filter "(exists clinvar.CLNSIG) & !((clinvar.CLNSIG =~ '.*Likely_benign.*') | (clinvar.CLNSIG =~ '.*Benign.*')) & isVariant( GEN[0].GT )" | eval $SNPSIFT_EXTRACTFIELDS > ~{sample_basename}.CLINVAR_ALL.tab

    ~{if defined(panel_gene_list) then "echo " + panel_gene_list + " | awk '{gsub(/\,/,\"\\n\",$0);print $0}' > panel_gene_list.txt" else ""}

    if [ -f panel_gene_list.txt ]; then
      echo "Generating panel filtered variant results..."
      cat panel_gene_list.txt

      zcat ~{input_vcf} | java -jar /home/biodocker/bin/snpEff/SnpSift.jar filter -s panel_gene_list.txt "(     (    (   (  ( gnomADexomes.AF < 0.05 ) | !( exists gnomADexomes.AC )  ) &  (  ( gnomADexomes.nhomalt < 10 ) | !( exists gnomADexomes.nhomalt )  ) & (  ( gnomAD.AF < 0.05 ) | !( exists gnomAD.AC )  ) & (  ( gnomAD.nhomalt < 10 ) | !( exists gnomAD.nhomalt )  ) & (  ( SLOpopulation.AC_Het < 100 ) | !( exists SLOpopulation.AC_Het )  ) & (  !( exists SLOpopulation.AC_Hom ) | ( SLOpopulation.AC_Hom <= 6 )  )   ) & (   ( ANN[*].IMPACT has 'MODERATE') | (ANN[*].IMPACT has 'HIGH') | (dbscSNV.ada_score > 0.5) | (dbscSNV.rf_score > 0.5) ) &  (!(clinvar.CLNSIG =~ '.*Likely_benign.*')) & (!(clinvar.CLNSIG =~ '.*Benign.*'))    )    |    ( clinvar.CLNSIG =~ '.*Likely_pathogenic.*' | clinvar.CLNSIG =~ '.*Pathogenic.*' | clinvar.CLNSIGCONF =~ '.*Likely_pathogenic.*' | clinvar.CLNSIGCONF =~ '.*Pathogenic.*' )    |    ( dbscSNV.ada_score > 0.5 | dbscSNV.rf_score > 0.5 )    |    (   ((( gnomADexomes.AF < 0.001 ) | !( exists gnomADexomes.AC )) & (( gnomADexomes.nhomalt < 4 ) | !( exists gnomADexomes.nhomalt )) & (( gnomAD.AF < 0.001 ) | !( exists gnomAD.AC )) & (( gnomAD.nhomalt < 4 ) | !( exists gnomAD.nhomalt )) & (( SLOpopulation.AC_Het < 10 ) | !( exists SLOpopulation.AC_Het )) & (!( exists SLOpopulation.AC_Hom ) | ( SLOpopulation.AC_Hom <= 4 ))) & isHom(GEN[0].GT)  )    )    &    isVariant( GEN[0].GT )  &  GEN[0].DP>5  & (ANN[*].GENE in SET[0])" | eval $SNPSIFT_EXTRACTFIELDS > ~{sample_basename}.PANEL_FILTERED.tab

      zcat ~{input_vcf} | java -jar /home/biodocker/bin/snpEff/SnpSift.jar filter -s panel_gene_list.txt "(isVariant( GEN[0].GT )) & (ANN[*].GENE in SET[0])" | eval $SNPSIFT_EXTRACTFIELDS > ~{sample_basename}.PANEL_ALL.tab
    fi

    echo $(date +"%Y_%m_%d_%I_%M_%p") > timestamp
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 6000
    cpu: 1
    runtime_minutes: 180
  }
  output {
    File RARE_FUNCTIONAL = "~{sample_basename}.RARE_FUNCTIONAL.tab"
    File HET_DOMINANT = "~{sample_basename}.HET_DOMINANT.tab"
    File COMPHET_RECESSIVE = "~{sample_basename}.COMPHET_RECESSIVE.tab"
    File HOM_RECESSIVE = "~{sample_basename}.HOM_RECESSIVE.tab"
    File CLINVAR_PATHOGENIC = "~{sample_basename}.CLINVAR_PATHOGENIC.tab"
    File CLINVAR_FILTERED = "~{sample_basename}.CLINVAR_FILTERED.tab"
    File CLINVAR_ALL = "~{sample_basename}.CLINVAR_ALL.tab"
    File? PANEL_FILTERED = "~{sample_basename}.PANEL_FILTERED.tab"
    File? PANEL_ALL = "~{sample_basename}.PANEL_ALL.tab"
    File? timestamp = "timestamp"
  }
}

task GetGenerateXLSXscript {
  input {
    String GenerateXLSXscriptUrl
    File? timestamp

    # Runtime parameters
    String docker
  }

  command {
    set -e
    wget ~{GenerateXLSXscriptUrl}
    cat ~{timestamp}
  }

  runtime {
    docker: docker
    requested_memory_mb_per_core: 2000
    cpu: 1
    runtime_minutes: 20
  }

  output {
    File makeXSLSXoutputs_Rscript = "makeXLSXoutputs.R"
  }
}

# Generate table of variants for interpretation
task GenerateXLSX {
  input {
    # Command parameters
    File RARE_FUNCTIONAL
    File HET_DOMINANT
    File COMPHET_RECESSIVE
    File HOM_RECESSIVE
    File CLINVAR_PATHOGENIC
    File CLINVAR_FILTERED
    File CLINVAR_ALL
    File? PANEL_FILTERED
    File? PANEL_ALL
    File? mitoResults_txt
    String sample_basename

    File makeXSLSXoutputs_Rscript

    # Runtime parameters
    String docker
  }

  command {
  Rscript ~{makeXSLSXoutputs_Rscript} --sample_basename=~{sample_basename} --RARE_FUNCTIONAL=~{RARE_FUNCTIONAL} --HET_DOMINANT=~{HET_DOMINANT} --COMPHET_RECESSIVE=~{COMPHET_RECESSIVE} --HOM_RECESSIVE=~{HOM_RECESSIVE} --CLINVAR_PATHOGENIC=~{CLINVAR_PATHOGENIC} --CLINVAR_FILTERED=~{CLINVAR_FILTERED} --CLINVAR_ALL=~{CLINVAR_ALL} ~{if defined(PANEL_FILTERED) then " --PANEL_FILTERED " + PANEL_FILTERED else ""}  ~{if defined(mitoResults_txt) then " --MITOMAP " + mitoResults_txt else ""}  ~{if defined(PANEL_ALL) then " --PANEL_ALL " + PANEL_ALL else ""} --XLSX_OUTPUT=~{sample_basename}.FinalReportNew.xlsx     
  }
  runtime {
    docker: docker
    requested_memory_mb_per_core: 6000
    cpu: 1
    runtime_minutes: 60
  }
  output {
    File XLSX_OUTPUT = "~{sample_basename}.FinalReportNew.xlsx"

  }
}

