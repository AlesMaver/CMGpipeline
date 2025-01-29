version 1.0
## Copyright CMG@KIGM, Ales Maver

# WORKFLOW DEFINITION 
workflow AnnotateVCF {
  input {
    File input_vcf
    File chromosome_list
    
    File gnomAD_vcf
    File gnomAD_vcf_index

    File gnomADexomes_vcf
    File gnomADexomes_vcf_index

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

    File SLOpopulation_vcf
    File SLOpopulation_vcf_index

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

    File gnomAD_pLi_bed
    File gnomAD_pLi_bed_index

    File gnomAD_misz_bed
    File gnomAD_misz_bed_index

    File metaDome_bed
    File metaDome_bed_index

    File pext_bed
    File pext_bed_index

    File fasta_reference
    File fasta_reference_index
    File fasta_reference_dict

    File dbNSFP
    File dbNSFP_index

    String? targetRegions

    #String bgzip_docker = "dockerbiotools/bcftools:latest"
    ##String bcftools_docker = "dceoy/bcftools:latest"
    String bcftools_docker = "dceoy/bcftools"
    #String bcftools_docker = "biocontainers/bcftools:v1.9-1-deb_cv1"
    String SnpEff_docker = "alesmaver/snpeff_v50:latest"
    String gatk_docker = "broadinstitute/gatk:latest"
    String gatk_path = "/gatk/gatk"
    String vcfanno_docker = "clinicalgenomics/vcfanno:0.3.2"

    String sample_basename
    String output_filename = sample_basename + ".annotated.vcf"

  }  

  Array[String] chromosomes = read_lines(chromosome_list)

  #String sample_basename = basename(input_vcf, ".vcf")
  #String output_filename = sample_basename + ".annotated.vcf"
  #String output_filename_gz = sample_basename + ".annotated.vcf.gz"
  String output_filename_gz = output_filename + ".gz"

  # Workflow calls begin here
  #call LeftAlignAndTrimVariants {
  #  input:
  #    input_vcf = input_vcf,
  #    sample_basename = sample_basename,
  #    fasta_reference = fasta_reference,
  #    fasta_reference_index = fasta_reference_index,
  #    fasta_reference_dict = fasta_reference_dict,
  #    gatk_path = gatk_path,
  #    docker = gatk_docker
  #}

  call NormalizeVCF {
    input:
      input_vcf = input_vcf,
      sample_basename = sample_basename,
      fasta_reference = fasta_reference,
      docker = bcftools_docker
  }

  call CompressAndIndexVCF as CompressAndIndexVCF1 {
    input:
      input_vcf = NormalizeVCF.output_vcf,
      sample_basename = sample_basename,
      docker = bcftools_docker
  }

  if( defined(targetRegions) ) {
    call StringToArray {
      input:
        input_string = select_first([targetRegions, ""]),
        separator = ";"
    }
  }

  # Call variants in parallel over grouped calling intervals
  scatter (chromosome in select_first([StringToArray.values, chromosomes]) ) {
    call bcftoolsAnnotate {
      input:
        input_vcf = CompressAndIndexVCF1.output_vcfgz,
        input_vcf_index = CompressAndIndexVCF1.output_vcfgz_index,
        sample_basename=sample_basename,

        HPO = HPO,
        HPO_index = HPO_index,
        OMIM = OMIM,
        OMIM_index = OMIM_index,
        gnomadConstraints = gnomadConstraints,
        gnomadConstraints_index = gnomadConstraints_index,
        CGD = CGD,
        CGD_index = CGD_index,

        bcftools_annotation_header = bcftools_annotation_header,
        
        chromosome = chromosome,
        output_filename = sample_basename + ".annotated.vcf.gz",
        docker = bcftools_docker
    }

    call VCFANNO {
      input:
        input_vcf = bcftoolsAnnotate.output_vcf,
        input_vcf_index = bcftoolsAnnotate.output_vcf_index,
        sample_basename = sample_basename,

        gnomAD_vcf = gnomAD_vcf,
        gnomAD_vcf_index = gnomAD_vcf_index,

        gnomADexomes_vcf = gnomADexomes_vcf,
        gnomADexomes_vcf_index = gnomADexomes_vcf_index,

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

        pext_bed = pext_bed,
        pext_bed_index = pext_bed_index,

        gnomAD_pLi_bed = gnomAD_pLi_bed,
        gnomAD_pLi_bed_index = gnomAD_pLi_bed_index,

        gnomAD_misz_bed = gnomAD_misz_bed,
        gnomAD_misz_bed_index = gnomAD_misz_bed_index,

        metaDome_bed = metaDome_bed,
        metaDome_bed_index = metaDome_bed_index,

        docker = vcfanno_docker
      }

   
    # Get snpEff and dbNSFP annotations
      call runSnpEff {
      input:
        input_vcf = VCFANNO.output_vcfgz,
        sample_basename = sample_basename,
        dbNSFP = dbNSFP,
        dbNSFP_index = dbNSFP_index,
        docker = SnpEff_docker
    }

    call CompressAndIndexVCF as CompressAndIndexVCF2 {
    input:
      input_vcf = runSnpEff.output_vcf,
      sample_basename = sample_basename,
      docker = bcftools_docker
  }
  }

  # Merge per-interval GVCFs
  call MergeVCFs {
    input:
      input_vcfs = CompressAndIndexVCF2.output_vcfgz,
      input_vcfs_indexes = CompressAndIndexVCF2.output_vcfgz_index,
      sample_basename = sample_basename,
      output_filename_gz = output_filename_gz,
      docker = gatk_docker,
      gatk_path = gatk_path
  }

  # Generate variant table
  call GenerateVariantTable {
    input:
      input_vcf = MergeVCFs.output_vcfgz,
      input_vcf_index = MergeVCFs.output_vcfgz_index,
      sample_basename = sample_basename,
      docker = SnpEff_docker
  }

  # Outputs that will be retained when execution is complete
  output {
    File output_vcf = MergeVCFs.output_vcfgz
    File output_vcf_index = MergeVCFs.output_vcfgz_index
    #File output_table = GenerateVariantTable.output_table
  }
}


##################
# TASK DEFINITIONS
##################

task LeftAlignAndTrimVariants {
  input {
    # Command parameters
    File input_vcf
    String sample_basename
    
    File fasta_reference
    File fasta_reference_index
    File fasta_reference_dict
    String gatk_path

    # Runtime parameters
    String docker
  }
  
  command {
  set -e
  
  ~{gatk_path} --java-options -Xmx8G  \
      LeftAlignAndTrimVariants -V ~{input_vcf} -R  ~{fasta_reference} --split-multi-allelics --output ~{sample_basename}.fixed.vcf.gz
  }
  runtime {
    docker: docker
    maxRetries: 3
  }
  output {
    File output_vcf = "~{sample_basename}.fixed.gz"
  }
}

task NormalizeVCF {
  input {
    # Command parameters
    File input_vcf
    String sample_basename
    
    File fasta_reference

    # Runtime parameters
    String docker
  }
  
  command {
  set -e
    bcftools view -Oz ~{input_vcf} > ~{sample_basename}.vcf.gz
    bcftools index -t ~{sample_basename}.vcf.gz
    bcftools norm -m-any -f ~{fasta_reference} ~{sample_basename}.vcf.gz > ~{sample_basename}.normalized.vcf.gz
  }
  runtime {
    docker: docker
    maxRetries: 3
    #requested_memory_mb_per_core: 5000
    #cpu: 1
    #runtime_minutes: 90
    requested_memory_mb_per_core: 2000
    cpu: 4
    runtime_minutes: 120
  }
  output {
    File output_vcf = "~{sample_basename}.normalized.vcf.gz"
  }
}

# Compress VCF if not already and index
task CompressAndIndexVCF {
  input {
    # Command parameters
    File input_vcf
    String sample_basename

    # Runtime parameters
    String docker
  }

  command {
  set -e
    bcftools view -Oz ~{input_vcf} > ~{sample_basename}.vcf.gz
    bcftools index -t ~{sample_basename}.vcf.gz
  }
  runtime {
    docker: docker
    maxRetries: 3
    #requested_memory_mb_per_core: 5000
    #cpu: 1
    #runtime_minutes: 90
    requested_memory_mb_per_core: 1000
    cpu: 8
    runtime_minutes: 120
  }
  output {
    File output_vcfgz = "~{sample_basename}.vcf.gz"
    File output_vcfgz_index = "~{sample_basename}.vcf.gz.tbi"
  }
}

# BCFtools annotate and return gzipped VCF file
task bcftoolsAnnotate {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    String sample_basename

    File HPO
    File HPO_index
    File OMIM
    File OMIM_index
    File gnomadConstraints
    File gnomadConstraints_index
    File CGD
    File CGD_index

    File bcftools_annotation_header

    String? chromosome

    String output_filename

    # Runtime parameters
    String docker
  }

  command {
    set -e
    bcftools view ~{ if defined(chromosome) then " -r " + chromosome else " "} ~{input_vcf} | \
    bcftools annotate -a ~{HPO} -h ~{bcftools_annotation_header} -c CHROM,POS,TO,-,HPO | \
    bcftools annotate -a ~{OMIM}  -h ~{bcftools_annotation_header} -c CHROM,POS,TO,-,OMIM --merge-logic OMIM:unique | \
    bcftools annotate -a ~{gnomadConstraints} -h ~{bcftools_annotation_header} -c CHROM,POS,TO,-,oe_mis,pLI,pRec | \
    bcftools annotate -a ~{CGD} -h ~{bcftools_annotation_header} -c CHROM,POS,TO,-,-,Disease_name,Inheritance,Age,Categorization -Oz -o ~{output_filename}
    
    bcftools index -t ~{output_filename}
  }
  runtime {
    docker: docker
    maxRetries: 3
    requested_memory_mb_per_core: 1000
    cpu: 12
    runtime_minutes: 120
  }
  output {
    File output_vcf = "~{output_filename}"
    File output_vcf_index = "~{output_filename}.tbi"
  }
}

task VCFANNO {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    String sample_basename

    File gnomADexomes_vcf
    File gnomADexomes_vcf_index

    File gnomAD_vcf
    File gnomAD_vcf_index

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

    File SLOpopulation_vcf
    File SLOpopulation_vcf_index

    File ClinVar_vcf
    File ClinVar_vcf_index

    File SpliceAI
    File SpliceAI_index

    File dbscSNV
    File dbscSNV_index

    File pext_bed
    File pext_bed_index

    File gnomAD_pLi_bed
    File gnomAD_pLi_bed_index

    File gnomAD_misz_bed
    File gnomAD_misz_bed_index

    File metaDome_bed
    File metaDome_bed_index

    # Runtime parameters
    String docker
  }
  
  command {
  set -e

  wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/references/missenseConstrained_oe.bed.gz
  wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/references/missenseConstrained_oe.bed.gz.tbi
  wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/references/missenseConstrained_chisq.bed.gz
  wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/references/missenseConstrained_chisq.bed.gz.tbi

  echo [[annotation]] > conf.toml
  echo file=\"~{gnomAD_vcf}\" >> conf.toml
  echo fields = [\"AC\",\"AF\",\"nhomalt\",\"AC_male\",\"nhomalt_male\",\"AC_female\",\"nhomalt_female\",\"AC_nfe_seu\",\"AC_raw\",\"AF_raw\"] >> conf.toml
  echo ops=[\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\"] >> conf.toml
  echo names=[\"gnomAD.AC\",\"gnomAD.AF\",\"gnomAD.nhomalt\",\"gnomAD.AC_male\",\"gnomAD.nhomalt_male\",\"gnomAD.AC_female\",\"gnomAD.nhomalt_female\",\"gnomAD.AC_nfe_seu\",\"gnomAD.AC_raw\",\"gnomAD.AF_raw\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{gnomADexomes_vcf}\" >> conf.toml
  echo fields = [\"AC\",\"AF\",\"nhomalt\",\"AC_male\",\"nhomalt_male\",\"AC_female\",\"nhomalt_female\",\"AC_nfe_seu\",\"AC_raw\",\"AF_raw\"] >> conf.toml
  echo ops=[\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\",\"self\"] >> conf.toml
  echo names=[\"gnomADexomes.AC\",\"gnomADexomes.AF\",\"gnomADexomes.nhomalt\",\"gnomADexomes.AC_male\",\"gnomADexomes.nhomalt_male\",\"gnomADexomes.AC_female\",\"gnomADexomes.nhomalt_female\",\"gnomADexomes.AC_nfe_seu\",\"gnomADexomes.AC_raw\",\"gnomADexomes.AF_raw\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{TopMed_vcf}\" >> conf.toml
  echo fields = [\"AC\",\"AF\",\"AN\",\"Het\",\"Hom\"] >> conf.toml
  echo ops=[\"self\",\"self\",\"self\",\"self\",\"self\"] >> conf.toml
  echo names=[\"TopMed_AC\",\"TopMed_AF\",\"TopMed_AN\",\"TopMed_Het\",\"TopMed_Hom\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{Regeneron_vcf}\" >> conf.toml
  echo fields = [\"ALL_AC\",\"ALL_AF\",\"ALL_AN\"] >> conf.toml
  echo ops=[\"self\",\"self\",\"self\"] >> conf.toml
  echo names=[\"Regeneron_AC\",\"Regeneron_AF\",\"Regeneron_AN\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{GnomAD4_exomes_vcf}\" >> conf.toml
  echo fields = [\"AC_joint\",\"AF_joint\",\"AN_joint\", \"nhomalt_joint\"] >> conf.toml
  echo ops=[\"self\",\"self\",\"self\",\"self\"] >> conf.toml
  echo names=[\"GnomAD4_exomes_AC\",\"GnomAD4_exomes_AF\",\"GnomAD4_exomes_AN\",\"GnomAD4_exomes_Hom\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{GnomAD4_genomes_vcf}\" >> conf.toml
  echo fields = [\"AC_joint\",\"AF_joint\",\"AN_joint\", \"nhomalt_joint\"] >> conf.toml
  echo ops=[\"self\",\"self\",\"self\",\"self\"] >> conf.toml
  echo names=[\"GnomAD4_genomes_AC\",\"GnomAD4_genomes_AF\",\"GnomAD4_genomes_AN\",\"GnomAD4_genomes_Hom\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{SLOpopulation_vcf}\" >> conf.toml
  echo fields = [\"AC\",\"AC_Het\",\"AC_Hom\",\"AC_Hemi\",\"AC_Het\"] >> conf.toml
  echo ops=[\"self\",\"self\",\"self\",\"self\",\"self\"] >> conf.toml
  echo names=[\"SLOpopulation.AC\",\"SLOpopulation.AC_Het\",\"SLOpopulation.AC_Hom\",\"SLOpopulation.AC_Hemi\",\"SLOpopulation.AC_Het\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{ClinVar_vcf}\" >> conf.toml
  echo fields = [\"ID\",\"CLNDN\",\"CLNSIG\",\"CLNHGVS\",\"CLNSIGCONF\",\"CLNSIGINCL\"] >> conf.toml
  echo ops=[\"self\",\"self\",\"self\",\"self\",\"self\",\"self\"] >> conf.toml
  echo names=[\"clinvar.ID\",\"clinvar.CLNDN\",\"clinvar.CLNSIG\",\"clinvar.CLNHGVS\",\"clinvar.CLNSIGCONF\",\"clinvar.CLNSIGINCL\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{SpliceAI}\" >> conf.toml
  echo fields = [\"SpliceAI\"] >> conf.toml
  echo ops=[\"self\"] >> conf.toml
  echo names=[\"SpliceAI.SpliceAI\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{dbscSNV}\" >> conf.toml
  echo columns=[17,18] >> conf.toml
  echo ops=[\"self\",\"self\"] >> conf.toml
  echo names=[\"dbscSNV.ada_score\", \"dbscSNV.rf_score\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{pext_bed}\" >> conf.toml
  echo columns=[4] >> conf.toml
  echo ops=[\"self\"] >> conf.toml
  echo names=[\"pext_score\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"missenseConstrained_oe.bed.gz\" >> conf.toml
  echo columns=[4] >> conf.toml
  echo ops=[\"self\"] >> conf.toml
  echo names=[\"RMC_OE\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{gnomAD_pLi_bed}\" >> conf.toml
  echo columns=[4] >> conf.toml
  echo ops=[\"self\"] >> conf.toml
  echo names=[\"gnomAD_pLi\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{gnomAD_misz_bed}\" >> conf.toml
  echo columns=[4] >> conf.toml
  echo ops=[\"self\"] >> conf.toml
  echo names=[\"gnomAD_misz\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"missenseConstrained_chisq.bed.gz\" >> conf.toml
  echo columns=[4] >> conf.toml
  echo ops=[\"self\"] >> conf.toml
  echo names=[\"RMC_CHISQ\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{blacklisted_regions_bed}\" >> conf.toml
  echo columns=[4] >> conf.toml
  echo ops=[\"lua:replaceSpacesWithUnderscores\(vals\)\"] >> conf.toml
  echo names=[\"blacklisted_regions\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{GERP_bed}\" >> conf.toml
  echo columns=[5] >> conf.toml
  echo ops=[\"lua:calculateMean\(vals\)\"] >> conf.toml
  echo names=[\"GERP_RS\"] >> conf.toml

  echo [[annotation]] >> conf.toml
  echo file=\"~{metaDome_bed}\" >> conf.toml
  echo columns=[4] >> conf.toml
  echo ops=[\"self\"] >> conf.toml
  echo names=[\"metaDome\"] >> conf.toml

  wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/common/annotation/custom.lua
  vcfanno -lua custom.lua -p 16 conf.toml ~{input_vcf} | gzip > ~{sample_basename}.vcf.gz
  }
  
  runtime {
    docker: docker
    maxRetries: 3
    #requested_memory_mb_per_core: 2000
    #cpu: 4
    #runtime_minutes: 240
    requested_memory_mb_per_core: 1000
    cpu: 8
    runtime_minutes: 300
  }
  output {
    File output_vcfgz = "~{sample_basename}.vcf.gz"
  }
}

# SNPeff task
task runSnpEff {
  input {
    # Command parameters
    File input_vcf
    String sample_basename

    File dbNSFP
    File dbNSFP_index

    # Runtime parameters
    String docker
  }

  command { 
  wget -t 1 -T 20 https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/transcripts/canon_transcripts_hg19_refseq.txt
  unset https_proxy
  wget  -t 1 -T 20 https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/transcripts/canon_transcripts_hg19_refseq.txt
  
  java -jar /home/biodocker/bin/snpEff/snpEff.jar -lof -noInteraction -spliceRegionIntronMax 20 -nodownload hg19 ~{input_vcf} > ~{sample_basename}.temp1.vcf 
  java -jar /home/biodocker/bin/snpEff/SnpSift.jar dbNSFP -db ~{dbNSFP} -f REVEL_score,REVEL_rankscore,SIFT_pred,SIFT4G_pred,Polyphen2_HDIV_pred,MutationTaster_pred,MetaSVM_pred,M-CAP_pred,PrimateAI_pred,Aloft_Fraction_transcripts_affected,Aloft_prob_Tolerant,Aloft_prob_Recessive,Aloft_prob_Dominant,Aloft_pred,Aloft_Confidence,CADD_phred,DANN_rankscore,GERP++_NR,GERP++_RS,Interpro_domain,GTEx_V7_gene,GTEx_V7_tissue,Geuvadis_eQTL_target_gene,Polyphen2_HDIV_score -v ~{sample_basename}.temp1.vcf | sed "s/dbNSFP_GERP++/dbNSFP_GERP/g" > ~{sample_basename}.snpEff.vcf

  }
  runtime {
    docker: docker
    maxRetries: 3
    requested_memory_mb_per_core: 1000
    cpu: 16
    runtime_minutes: 120
  }
  output {
    File output_vcf = "~{sample_basename}.snpEff.vcf"
  }
}

# Merge GVCFs generated per-interval for the same sample
task MergeVCFs {
  input {
    # Command parameters
    Array[File] input_vcfs
    Array[File] input_vcfs_indexes
    String sample_basename
    String output_filename_gz

    String gatk_path

    # Runtime parameters
    String docker
  }
  
  command {
  set -e

    ~{gatk_path} --java-options -Xmx4G  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{output_filename_gz}
  }
  runtime {
    docker: docker
    maxRetries: 3
    #requested_memory_mb_per_core: 5000
    #cpu: 1
    #runtime_minutes: 60
    requested_memory_mb_per_core: 1000
    cpu: 8
    runtime_minutes: 120
  }
  output {
    File output_vcfgz = "~{output_filename_gz}"
    File output_vcfgz_index = "~{output_filename_gz}.tbi"
  }
}


# Generate table of variants for interpretation
task GenerateVariantTable {
  input {
    # Command parameters
    File input_vcf
    File input_vcf_index
    String sample_basename

    # Runtime parameters
    String docker
  }

  command {
    set -e
    java -jar /home/biodocker/bin/snpEff/SnpSift.jar extractFields -s "," -e "." ~{input_vcf} CHROM POS REF ALT GEN[0].GT GEN[0].AD GEN[0].DP GEN[0].GQ "ANN[*].GENE" "Disease_name" "Categorization" Inheritance Age HPO "OMIM" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].RANK" "ANN[*].HGVS_P" "ANN[*].IMPACT" "ANN[*].EFFECT"  SLOpopulation.AC_Het SLOpopulation.AC_Hom SLOpopulation.AC_Hemi gnomAD.AC gnomAD.AF gnomAD.nhomalt gnomAD.AC_male gnomAD.nhomalt_male gnomAD.AC_female gnomAD.nhomalt_female gnomAD.AC_nfe_seu gnomAD.AC_raw gnomAD.AF_raw gnomADexomes.AC gnomADexomes.AF gnomADexomes.nhomalt gnomADexomes.AC_male gnomADexomes.nhomalt_male gnomADexomes.AC_female gnomADexomes.nhomalt_female gnomADexomes.AC_nfe_seu gnomADexomes.AC_raw gnomADexomes.AF_raw clinvar.CLNSIG clinvar.CLNDN clinvar.CLNHGVS clinvar.CLNSIGINCL dbNSFP_REVEL_rankscore  dbNSFP_MetaSVM_pred dbNSFP_CADD_phred  dbNSFP_DANN_rankscore dbNSFP_SIFT_pred  dbNSFP_SIFT4G_pred  dbNSFP_Polyphen2_HDIV_pred  dbNSFP_MutationTaster_pred dbNSFP_PrimateAI_pred  dbNSFP_Polyphen2_HDIV_score SpliceAI.SpliceAI dbscSNV.ada_score dbscSNV.rf_score dbNSFP_GERP___NR dbNSFP_GERP___RS dbNSFP_GTEx_V7_gene dbNSFP_GTEx_V7_tissue dbNSFP_Interpro_domain dbNSFP_Geuvadis_eQTL_target_gene dbNSFP_Aloft_Fraction_transcripts_affected  dbNSFP_Aloft_prob_Tolerant  dbNSFP_Aloft_prob_Recessive  dbNSFP_Aloft_prob_Dominant  dbNSFP_Aloft_pred  dbNSFP_Aloft_Confidence pLI oe_mis pRec "LOF[*].GENE" "LOF[*].GENEID" "LOF[*].NUMTR" "LOF[*].PERC" "NMD[*].GENE" "NMD[*].GENEID" "NMD[*].NUMTR" "NMD[*].PERC"  > ~{sample_basename}.tab
  }
  runtime {
    docker: docker
    maxRetries: 3
    #requested_memory_mb_per_core: 4000
    #cpu: 2
    #runtime_minutes: 60
    requested_memory_mb_per_core: 1000
    cpu: 8
    runtime_minutes: 120
  }
  output {
    File output_table = " ~{sample_basename}.tab"
  }
}

task StringToArray {
  input {
    String input_string
    String separator
  }
  command <<<
    echo '~{input_string}' | tr '~{separator}' \\n | tr -d "[:blank:]" > intervals.list
    echo '~{input_string}' | tr '~{separator}' \\n | tr -d "[:blank:]"
  >>>
  runtime {
    docker:"biocontainers/bcftools:v1.9-1-deb_cv1"
    requested_memory_mb_per_core: 500
    cpu: 1
    runtime_minutes: 5
  }
  output {
    Array[String] values = read_lines(stdout())
    File intervals_list = "intervals.list"
  }
}

