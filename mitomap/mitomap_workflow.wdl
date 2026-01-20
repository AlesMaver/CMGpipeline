version 1.0
## Copyright CMG@KIGM, Ales Maver

## Usage
# This workflow accepts three types of inputs: Illumina FASTQ files, a BAM file, or CRAM files
# Currently, the input CRAM files should be aligned to the hg19 reference genome assembly, we will implement support for other genome formats in the future
# The CRAM output is optional and disabled by default at the moment, until production switches to CRAM

## Subworkflows
## import "https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/MitoMap.wdl" as MitoMap

# WORKFLOW DEFINITION 
workflow MitoMapWorkflow {
  input {
    File input_vcf
    String sample_basename

    File reference_fa
    File reference_fai
    File reference_dict
  }  
  
  call AnalyseInputVcf {
    input:
      input_vcf = input_vcf  
  }  
  
  call CreateMitoFasta {
    input:
      input_vcf = input_vcf,
      sample_basename = sample_basename,

      reference_fa = reference_fa,
      reference_fai = reference_fai,
      reference_dict = reference_dict,

      docker = "broadinstitute/gatk3:3.8-1",
      
      variant_exists = AnalyseInputVcf.variant_exists
  }

  call MitoMap {
    input:
      mtDNA_fasta = CreateMitoFasta.mtDNA_fasta,
      sample_basename = sample_basename,
      variant_exists = AnalyseInputVcf.variant_exists
  }

  output {
    File? mitoResults_xls = MitoMap.mitoResults_xls
    File? mitoResults_txt = MitoMap.mitoResults_txt

  }
}


# Let's check if input VCF file has any mitochondrial variants.
task AnalyseInputVcf {
    input {
        File input_vcf
    }   
 
    command <<<
        echo Analysing ~{input_vcf}      
        grep '^chrM' ~{input_vcf} | wc -l > variant_count.txt
        echo Number of mitochondrial variants:
        cat variant_count.txt
        echo ------
        vcf_file=~{input_vcf}
        mito_variants=$(grep '^chrM' $vcf_file)
        variant_count=$(echo "$mito_variants" | wc -l)
        echo Number of mitochondrial variants: $variant_count
        if [ $variant_count -gt 1 ]
        then
            echo True > variant_exists.txt
        else
            echo False > variant_exists.txt
        fi
        cat variant_exists.txt
    >>>
    
    output {
        Int variant_count = read_int("variant_count.txt")
        Boolean variant_exists = read_boolean("variant_exists.txt")
    }

    runtime {
        docker: "broadinstitute/gatk3:3.8-1"
        cpu: 2
        requested_memory_mb_per_core: 1000
        runtime_minutes: 5
    }            
}



## An additional option for calculating coverage using GATK
task CreateMitoFasta {
    input {
        File input_vcf
        String sample_basename

        File reference_fa
        File reference_fai
        File reference_dict

        String docker
        
        Boolean variant_exists
    }

    command {
      set -e
      echo ~{variant_exists}
      do_GenomeAnalysis=~{variant_exists}
      echo Perform GenomeAnalysis: $do_GenomeAnalysis
      if [ "$do_GenomeAnalysis" = true ]
      then
        java -Xmx3g -jar /usr/GenomeAnalysisTK.jar \
          -T FastaAlternateReferenceMaker \
          -R ~{reference_fa} \
          --variant ~{input_vcf} \
          -o mtDNA.fasta \
          -L chrM
      else
        touch mtDNA.fasta
      fi
      ls -ls mtDNA.fasta
    }
    
    output {
    File mtDNA_fasta = "mtDNA.fasta"
    }

    runtime {
        docker: "~{docker}"
        maxRetries: 3
        requested_memory_mb_per_core: 1000
        cpu: 8
        runtime_minutes: 60
    }
}

task MitoMap {
    input {
        File mtDNA_fasta
        String sample_basename
        Boolean variant_exists
    }

    command <<<
      echo ~{variant_exists}
      do_mitomap=~{variant_exists}
      echo Perform MitoMap: $do_mitomap
      if [ "$do_mitomap" = true ]
      then
        #cp /usr/src/app/mitomap.py ./
        wget --no-check-certificate https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/mitomap/mitomap.py
        cp ~{mtDNA_fasta} ./
        python mitomap.py > ~{sample_basename}_mitoResults.txt
      else
        touch ~{sample_basename}_mitoResults.txt
        echo "query	tpos	qpos	tnt	qnt	ntchange	allele	calc_locus	calc_aachange	conservation	haplogroup	verbose_haplogroup	patientphenotype	mmutid	rtmutid	polyid	subvar_cnt	is_polymorphism	is_mmut	is_rtmut	is_submitted	gb_cnt	gb_perc	hap_cnt	hap_perc" > ~{sample_basename}_mitoResults.txt 
      fi
      cp ~{sample_basename}_mitoResults.txt ~{sample_basename}_mitoResults.xls
      ls -ls ~{sample_basename}_mitoResults.*
    >>>

    output {
      File mitoResults_txt = "~{sample_basename}_mitoResults.txt"
      File mitoResults_xls = "~{sample_basename}_mitoResults.xls"
    }

    runtime {
        docker: "alesmaver/mitomap"
        maxRetries: 3
        requested_memory_mb_per_core: 1000
        cpu: 3
        runtime_minutes: 60
    }
}

