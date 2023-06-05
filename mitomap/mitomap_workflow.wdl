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

      docker = "broadinstitute/gatk3:3.8-1"
  }
  
  call MitoMap {
    input:
    mtDNA_fasta = CreateMitoFasta.mtDNA_fasta,
    sample_basename = sample_basename
  }

  output {
    File? mitoResults_xls = MitoMap.mitoResults_xls
    File? mitoResults_txt = MitoMap.mitoResults_txt

  }
}


# Let's check if input VCF file has any mitochondrial variants.
# number greater than 0 means we have variants

task AnalyseInputVcf {
    input {
        File input_vcf
    }   
 
    command {
        echo Analysing ~{input_vcf}      
        grep '^chrM' ~{input_vcf} | wc -l > variant_count.txt
        echo Number of mitochondrial variants:
        cat variant_count.txt
        echo ---
        filename=~{input_vcf}
        variant_count=$(grep '^chrM' $filename | wc -l)
        echo $variant_count
        if [ $variant_count -gt 0 ]
        then
          echo True > variant_exists.txt
        else
          echo False > variant_exists.txt
        fi
        cat variant_exists.txt
        
        echo "query   tpos    qpos    tnt     qnt     ntchange        allele  calc_locus      calc_aachange   conservation    haplogroup      verbose_haplogroup    patientphenotype        mmutid  rtmutid polyid  subvar_cnt      is_polymorphism is_mmut is_rtmut        is_submitted gb_cnt   gb_perc hap_cnt hap_perc" > header.txt
    }
    
    output {
        Int variant_count = read_int("variant_count.txt")
        Boolean variant_exists = read_boolean("variant_exists.txt")
    }

    runtime {
        docker: "broadinstitute/gatk3:3.8-1"
        cpu: 1
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
    }

    command {
    set -e
    java -Xmx3g -jar /usr/GenomeAnalysisTK.jar \
      -T FastaAlternateReferenceMaker \
      -R ~{reference_fa} \
      --variant ~{input_vcf} \
      -o mtDNA.fasta \
      -L chrM
     }
    
    output {
    File mtDNA_fasta = "mtDNA.fasta"
    }

    runtime {
        docker: "~{docker}"
        maxRetries: 3
        requested_memory_mb_per_core: 9000
        cpu: 1
        runtime_minutes: 60
    }
}

task MitoMap {
    input {
        File mtDNA_fasta
        String sample_basename
    }

    command <<<
    
      #cp /usr/src/app/mitomap.py ./
      wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/mitomap/mitomap.py
      cp ~{mtDNA_fasta} ./
      python mitomap.py > ~{sample_basename}_mitoResults.txt
      cp ~{sample_basename}_mitoResults.txt ~{sample_basename}_mitoResults.xls
    >>>

    output {
      File mitoResults_txt = "~{sample_basename}_mitoResults.txt"
      File mitoResults_xls = "~{sample_basename}_mitoResults.xls"
    }

    runtime {
        docker: "alesmaver/mitomap"
        maxRetries: 3
        requested_memory_mb_per_core: 3000
        cpu: 1
        runtime_minutes: 60
    }
}

