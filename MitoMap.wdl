## # QualiMap
##
## This WDL tool wraps the [QualiMap](http://qualimap.bioinfo.cipf.es/) tool.
## QualiMap computes metrics to facilitate evaluation of sequencing data. 
version 1.0

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
        # requested_memory_mb_per_core: 9000
        # cpu: 1
        requested_memory_mb_per_core: 1000
        cpu: 10
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
        requested_memory_mb_per_core: 1000
        cpu: 3
        runtime_minutes: 60
    }
}
