## # QualiMap
##
## This WDL tool wraps the [QualiMap](http://qualimap.bioinfo.cipf.es/) tool.
## QualiMap computes metrics to facilitate evaluation of sequencing data. 
version 1.0

task bamqc {
    input {
        File bam
        File? enrichment_bed
        String sample_basename
        Int ncpu = 5
        Int max_retries = 1
        Int memory_gb = 100            # 25
        Int? disk_size_gb
    }

    String out_directory = basename(bam, ".bam") + '_qualimap_bamqc_results'
    String out_tar_gz_file = sample_basename + ".Qualimap.tar.gz"

    Int java_heap_size = ceil(memory_gb * 0.9)
    Float bam_size = size(bam, "GiB")
    Int disk_size = select_first([disk_size_gb, ceil((bam_size * 2) + 10)])

    command {
        set -euo pipefail
        
        unset DISPLAY
        
        qualimap bamqc -bam ~{bam} \
            ~{if defined(enrichment_bed) then "-gff " + enrichment_bed else ""} \
            -outdir ~{out_directory} \
            -nt ~{ncpu} \
            -nw 400 \
            --java-mem-size=~{java_heap_size}g
        
        # Check if qualimap succeeded
        if [ ! -d "~{out_directory}/raw_data_qualimapReport/" ]; then
            exit 1
        fi

        tar -czf ~{out_tar_gz_file} ~{out_directory}
    }

    runtime {
        memory: memory_gb + " GB"
        disk: disk_size + " GB"
        docker: 'stjudecloud/qualimap:1.0.0'
        maxRetries: max_retries
        requested_memory_mb_per_core: 1000   ## 5000
        cpu: 64
        ## runtime_minutes: 180
        runtime_minutes: 900
    }

    output {
        File results = out_tar_gz_file
    }

    meta {
        author: "Andrew Thrasher, Andrew Frantz"
        email: "andrew.thrasher@stjude.org, andrew.frantz@stjude.org"
        description: "This WDL tool runs QualiMap's bamqc tool on the input BAM file."
    }

    parameter_meta {
        bam: "Input BAM format file to generate coverage for"
    }
}

## GATK coverage calculation can output per-gene coverage based on the refSeq models in the UCSC format 
## For this purpose a refSeqFile is needed and is created in the following way
## 1. Download refGene models from using the table option in the UCSC browser, use the standard UCSC output format 
## 2. Then run the following script
## rm -rf tmp
## mkdir -p tmp
## # Divide all chromosomes into separate files.
## gawk '($3 ~ /_/ || $0 ~ /^#/){next;}{print > "./tmp/tmp.sort."$3""}' refGene.select.refseq
## 
## # initiate an empty output
## echo -n "" > refGene.select.sorted.refseq
## head -n1 refGene.select.refseq > refGene.select.sorted.refseq
## 
## for chr in chr{1..22} chr{X,Y} ; do
##    # sort each chromsome file by transcript starts (field number 5)
##     sort  -k5,5g -s  ./tmp/tmp.sort.$chr >> refGene.select.sorted.refseq
## done 

## An additional option for calculating coverage using GATK
task DepthOfCoverage {
    input {
        File input_bam
        File input_bam_index
        String sample_basename

        File reference_fa
        File reference_fai
        File reference_dict
        
        File? enrichment_bed
        File refSeqFile
                
        # Runtime params
        String gatk_path
        Int threads
        String docker
    }

    command <<<
    set -e
    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      DepthOfCoverage \
      -R ~{reference_fa} \
      -I ~{input_bam} \
      -O targetGenes.coverage \
      ~{if defined(enrichment_bed) then "-L " + enrichment_bed else ""} \
       --omit-depth-output-at-each-base \
       -ip 2 \
       --summary-coverage-threshold 5 \
       --summary-coverage-threshold 10 \
       --summary-coverage-threshold 15 \
       --summary-coverage-threshold 20 \
       --summary-coverage-threshold 30 \
       --summary-coverage-threshold 40 \
       --summary-coverage-threshold 50 \
       --summary-coverage-threshold 60 \
       --summary-coverage-threshold 70 \
       --summary-coverage-threshold 80 \
       --summary-coverage-threshold 90 \
       --summary-coverage-threshold 100 \
       --gene-list ~{refSeqFile}

    cat targetGenes.coverage.sample_interval_summary | grep -v "Target" | awk -F '[\t:-]' '{print $1,$2,$3,$5}' OFS='\t' > ~{sample_basename}.coverage_mean.wig
    cat targetGenes.coverage.sample_interval_summary | grep -v "Target" | awk -F '[\t:-]' '{print $1,$2,$3,$11}' OFS='\t' > ~{sample_basename}.coverage.wig
    cat targetGenes.coverage.sample_interval_summary | grep -v "Target" | awk -F '[\t:-]' '{print $1,$2,$3,($11-100)}' OFS='\t' > ~{sample_basename}.coverage_neg.wig

    ~{gatk_path} --java-options "-Xmx8g -XX:ParallelGCThreads=~{threads}"  \
      DepthOfCoverage \
       -R ~{reference_fa} \
       -I ~{input_bam} \
       -O mitochondrial.coverage \
       -L chrM \
       --omit-depth-output-at-each-base \
       -ip 2 \
       --summary-coverage-threshold 5 \
       --summary-coverage-threshold 10 \
       --summary-coverage-threshold 15 \
       --summary-coverage-threshold 20 \
       --summary-coverage-threshold 30 \
       --summary-coverage-threshold 40 \
       --summary-coverage-threshold 50 \
       --summary-coverage-threshold 60 \
       --summary-coverage-threshold 70 \
       --summary-coverage-threshold 80 \
       --summary-coverage-threshold 90 \
       --summary-coverage-threshold 100

    cat mitochondrial.coverage.sample_interval_summary | grep -v "Target" | awk -F '[\t:-]' '{print $1,$2,$3,$5}' OFS='\t' > ~{sample_basename}.mitochondrial.coverage_mean.wig
    cat mitochondrial.coverage.sample_interval_summary | grep -v "Target" | awk -F '[\t:-]' '{print $1,$2,$3,$11}' OFS='\t' > ~{sample_basename}.mitochondrial.coverage.wig
    cat mitochondrial.coverage.sample_interval_summary | grep -v "Target" | awk -F '[\t:-]' '{print $1,$2,$3,($11-100)}' OFS='\t' > ~{sample_basename}.mitochondrial.coverage_neg.wig

    tar -czf ~{sample_basename}.DepthOfCoverage.tar.gz *coverage*
    >>>
    
    output {
    File DepthOfCoverage_output = "~{sample_basename}.DepthOfCoverage.tar.gz"
    }

    runtime {
        docker: "~{docker}"
        maxRetries: 3
        requested_memory_mb_per_core: 1000
        cpu: 12
        #runtime_minutes: 6000
        runtime_minutes: 9000
    }
}

## An additional option for calculating coverage using GATK
task DepthOfCoverage34 {
    input {
        File input_bam
        File input_bam_index
        String sample_basename

        File reference_fa
        File reference_fai
        File reference_dict
        
        File? enrichment_bed
        File refSeqFile
                
        # Runtime params
        String gatk_path
        Int threads
        String docker
    }

    command <<<
    set -e
    java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \
      -T DepthOfCoverage \
      -R ~{reference_fa} \
      -I ~{input_bam} \
      -o targetGenes.coverage \
      ~{if defined(enrichment_bed) then "-L " + enrichment_bed else ""} \
       -omitBaseOutput \
       -ip 2 \
       -allowPotentiallyMisencodedQuals \
       -ct 5 \
       -ct 10 \
       -ct 20 \
       -ct 50 \
       -ct 100 \
       -geneList ~{refSeqFile}

    cat targetGenes.coverage.sample_interval_summary | grep -v "Target" | awk -F '[\t:-]' '{print $1,$2,$3,$5}' OFS='\t' | sed 's/NaN/0.00/g' > ~{sample_basename}.coverage_mean.wig
    cat targetGenes.coverage.sample_interval_summary | grep -v "Target" | awk -F '[\t:-]' '{print $1,$2,$3,$11}' OFS='\t' | sed 's/NaN/0.00/g' > ~{sample_basename}.coverage.wig
    cat targetGenes.coverage.sample_interval_summary | grep -v "Target" | awk -F '[\t:-]' '{print $1,$2,$3,($11-100)}' OFS='\t' | sed 's/NaN/0.00/g' > ~{sample_basename}.coverage_neg.wig

    java -Xmx8g -jar /usr/GenomeAnalysisTK.jar \
      -T DepthOfCoverage \
      -R ~{reference_fa} \
      -I ~{input_bam} \
      -o mitochondrial.coverage \
       -L chrM \
       -omitBaseOutput \
       -ip 2 \
        -allowPotentiallyMisencodedQuals \
       -ct 5 \
       -ct 10 \
       -ct 20 \
       -ct 50 \
       -ct 100

    cat mitochondrial.coverage.sample_interval_summary | grep -v "Target" | awk -F '[\t:-]' '{print $1,$2,$3,$5}' OFS='\t' | sed 's/NaN/0.00/g' > ~{sample_basename}.mitochondrial.coverage_mean.wig
    cat mitochondrial.coverage.sample_interval_summary | grep -v "Target" | awk -F '[\t:-]' '{print $1,$2,$3,$11}' OFS='\t' | sed 's/NaN/0.00/g' > ~{sample_basename}.mitochondrial.coverage.wig
    cat mitochondrial.coverage.sample_interval_summary | grep -v "Target" | awk -F '[\t:-]' '{print $1,$2,$3,($11-100)}' OFS='\t' | sed 's/NaN/0.00/g' > ~{sample_basename}.mitochondrial.coverage_neg.wig

    tar -czf ~{sample_basename}.DepthOfCoverage.tar.gz *coverage*
    >>>
    
    output {
    File DepthOfCoverage_output = "~{sample_basename}.DepthOfCoverage.tar.gz"
    File coverage_mean_wig = "~{sample_basename}.coverage_mean.wig"
    File coverage_wig = "~{sample_basename}.coverage.wig"
    File coverage_neg_wig = "~{sample_basename}.coverage_neg.wig"
    }

    runtime {
        docker: "~{docker}"
        maxRetries: 3
        requested_memory_mb_per_core: 1000
        cpu: 12
        runtime_minutes: 2400
    }
}


task DownsampleBED {
  input {
    File? bed_file
    File? reference_fai

    Int select_every_nth_line = 2
    String docker = "pegi3s/bedtools"
  }
  
  String bed_filename = basename(select_first([bed_file, ""]))
  File reference_fai_filename = select_first([reference_fai, ""])

  command <<<
  set -e
  #awk {'print $1, $2'} OFS='\t' ~{reference_fai_filename} |head -n 26 | grep -v '_' | grep -v 'chrM' | grep -v 'chrY' > hg19.genome
  awk {'print $1, $2'} OFS='\t' ~{reference_fai_filename} |head -n 26 | grep -v '_'  > hg19.genome
  bedtools makewindows -w 100000 -g hg19.genome  > ~{bed_filename}
  #cat ~{bed_file} | awk 'NR % ~{select_every_nth_line} == 0' > downsampled_~{bed_filename}
  cat ~{bed_filename} | awk -F'\t' '{print $1,$2+5,$3-5}' OFS='\t' > downsampled_~{bed_filename}
  >>>

  runtime {
    docker: docker
    requested_memory_mb_per_core: 1000
    cpu: 2
    runtime_minutes: 60
  }
  output {
    File? downsampled_bed_file = "downsampled_~{bed_filename}"
  }
}
