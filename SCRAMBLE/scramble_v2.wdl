version 1.0
## Copyright CMG@KIGM, Ales Maver

### import "../FastqToVCFPipeline_3.wdl" as FastqToVcf
import "../CRAM_conversions.wdl" as CramToBam
import "../manta/manta_workflow.wdl" as Manta

# WORKFLOW DEFINITION
workflow SCRAMBLE_workflow {
    input {
        # input files are either bam or cram
        File? input_bam
        File? input_bam_index
        File? input_cram
        File? input_cram_index
        # Fields required for converting CRAM TO BAM  -- mandatory because of the definition in original CramToBam task
        File reference_fa
        File reference_fai
        File reference_dict
    }
    
    String gitc_docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"
    String samtools_path = "samtools" # Path to samtools command within GITC docker
    
    String sample_basename = sub(basename(select_first([input_bam, input_cram])), "[\_,\.].*", "" )
    
    if ( defined(input_cram) ) {
        call CramToBam.CramToBam as Cram_hg19_ToBam {
            input:
                input_cram = input_cram,
                sample_name = sample_basename,
                ref_dict = reference_dict,
                ref_fasta = reference_fa,
                ref_fasta_index = reference_fai,
                docker = gitc_docker,
                samtools_path = samtools_path
        }
    }

    call SCRAMBLE {
        input:
            #input_bam=Cram_hg19_ToBam.output_bam,
            #input_bam_index=Cram_hg19_ToBam.output_bai,
            sample_basename = sample_basename,
            input_bam = select_first([input_bam, Cram_hg19_ToBam.output_bam]),
            input_bam_index = select_first([input_bam_index, Cram_hg19_ToBam.output_bai])
    }

    call Manta.annotSV as MEI_annotSV {
      input:
        genome_build = "GRCh37",
        input_vcf = SCRAMBLE.output_meis_bed,
        output_tsv_name = sample_basename + ".MEIs.annotSV.tsv"
    }

    output {
        File MEI_output = SCRAMBLE.output_meis
        File? MEI_output_annotated = MEI_annotSV.sv_variants_tsv
    }

}

task SCRAMBLE {
  input {
    # Command parameters
    File input_bam
    File input_bam_index
    String sample_basename

    # Runtime parameters
    String docker = "alesmaver/scramble2"
  }

  command <<<
    cluster_identifier ~{input_bam} > clusters.txt

    Rscript --vanilla /app/cluster_analysis/bin/SCRAMble.R \
        --out-name $PWD/output \
        --cluster-file "$(realpath clusters.txt)" \
        --install-dir /app/cluster_analysis/bin \
        --mei-refs /app/cluster_analysis/resources/MEI_consensus_seqs.fa \
        --eval-meis \
        --no-vcf

    cp $PWD/output_MEIs.txt ~{sample_basename}.output_MEIs.txt
    awk -F'[:\t]' -v OFS="\t" '{print $1,$2,$2+1,"INS",$0}' ~{sample_basename}.output_MEIs.txt > ~{sample_basename}.MEIs.bed
  >>>
  
  runtime {
    docker: docker
    requested_memory_mb_per_core: 2000
    cpu: 6
    runtime_minutes: 300
    continueOnReturnCode: true
  }
  output {
    File output_meis = "~{sample_basename}.output_MEIs.txt"
    File output_meis_bed = "~{sample_basename}.MEIs.bed"
  }
}
