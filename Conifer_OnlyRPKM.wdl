version 1.0
## Copyright CMG@KIGM, Ales Maver
import "./Conifer.wdl" as Conifer
import "./CRAM_conversions.wdl" as CramConversions

# WORKFLOW DEFINITION 
workflow Conifer_OnlyRPKM {
  input {
    String? sample_basename
    File? input_bam
    File? input_bam_index
    File? input_cram
    File? input_cram_index
    File? reference_fa
    File? reference_fai
    File? reference_dict

    String? enrichment
    File? enrichment_bed
  }  

  String sample_name = select_first([sample_basename, sub(basename(select_first([input_bam, input_cram, [""]])), "[\_,\.].*", "")  ])

  if (defined(input_cram)) {
    call CramConversions.CramToBam as CramToBam {
        input:
          sample_name = sample_name,
          input_cram = input_cram,
          ref_fasta = reference_fa,
          ref_fasta_index = reference_fai,
          ref_dict = reference_dict,
          docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817",
          samtools_path = "samtools"
    }
  }
  
  call Conifer.MakeRPKM as MakeRPKM {
      input:
        input_bam=select_first([input_bam, CramToBam.output_bam]),
        input_bam_index=select_first([input_bam_index, CramToBam.output_bai]),
        sample_basename=sample_name,
        enrichment=enrichment,
        enrichment_bed=enrichment_bed
  }

  output {
    File output_rpkm = MakeRPKM.output_rpkm
  }
}
