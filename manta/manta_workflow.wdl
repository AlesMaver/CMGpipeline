version 1.0

# MIT License
#
# Copyright (c) 2018 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


import "./manta.wdl" as manta
import "../AnnotationPipeline.wdl" as Annotation
import "../CRAM_conversions.wdl" as CramConversions

workflow SVcalling {
    input {
        File? bamFile
        File? bamIndex
        File? cramFile
        File? cramIndex
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        String sample
        String newId = "\'%CHROM\\_%POS\'"
        Boolean excludeMisHomRef = false
        Boolean excludeFpDupDel = false
        String outputDir = "."
        
        # when calling manta germline: genome (exome=false), exome (exome=true)
        Boolean exome = false

        Array[File] input_manta_reference_vcfs = select_first([input_manta_reference_vcfs, [""]])

        Map[String, String] dockerImages = {
            "bcftools": "quay.io/biocontainers/bcftools:1.10.2--h4f4756c_2",
            "clever": "quay.io/biocontainers/clever-toolkit:2.4--py36hcfe0e84_6",
            "delly": "quay.io/biocontainers/delly:0.8.1--h4037b6b_1",
            "gridss": "quay.io/biocontainers/gridss:2.9.4--0",
            "manta": "quay.io/biocontainers/manta:1.4.0--py27_1",
            "picard":"quay.io/biocontainers/picard:2.23.2--0",
            "samtools": "quay.io/biocontainers/samtools:1.10--h9402c20_2",
            "survivor": "quay.io/biocontainers/survivor:1.0.6--h6bb024c_0",
            "smoove": "quay.io/biocontainers/smoove:0.2.5--0",
            "duphold": "quay.io/biocontainers/duphold:0.2.1--h516909a_1"
        }
    }

    meta {allowNestedInputs: true}

    if (defined(cramFile)) {
      call CramConversions.CramToBam as CramToBam {
        input:
          sample_name = sample,
          input_cram = cramFile,
          ref_fasta = referenceFasta,
          ref_fasta_index = referenceFastaFai,
          ref_dict = referenceFastaDict,
          docker = "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817",
          samtools_path = "samtools"
      }
    }

    String SVdir = outputDir + '/structural-variants/'

    call DownloadFiles

    call manta.Germline as manta {
        input:
            dockerImage = dockerImages["manta"],
            bamFile = select_first([bamFile, CramToBam.output_bam]),
            bamIndex = select_first([bamIndex, CramToBam.output_bai]),
            sample = sample,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            callRegions = DownloadFiles.callRegions,
            callRegionsIndex = DownloadFiles.callRegionsIndex,
            runDir = SVdir + 'manta/',
            exome = exome
    }

    call MergeMantaFiles {
        input:
            input_manta_vcf = manta.mantaVCF,
            input_manta_reference_vcfs = input_manta_reference_vcfs,
            sample_basename = sample
    }    

    call AnnotateMantaVCF {
        input:
            input_vcf = MergeMantaFiles.merged_vcf,
            sample_basename = sample
    }

    call annotSV {
        input:
            genome_build = "GRCh37",
            input_vcf = AnnotateMantaVCF.output_manta_filtered_vcf,
            output_tsv_name = sample + ".manta.AnnotSV.tsv"
    }

    output {
        File mantaVcf = manta.mantaVCF
        File mantaVcfindex = manta.mantaVCFindex
        File output_sv_table = AnnotateMantaVCF.output_sv_table
        File output_manta_filtered_vcf = AnnotateMantaVCF.output_manta_filtered_vcf
        File? annotSV_tsv = annotSV.sv_variants_tsv
    }

    parameter_meta {
        outputDir: {description: "The directory the output should be written to.", category: "common"}
        referenceFasta: { description: "The reference fasta file", category: "required" }
        referenceFastaFai: { description: "Fasta index (.fai) file of the reference", category: "required" }
        referenceFastaDict: { description: "Sequence dictionary (.dict) file of the reference", category: "required" }
        bamFile: {description: "sorted BAM file", category: "required"}
        bamIndex: {description: "BAM index(.bai) file", category: "required"}
        sample: {description: "The name of the sample", category: "required"}
        newId: {description: "Assign ID on the fly (e.g. --set-id +'%CHROM\_%POS').", category: "advanced"}
        excludeMisHomRef: {description: "Option to exclude missing and homozygous reference genotypes.", category: "advanced"}
        excludeFpDupDel: {description: "Option to exclude false positive duplications and deletions according to DUPHOLD.", category: "advanced"}
        dockerImages: {description: "A map describing the docker image used for the tasks.",
                           category: "advanced"}
    }
}

###
###
###
task DownloadFiles {
  command {
    wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/references/blacklisted_regions/callable_regions.hg19.bed.gz
    wget https://raw.githubusercontent.com/AlesMaver/CMGpipeline/master/references/blacklisted_regions/callable_regions.hg19.bed.gz.tbi
  }
  output {
    File callRegions = "callable_regions.hg19.bed.gz"
    File callRegionsIndex = "callable_regions.hg19.bed.gz.tbi"
  }
  runtime {
    docker: "alesmaver/exome_depth"
  }
}

task MergeMantaFiles {
  input {
    # Command parameters
    File input_manta_vcf
    Array[File] input_manta_reference_vcfs 
    String sample_basename
  }

  command <<<
  MANTA_VCFS_DIR=$(dirname ~{input_manta_reference_vcfs[0]})
  cp ~{input_manta_vcf} $MANTA_VCFS_DIR
  cp $MANTA_VCFS_DIR/*.vcf ./
  cp $MANTA_VCFS_DIR/*.vcf.gz ./

  gunzip ./*manta.vcf.gz
  for file in *manta.vcf; do
    awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' $file > "${file/vcf/filtered.vcf}"
  done

  ls -d ./*manta.filtered.vcf > fileList

  SURVIVOR merge fileList 500 1 1 0 0 20 ./merged.vcf
  >>>

  runtime {
    docker: "quay.io/biocontainers/survivor:1.0.6--h6bb024c_0"
    requested_memory_mb_per_core: 1000
    cpu: 6
    runtime_minutes: 120
    docker_user: "root"
  }
  output {
    File merged_vcf = "merged.vcf"
  }
}

# SNPeff task
task AnnotateMantaVCF {
  input {
    # Command parameters
    File input_vcf
    String sample_basename
  }

  command { 
  java -jar /home/biodocker/bin/snpEff/snpEff.jar -noInteraction -noDownload hg19 ~{input_vcf} > ~{sample_basename}.merged.annotated.vcf
  java -jar /home/biodocker/bin/snpEff/SnpSift.jar filter "(ANN[*].IMPACT has 'MODERATE' | ANN[*].IMPACT has 'HIGH') & SUPP < 3 & isVariant(GEN[~{sample_basename}].GT)" ~{sample_basename}.merged.annotated.vcf > ~{sample_basename}.merged.annotated.filtered.vcf
  java -jar /home/biodocker/bin/snpEff/SnpSift.jar extractFields -s "," -e "." ~{sample_basename}.merged.annotated.filtered.vcf CHROM POS REF ALT SVLEN QUAL ANN[*].GENE ANN[1].IMPACT ANN[1].EFFECT GEN[*].GT > ~{sample_basename}.mantaSVs.txt

  # Create an unannotated snpEff file for annotSV and other applications
  java -jar /home/biodocker/bin/snpEff/SnpSift.jar filter "SUPP < 3 & isVariant(GEN[~{sample_basename}].GT)" ~{input_vcf} > ~{sample_basename}.merged.filtered.vcf

  }
  runtime {
    docker: "alesmaver/snpeff_v50"
    maxRetries: 3
    requested_memory_mb_per_core: 2000
    cpu: 8
    runtime_minutes: 120
  }
  output {
    File output_sv_table = "~{sample_basename}.mantaSVs.txt"
    File output_manta_filtered_vcf = "~{sample_basename}.merged.filtered.vcf"
    File output_manta_annotated_filtered_vcf = "~{sample_basename}.merged.annotated.filtered.vcf"
  }
}

task annotSV {
  input {
    String genome_build
    File input_vcf
    String output_tsv_name = "AnnotSV.tsv"
  }

  runtime {
    requested_memory_mb_per_core: 1000
    #docker: "mgibio/annotsv-cwl:2.1" # Disabling the old image and trying the 3.1 image due to inclusion of cytobands in the annotation
    # The image was created by:
    # 1. Pulling trinhanne/annotsv:3.1 image
    # 2. cd-ing into the /opt/AnnotSV_3.1 directory
    # 3. Running the make PREFIX=. install-human-annotation
    docker: "alesmaver/annotsv"
    # continueOnReturnCode: true
    cpu: 8
  }

  command <<<
    touch ~{output_tsv_name}
    /opt/AnnotSV_3.1/bin/AnnotSV -bedtools /usr/bin/bedtools -outputDir "$PWD" \
    -genomeBuild ~{genome_build} \
    -SVinputFile ~{input_vcf} \
    -outputFile ~{output_tsv_name}
    echo "Done."
  >>>

  output {
    File? sv_variants_tsv = output_tsv_name
  }
}
