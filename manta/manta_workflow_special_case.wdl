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


#import "./manta.wdl" as manta
import "./manta_workflow.wdl" as manta_wf
import "../AnnotationPipeline.wdl" as Annotation

# T.M. 15.11.2002
# This is a special case of manta_workflow.wdl
# we want to produce manta annotSV file from the manta.vcf of the previous run
#  


workflow SVcalling_wo_Germline {
    input {
        # File bamFile
        # File bamIndex
        File input_manta_vcf
        
        #File referenceFasta
        #File referenceFastaFai
        #File referenceFastaDict
        
        String sample
        
        #String newId = "\'%CHROM\\_%POS\'"
        #Boolean excludeMisHomRef = false
        #Boolean excludeFpDupDel = false
        
        String outputDir = "."
        
        ## when calling manta germline: genome (exome=false), exome (exome=true)
        #Boolean exome = false

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

    String SVdir = outputDir + '/structural-variants/'

    #call manta.Germline as manta {
    #    input:
    #        dockerImage = dockerImages["manta"],
    #        bamFile = bamFile,
    #        bamIndex = bamIndex,
    #        sample = sample,
    #        referenceFasta = referenceFasta,
    #        referenceFastaFai = referenceFastaFai,
    #        runDir = SVdir + 'manta/',
    #        exome = exome
    #}

    call manta_wf.MergeMantaFiles as MergeMantaFiles {
        input:
            input_manta_vcf = input_manta_vcf,
            input_manta_reference_vcfs = input_manta_reference_vcfs,
            sample_basename = sample
    }    

    call manta_wf.AnnotateMantaVCF as AnnotateMantaVCF {
        input:
            input_vcf = MergeMantaFiles.merged_vcf,
            sample_basename = sample
    }

    call manta_wf.annotSV as annotSV{
        input:
            genome_build = "GRCh37",
            input_vcf = AnnotateMantaVCF.output_manta_filtered_vcf,
            output_tsv_name = sample + ".manta.AnnotSV.tsv"
    }

    output {
        ## File mantaVcf = manta.mantaVCF
        ## File mantaVcfindex = manta.mantaVCFindex
        #File output_sv_table = AnnotateMantaVCF.output_sv_table
        #File output_manta_filtered_vcf = AnnotateMantaVCF.output_manta_filtered_vcf
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


