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

import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/bcftools.wdl" as bcftools
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/bwa.wdl" as bwa
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/clever.wdl" as clever
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/common.wdl" as common
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/delly.wdl" as delly
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/duphold.wdl" as duphold
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/gridss.wdl" as gridss
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/manta.wdl" as manta
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/picard.wdl" as picard
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/samtools.wdl" as samtools
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/smoove.wdl" as smoove
import "https://raw.githubusercontent.com/biowdl/tasks/2b158fa89c1541bfa8604a855607a135e26906d3/survivor.wdl" as survivor

workflow SVcalling {
    input {
        File bamFile
        File bamIndex
        File referenceFasta
        File referenceFastaFai
        File referenceFastaDict
        BwaIndex bwaIndex
        String sample
        String newId = "\'%CHROM\\_%POS\'"
        Boolean excludeMisHomRef = false
        Boolean excludeFpDupDel = false
        String outputDir = "."
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
	
    call smoove.Call as smoove {
        input:
            dockerImage = dockerImages["smoove"],
            bamFile = bamFile,
            bamIndex = bamIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            sample = sample,
            outputDir = SVdir + 'smoove'
    }   

    call delly.CallSV as delly {
        input:
            dockerImage = dockerImages["delly"],
            bamFile = bamFile,
            bamIndex = bamIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            outputPath = SVdir + 'delly/' + sample + ".delly.bcf"
    }   

    call bcftools.View as delly2vcf {
        input:
            dockerImage = dockerImages["bcftools"],
            inputFile = delly.dellyBcf,
            outputPath = SVdir + 'delly/' + sample + ".delly.vcf"
    } 

    call clever.Prediction as clever {
        input:
            dockerImage = dockerImages["clever"],
            bamFile = bamFile,
            bamIndex = bamIndex,
            bwaIndex = bwaIndex,
            outputPath = SVdir + 'clever/'
    } 
    
    call samtools.FilterShortReadsBam {
        input:
            dockerImage = dockerImages["samtools"],
            bamFile = bamFile,
            outputPathBam = SVdir + 'filteredBam/' + sample + ".filtered.bam"
    }

    call clever.Mateclever as mateclever {
        input:
            dockerImage = dockerImages["clever"],
            fiteredBam = FilterShortReadsBam.filteredBam,
            indexedFiteredBam = FilterShortReadsBam.filteredBamIndex,
            bwaIndex = bwaIndex,
            predictions = clever.predictions,
            outputPath = SVdir + 'mateclever/'
    }

    call manta.Germline as manta {
        input:
            dockerImage = dockerImages["manta"],
            bamFile = bamFile,
            bamIndex = bamIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            runDir = SVdir + 'manta/'
    }

    call gridss.GRIDSS as gridss {
        input:
            tumorBam = bamFile,
            tumorBai = bamIndex,
            tumorLabel = sample,
            reference = bwaIndex,
            outputPrefix = SVdir + 'GRIDSS/' + sample,
            dockerImage = dockerImages["gridss"]
    }

    Array[Pair[File,String]] vcfAndCaller = [(delly2vcf.outputVcf, "delly"),(manta.mantaVCF,"manta"), 
       (mateclever.matecleverVcf, "clever"),(smoove.smooveVcf,"smoove")]

    scatter (pair in vcfAndCaller) {
        String modifiedPrefix = SVdir + 'modifiedVCFs/' + sample + "." + pair.right

        call picard.RenameSample as renameSample {
            input:
                dockerImage = dockerImages["picard"],
                inputVcf = pair.left,
                outputPath =  modifiedPrefix + '.renamed.vcf',
                newSampleName = sample + "." + pair.right
        }

        call bcftools.Sort as sort {
            input:
                dockerImage = dockerImages["bcftools"],
                inputFile = renameSample.renamedVcf,
                outputPath = modifiedPrefix + '.sorted.vcf',
                tmpDir = SVdir + 'tmp-' + sample + "." + pair.right
        }

        call bcftools.Annotate as setId {
            input:
                dockerImage = dockerImages["bcftools"],
                inputFile = sort.outputVcf,
                outputPath = modifiedPrefix + '.changed_id.vcf',
                newId = newId
        }

        call duphold.Duphold as annotateDH {
            input:
            dockerImage = dockerImages["duphold"],
            inputVcf = setId.outputVcf,
            bamFile = bamFile,
            bamIndex = bamIndex,
            referenceFasta = referenceFasta,
            referenceFastaFai = referenceFastaFai,
            outputPath = modifiedPrefix + '.duphold.vcf',
            sample = sample + "." + pair.right
        }

        if (excludeMisHomRef) {
            call bcftools.View as removeMisHomRR {
                input:
                dockerImage = dockerImages["bcftools"],
                inputFile = annotateDH.outputVcf,
                outputPath = modifiedPrefix + '.MisHomRef_filtered.vcf',
                excludeUncalled = true,
                exclude = "'GT=\"RR\"'"
            }
        }

        if (excludeFpDupDel) {
            call bcftools.View as removeFpDupDel {
                input:
                    dockerImage = dockerImages["bcftools"],
                    inputFile = select_first([removeMisHomRR.outputVcf,annotateDH.outputVcf]),
                    outputPath = modifiedPrefix + '.FpDelDup_filtered.vcf',
                    include = "'(SVTYPE = \"DEL\" & FMT/DHFFC[0] < 0.7) | (SVTYPE = \"DUP\" & FMT/DHBFC[0] > 1.3)'"
            }
        }

        File toBeMergedVcfs = select_first([removeFpDupDel.outputVcf,removeMisHomRR.outputVcf, annotateDH.outputVcf])
    }

    call survivor.Merge as survivor {
        input:
            dockerImage = dockerImages["survivor"],
            filePaths = toBeMergedVcfs,
            outputPath = SVdir + 'survivor/' + sample + '.merged.vcf'
    }

    output {
        File cleverPredictions = clever.predictions
        File cleverVcf = mateclever.matecleverVcf
        File mantaVcf = manta.mantaVCF
        File dellyBcf = delly.dellyBcf
        File dellyVcf = delly2vcf.outputVcf
        File? dellyVcfIndex = delly2vcf.outputVcfIndex
        File gridssVcf = gridss.vcf
        File gridssVcfIndex = gridss.vcfIndex
        File survivorVcf = survivor.mergedVcf
        File smooveVcf = smoove.smooveVcf
        Array[File] modifiedVcfs = toBeMergedVcfs
    }

    parameter_meta {
        outputDir: {description: "The directory the output should be written to.", category: "common"}
        referenceFasta: { description: "The reference fasta file", category: "required" }
        referenceFastaFai: { description: "Fasta index (.fai) file of the reference", category: "required" }
        referenceFastaDict: { description: "Sequence dictionary (.dict) file of the reference", category: "required" }
        bamFile: {description: "sorted BAM file", category: "required"}
        bamIndex: {description: "BAM index(.bai) file", category: "required"}
        bwaIndex: {description: "Struct containing the BWA reference files", category: "required"}
        sample: {description: "The name of the sample", category: "required"}
        newId: {description: "Assign ID on the fly (e.g. --set-id +'%CHROM\_%POS').", category: "advanced"}
        excludeMisHomRef: {description: "Option to exclude missing and homozygous reference genotypes.", category: "advanced"}
        excludeFpDupDel: {description: "Option to exclude false positive duplications and deletions according to DUPHOLD.", category: "advanced"}
        dockerImages: {description: "A map describing the docker image used for the tasks.",
                           category: "advanced"}
    }
}
