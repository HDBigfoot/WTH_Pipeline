#!/usr/bin/env nextflow

process Preprocessing {

    conda 'gatk4'

    publishDir params.outdir + "/Dedup", mode: 'copy', saveAs: { filename -> if (filename.endsWith("_dedup.bam")) {"${sampleName}_rg.bam"}
                                                               else if (filename.endsWith("_dedup.bai")) {"${sampleName}_rg.bai"}}

    input:
        val sampleName
        path bwa_aligned
        path ref
        path ref_index
        path ref_dict

    output:
        path "${bwa_aligned}_rg.bam", emit: bam_processed

    script:
    """
    gatk FixMateInformation --ASSUME_SORTED false --I ${bwa_aligned} --O ${bwa_aligned}_fixed.bam
    gatk SortSam --I ${bwa_aligned}_fixed.bam --O ${bwa_aligned}_sorted.bam --SO coordinate
    gatk AddOrReplaceReadGroups --I ${bwa_aligned}_sorted.bam --RGLB lib1 --RGPL illumina --RGPU unit1 --RGSM ${sampleName} --O ${bwa_aligned}_rg.bam
    """

}
