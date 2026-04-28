#!/usr/bin/env nextflow

process SNPStatistics {

    conda 'gatk4'

    input:
        val sampleName
        path low_vcf
        path unfixed_vcf
        path fixed_vcf

    output:
        path "${low_vcf}_snpstats.tsv", emit: low_snpstats
        path "${unfixed_vcf}_snpstats.tsv", emit: unfixed_snpstats
        path "${fixed_vcf}_snpstats.tsv", emit: fixed_snpstats

    script:
    """
    gatk VariantsToTable --V ${low_vcf} --F AD --F DP --F FREQ --O ${low_vcf}_snpstats.tsv
    gatk VariantsToTable --V ${unfixed_vcf} --F AD --F DP --F FREQ --O ${unfixed_vcf}_snpstats.tsv
    gatk VariantsToTable --V ${fixed_vcf} --F AD --F DP --F FREQ --O ${fixed_vcf}_snpstats.tsv
    """

}
