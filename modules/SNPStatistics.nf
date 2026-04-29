#!/usr/bin/env nextflow

process SNPStatistics {

    conda 'gatk4'

    publishDir params.outdir + "/SNPStats", mode: 'copy', saveAs: {filename -> if (filename.endsWith("_low_snpstats.tsv")) {"${sampleName}_low_snpstats.tsv"}
                                                                  else if (filename.endsWith("_unfixed_snpstats.tsv")) {"${sampleName}_unfixed_snpstats.tsv"}
                                                                  else if (filename.endsWith("_fixed_snpstats.tsv")) {"${sampleName}_fixed_snpstats.tsv"}}

    input:
        val sampleName
        path low_vcf
        path unfixed_vcf
        path fixed_vcf

    output:
        path "${low_vcf}_low_snpstats.tsv", emit: low_snpstats
        path "${unfixed_vcf}_unfixed_snpstats.tsv", emit: unfixed_snpstats
        path "${fixed_vcf}_fixed_snpstats.tsv", emit: fixed_snpstats

    script:
    """
    gatk VariantsToTable --V ${low_vcf} --F POS --GF AD --GF DP --GF FREQ --O ${low_vcf}_low_snpstats.tsv
    gatk VariantsToTable --V ${unfixed_vcf} --F POS --GF AD --GF DP --GF FREQ --O ${unfixed_vcf}_unfixed_snpstats.tsv
    gatk VariantsToTable --V ${fixed_vcf} --F POS --GF AD --GF DP --GF FREQ --O ${fixed_vcf}_fixed_snpstats.tsv
    """

}
