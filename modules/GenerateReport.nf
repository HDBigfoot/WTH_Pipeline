#!/usr/bin/env nextflow

process GenerateReport {

    conda 'tbvcfreport'

    publishDir params.outdir + "/tbvcfreport", mode: 'copy', saveAs: {filename -> if (filename.endsWith(".low_variants_report.txt")) {"${sampleName}.low_variants_report.txt"}
                                                                     else if (filename.endsWith(".unfixed_variants_report.txt")) {"${sampleName}.unfixed_variants_report.txt"}
                                                                     else if (filename.endsWith(".fixed_variants_report.txt")) {"${sampleName}.fixed_variants_report.txt"}}

    input:
        val sampleName
        path ann_low_vcf
        path ann_unfixed_vcf
        path ann_fixed_vcf

    output:
        path(*.low_variants_report.txt), emit: low_report
        path(*.unfixed_variants_report.txt), emit: unfixed_report
        path(*.fixed_variants_report.txt), emit: fixed_report

    script:
    """
    tbvcfreport generate ${ann_low_vcf}
    tbvcfreport generate ${ann_unfixed_vcf}
    tbvcfreport generate ${ann_fixed_vcf}
    """

}
