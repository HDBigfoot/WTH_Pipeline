#!/usr/bin/env nextflow

process ReportCleanup {

    conda 'tbvcfreport'

    publishDir params.outdir + "/tbvcfreport", mode: 'copy', saveAs: {filename -> if (filename.endsWith(".low_variants_report.txt")) {"${sampleName}.low_variants_report.txt"}
                                                                     else if (filename.endsWith(".unfixed_variants_report.txt")) {"${sampleName}.unfixed_variants_report.txt"}
                                                                     else if (filename.endsWith(".fixed_variants_report.txt")) {"${sampleName}.fixed_variants_report.txt"}}

    input:
        val sampleName
        path low_report
        path unfixed_report
        path fixed_report

    output:
        path "*.low_variants_report.txt", emit: low_report
        path "*.unfixed_variants_report.txt", emit: unfixed_report
        path "*.fixed_variants_report.txt", emit: fixed_report

    script:
    """
    tail -n +3 ${low_report} | sed 's/#//' > ${low_report}.low.tsv
    tail -n +3 ${unfixed_report} | sed 's/#//' > ${unfixed_report}.unfixed.tsv
    tail -n +3 ${fixed_report} | sed 's/#//' > ${fixed_report}.fixed.tsv

    """

}
