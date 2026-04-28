#!/usr/bin/env nextflow

process ReportCleanup {

    publishDir params.outdir + "/Complete_Report", mode: 'copy', saveAs: {filename -> if (filename.endsWith(".complete_report.tsv")) {"${sampleName}.complete_report.tsv"}}

    input:
        val sampleName
        path low_report
        path unfixed_report
        path fixed_report
        path low_snpstats
        path unfixed_snpstats
        path fixed_snpstats

    output:
        path "${fixed_report}.complete_report.tsv", emit: complete_report

    script:
    """
    tail -n +3 ${low_report} | sed 's/#//' > ${low_report}.low.tsv
    tail -n +3 ${unfixed_report} | sed 's/#//' > ${unfixed_report}.unfixed.tsv
    tail -n +3 ${fixed_report} | sed 's/#//' > ${fixed_report}.fixed.tsv
    Rscript ${projectDir}/Scripts/merge_low.R ${low_snpstats} ${low_report}.low.tsv > ${low_report}.low.report.tsv
    Rscript ${projectDir}/Scripts/merge_unfixed.R ${unfixed_snpstats} ${unfixed_report}.unfixed.tsv > ${unfixed_report}.unfixed.final.tsv
    Rscript ${projectDir}/Scripts/merge_fixed.R ${fixed_snpstats} ${fixed_report}.fixed.tsv > ${fixed_report}.fixed.final.tsv
    cat ${fixed_report}.fixed.final.tsv ${unfixed_report}.unfixed.final.tsv ${low_report}.low.final.tsv > ${fixed_report}.complete_report.tsv
    """

}
