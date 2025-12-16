#!/usr/bin/env nextflow

process Filtering {

    conda 'gatk4'

    publishDir params.outdir + "/VCF", mode: 'copy', saveAs: {filename -> if (filename.endsWith(".low.vcf.gz")) {"${sampleName}.low.vcf.gz"}
                                                             else if (filename.endsWith(".low.vcf.gz.tbi")) {"${sampleName}.low.vcf.gz.tbi"}
                                                             else if (filename.endsWith(".unfixed.vcf.gz")) {"${sampleName}.unfixed.vcf.gz"}
                                                             else if (filename.endsWith(".unfixed.vcf.gz.tbi")) {"${sampleName}.unfixed.vcf.gz.tbi"}
                                                             else if (filename.endsWith(".fixed.vcf.gz")) {"${sampleName}.fixed.vcf.gz"}
                                                             else if (filename.endsWith(".fixed.vcf.gz.tbi")) {"${sampleName}.fixed.vcf.gz.tbi"}}

    input:
        val sampleName
        path called_low_vcf
        path called_unfixed_vcf
        path called_fixed_vcf
        path ref
        path ref_index
        path ref_dict
        path mask
        path mask_index

    output:
        path "${called_low_vcf}.low.vcf.gz", emit: low_vcf
        path "${called_low_vcf}.low.vcf.gz.tbi", emit: low_idx
        path "${called_unfixed_vcf}.unfixed.vcf.gz", emit: unfixed_vcf
        path "${called_unfixed_vcf}.unfixed.vcf.gz.tbi", emit: unfixed_idx
        path "${called_fixed_vcf}.fixed.vcf.gz", emit: fixed_vcf
        path "${called_fixed_vcf}.fixed.vcf.gz.tbi", emit: fixed_idx

    script:
    """
    gatk VariantFiltration --R ${ref} --V ${called_low_vcf} --mask ${mask} --filter-expression "HOM > 0" --filter-name "FAILED" --O ${called_low_vcf}.flagged.snp.vcf
    gatk SelectVariants --R ${ref} --V ${called_low_vcf}.flagged.snp.vcf --exclude-filtered --O ${called_low_vcf}.low.raw.vcf
    sed 's/NC_000962.3/Chromosome/' ${called_low_vcf}.low.raw.vcf > ${called_low_vcf}.low.vcf
    bgzip ${called_low_vcf}.low.vcf
    tabix ${called_low_vcf}.low.vcf.gz
    gatk VariantFiltration --R ${ref} --V ${called_unfixed_vcf} --mask ${mask} --filter-expression "HOM > 0" --filter-name "FAILED" --O ${called_unfixed_vcf}.flagged.vSNPs.vcf
    gatk SelectVariants --R ${ref} --V ${called_unfixed_vcf}.flagged.vSNPs.vcf --exclude-filtered --O ${called_unfixed_vcf}.unfixed.raw.vcf
    sed 's/NC_000962.3/Chromosome/' ${called_unfixed_vcf}.unfixed.raw.vcf > ${called_unfixed_vcf}.unfixed.vcf
    bgzip ${called_unfixed_vcf}.unfixed.vcf
    tabix ${called_unfixed_vcf}.unfixed.vcf.gz
    gatk VariantFiltration -R ${ref} --V ${called_fixed_vcf} --mask ${mask} --O ${called_fixed_vcf}.flagged.fSNPs.vcf
    gatk SelectVariants -R ${ref} --V ${called_fixed_vcf}.flagged.fSNPs.vcf --exclude-filtered --O ${called_fixed_vcf}.fixed.raw.vcf
    sed 's/NC_000962.3/Chromosome/' ${called_fixed_vcf}.fixed.raw.vcf > ${called_fixed_vcf}.fixed.vcf
    bgzip ${called_fixed_vcf}.fixed.vcf
    tabix ${called_fixed_vcf}.fixed.vcf.gz
    """

}
