#--------------------------------------------------------------------------------------
#
#    01GENOTYPES: European subset, PLINK conversion and genotype PCs
#
#--------------------------------------------------------------------------------------

GENO = config["genotypes"]

rule euro_samples:
    input:   GENO["sample_sheet"]
    output:  GENO["euro_samples"]["output"]
    params:  sheet        = GENO["euro_samples"]["sheet"],
             id_col       = GENO["euro_samples"]["id_col"],
             ancestry_col = GENO["euro_samples"]["ancestry_col"],
             ancestry     = GENO["euro_samples"]["ancestry"]
    singularity: config["containers"]["prs"]
    resources: time = "00:10:00"
    message: "Extract European-ancestry sample IDs from the sample sheet"
    benchmark: GENO["euro_samples"]["benchmark"]
    log:     GENO["euro_samples"]["log"]
    script:  GENO["euro_samples"]["script"]

rule subset_vcf:
    input:   vcf     = GENO["vcf"],
             samples = rules.euro_samples.output
    output:  GENO["subset_vcf"]["output"]
    params:  md5 = GENO["vcf_md5"]
    envmodules: "BCFtools"
    resources: time = "01:00:00"
    message: "Check the VCF md5 and subset it to European-ancestry samples"
    benchmark: GENO["subset_vcf"]["benchmark"]
    log:     GENO["subset_vcf"]["log"]
    shell:
             """
             (
             echo "{params.md5}  {input.vcf}" | md5sum -c -
             bcftools view --samples-file {input.samples} --output-type z --output {output} {input.vcf}
             bcftools index --tbi {output}
             echo "Samples in subset: $(bcftools query -l {output} | wc -l)"
             echo "Variants in subset: $(bcftools index --nrecords {output})"
             ) > {log} 2>&1
             """

rule vcf_to_plink:
    input:   rules.subset_vcf.output
    output:  GENO["vcf_to_plink"]["output"]
    params:  prefix = GENO["vcf_to_plink"]["prefix"]
    singularity: config["containers"]["prs"]
    resources: time = "01:00:00"
    message: "Convert the European-subset VCF to PLINK format"
    benchmark: GENO["vcf_to_plink"]["benchmark"]
    log:     GENO["vcf_to_plink"]["log"]
    shell:
             """
             plink --vcf {input} \
                   --double-id \
                   --set-missing-var-ids @:# \
                   --make-bed \
                   --out {params.prefix} > {log} 2>&1
             """

rule ld_prune:
    input:   rules.vcf_to_plink.output
    output:  GENO["ld_prune"]["output"]
    params:  bfile          = GENO["vcf_to_plink"]["prefix"],
             maf            = GENO["ld_prune"]["maf"],
             indep_pairwise = GENO["ld_prune"]["indep_pairwise"],
             prefix         = GENO["ld_prune"]["prefix"]
    singularity: config["containers"]["prs"]
    resources: time = "01:00:00"
    message: "LD prune common SNPs before the genotype PCA"
    benchmark: GENO["ld_prune"]["benchmark"]
    log:     GENO["ld_prune"]["log"]
    shell:
             """
             plink --bfile {params.bfile} \
                   --maf {params.maf} \
                   --indep-pairwise {params.indep_pairwise} \
                   --out {params.prefix} > {log} 2>&1
             """

rule genotype_pcs:
    input:   bed     = rules.vcf_to_plink.output,
             prune   = rules.ld_prune.output
    output:  GENO["genotype_pcs"]["output"]
    params:  bfile  = GENO["vcf_to_plink"]["prefix"],
             pcs    = GENO["genotype_pcs"]["pcs"],
             prefix = GENO["genotype_pcs"]["prefix"]
    singularity: config["containers"]["prs"]
    resources: time = "01:00:00"
    message: "Genotype PCA on LD-pruned SNPs in European-ancestry samples"
    benchmark: GENO["genotype_pcs"]["benchmark"]
    log:     GENO["genotype_pcs"]["log"]
    shell:
             """
             plink --bfile {params.bfile} \
                   --extract {input.prune} \
                   --pca {params.pcs} \
                   --out {params.prefix} > {log} 2>&1
             """

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
