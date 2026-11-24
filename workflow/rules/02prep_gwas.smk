#--------------------------------------------------------------------------------------
#
#    02PREP-GWAS: format and QC GWAS summary statistics for PRS-CS
#
#--------------------------------------------------------------------------------------

PREP = config["prep_gwas"]

def raw_gwas_file(wildcards):
    gwas_dir = config["setup"]["get_gwas_sumstats"]["dir"].format(gwas = wildcards.gwas)
    return f"{gwas_dir}/{config['gwas'][wildcards.gwas]['file']}"

rule format_gwas:
    input:   raw_gwas = raw_gwas_file,
             done     = config["setup"]["get_gwas_sumstats"]["output"],
             ref      = config["setup"]["get_prscs_ref"]["output"],
             bed      = config["genotypes"]["vcf_to_plink"]["output"]   # the .bim at this prefix is read
    output:  sumstats = PREP["format_gwas"]["sumstats"],
             n_gwas   = PREP["format_gwas"]["n_gwas"],
             qc       = PREP["format_gwas"]["qc"],
             qc_steps = PREP["format_gwas"]["qc_steps"]
    wildcard_constraints: gwas = "|".join(GWAS)
    params:  gwas_cfg    = lambda w: config["gwas"][w.gwas],
             filters     = PREP["format_gwas"]["filters"],
             n_summary   = PREP["format_gwas"]["n_summary"],
             auto_maf    = PREP["format_gwas"]["auto_maf"],
             ref_snpinfo = config["setup"]["get_prscs_ref"]["ld_dir"] + "/snpinfo_ukbb_hm3",
             bim         = config["genotypes"]["vcf_to_plink"]["prefix"] + ".bim"
    singularity: config["containers"]["prs"]
    resources: mem_mb = 30000, time = "02:00:00"
    message: "Format and QC {wildcards.gwas} GWAS summary statistics for PRS-CS"
    benchmark: PREP["format_gwas"]["benchmark"]
    log:     PREP["format_gwas"]["log"]
    script:  PREP["format_gwas"]["script"]

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
