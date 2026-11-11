#--------------------------------------------------------------------------------------
#
#    03PRS: PRS-CS posterior SNP weights and per-donor polygenic scores
#
#--------------------------------------------------------------------------------------

PRS = config["prs"]
RUN = PRS["run_prscs"]
CHROMS = list(range(1, 23))

# PRS-CS names its output from a, b and phi (mcmc_gtb.py); the same tag goes in every
# 03PRS filename, so changing a setting creates new files rather than reusing old ones
if RUN["phi"] is None:
    PRSCS_TAG = "a%d_b%.1f_phiauto" % (RUN["a"], RUN["b"])
    PHI_ARG   = ""
else:
    PRSCS_TAG = "a%d_b%.1f_phi%1.0e" % (RUN["a"], RUN["b"], RUN["phi"])
    PHI_ARG   = f"--phi={RUN['phi']}"

def prs_path(rule, key):
    # Fill {tag}; leave the {gwas} and {chr} wildcards in place
    return PRS[rule][key].replace("{tag}", PRSCS_TAG)

rule vcf_to_pgen:
    input:   config["genotypes"]["subset_vcf"]["output"]
    output:  PRS["vcf_to_pgen"]["output"]
    params:  prefix = PRS["vcf_to_pgen"]["prefix"]
    envmodules: "PLINK"
    resources: time = "01:00:00"
    message: "Convert the European-subset VCF to PLINK 2 format with imputation dosages"
    benchmark: PRS["vcf_to_pgen"]["benchmark"]
    log:     PRS["vcf_to_pgen"]["log"]
    shell:
             """
             plink2 --vcf {input} dosage=DS \
                    --set-missing-var-ids @:# \
                    --make-pgen \
                    --out {params.prefix} > {log} 2>&1
             """

rule run_prscs:
    input:   sumstats = config["prep_gwas"]["format_gwas"]["sumstats"],
             n_gwas   = config["prep_gwas"]["format_gwas"]["n_gwas"],
             bed      = config["genotypes"]["vcf_to_plink"]["output"],   # PRS-CS reads the .bim at this prefix
             prscs    = config["setup"]["get_prscs"]["output"],
             ld_ref   = config["setup"]["get_prscs_ref"]["output"]
    output:  prs_path("run_prscs", "output")
    wildcard_constraints: gwas = "|".join(GWAS), chr = "|".join(str(c) for c in CHROMS)
    params:  script     = config["setup"]["get_prscs"]["dir"] + "/PRScs.py",
             ref_dir    = config["setup"]["get_prscs_ref"]["ld_dir"],
             bim_prefix = config["genotypes"]["vcf_to_plink"]["prefix"],
             out_prefix = RUN["out_prefix"],
             a          = RUN["a"],
             b          = RUN["b"],
             phi_arg    = PHI_ARG,
             n_iter     = RUN["n_iter"],
             n_burnin   = RUN["n_burnin"],
             thin       = RUN["thin"],
             seed       = RUN["seed"]
    conda:   RUN["env"]
    threads: 4
    resources: mem_mb = 16000, time = "06:00:00"
    message: "PRS-CS posterior effects for {wildcards.gwas}, chromosome {wildcards.chr}"
    benchmark: prs_path("run_prscs", "benchmark")
    log:     prs_path("run_prscs", "log")
    shell:
             """
             # PRS-CS README: stop numpy using every core on the node
             export MKL_NUM_THREADS={threads} NUMEXPR_NUM_THREADS={threads} OMP_NUM_THREADS={threads}
             mkdir -p $(dirname {params.out_prefix})
             python {params.script} \
                    --ref_dir={params.ref_dir} \
                    --bim_prefix={params.bim_prefix} \
                    --sst_file={input.sumstats} \
                    --n_gwas=$(cat {input.n_gwas}) \
                    --out_dir={params.out_prefix} \
                    --chrom={wildcards.chr} \
                    --a={params.a} --b={params.b} {params.phi_arg} \
                    --n_iter={params.n_iter} --n_burnin={params.n_burnin} --thin={params.thin} \
                    --seed={params.seed} > {log} 2>&1
             """

rule combine_weights:
    input:   lambda w: expand(prs_path("run_prscs", "output"), gwas = w.gwas, chr = CHROMS)
    output:  prs_path("combine_weights", "output")
    wildcard_constraints: gwas = "|".join(GWAS)
    resources: time = "00:10:00"
    message: "Combine PRS-CS weights across chromosomes for {wildcards.gwas}"
    benchmark: prs_path("combine_weights", "benchmark")
    log:     prs_path("combine_weights", "log")
    shell:
             """
             cat {input} > {output}
             echo "SNPs with weights: $(wc -l < {output})" > {log}
             """

rule score:
    input:   weights = rules.combine_weights.output,
             pgen    = rules.vcf_to_pgen.output
    output:  prs_path("score", "output")
    wildcard_constraints: gwas = "|".join(GWAS)
    params:  pfile  = PRS["vcf_to_pgen"]["prefix"],
             prefix = prs_path("score", "prefix")
    envmodules: "PLINK"
    resources: time = "00:30:00"
    message: "Score donors on {wildcards.gwas} PRS-CS weights (dosages)"
    benchmark: prs_path("score", "benchmark")
    log:     prs_path("score", "log")
    shell:
             """
             # Weights columns: 2 = SNP ID, 4 = effect allele (A1), 6 = posterior effect
             plink2 --pfile {params.pfile} \
                    --score {input.weights} 2 4 6 cols=+scoresums \
                    --out {params.prefix} > {log} 2>&1
             """

rule combine_scores:
    input:   sscores = expand(prs_path("score", "output"), gwas = GWAS),
             pcs     = config["genotypes"]["genotype_pcs"]["output"],
             samples = config["genotypes"]["euro_samples"]["output"],
             sheet   = config["genotypes"]["sample_sheet"]
    output:  scores  = prs_path("combine_scores", "scores"),
             summary = prs_path("combine_scores", "summary")
    params:  gwas    = GWAS,
             sheet   = config["genotypes"]["euro_samples"]["sheet"],
             id_col  = config["genotypes"]["euro_samples"]["id_col"],
             sex_col = PRS["combine_scores"]["sex_col"],
             pcw_col = PRS["combine_scores"]["pcw_col"]
    singularity: config["containers"]["prs"]
    resources: time = "00:10:00"
    message: "Combine PGS across GWAS with genotype PCs, sex and PCW"
    benchmark: prs_path("combine_scores", "benchmark")
    log:     prs_path("combine_scores", "log")
    script:  PRS["combine_scores"]["script"]

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
