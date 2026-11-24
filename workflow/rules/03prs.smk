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

def n_desc(gwas):
    # Short description of how per-SNP N was set for a GWAS, for the report
    n = config["gwas"][gwas]["n"]
    if n["method"] == "column":
        return f"column: {n['col']} x {n['scale']}"
    if n["method"] == "case_control":
        return f"case_control: 4 / (1/{n['cases']} + 1/{n['controls']})"
    if n["method"] == "fixed":
        return f"fixed: {n['value']}"
    return f"auto ({config['gwas'][gwas]['trait_type']})"

def report_settings():
    # group|setting|value entries joined by ; (values come from the keys the rules use)
    qc   = config["prep_gwas"]["format_gwas"]
    geno = config["genotypes"]
    rows = [
        ("PRS-CS", "PRS-CS commit", config["setup"]["get_prscs"]["commit"][:7]),
        ("PRS-CS", "LD reference", config["setup"]["get_prscs_ref"]["ld_dir"].split("/")[-1]),
        ("PRS-CS", "a", RUN["a"]),
        ("PRS-CS", "b", RUN["b"]),
        ("PRS-CS", "phi", "auto" if RUN["phi"] is None else RUN["phi"]),
        ("PRS-CS", "MCMC iterations / burn-in / thin", f"{RUN['n_iter']} / {RUN['n_burnin']} / {RUN['thin']}"),
        ("PRS-CS", "Seed", RUN["seed"]),
        ("GWAS QC", "Invalid SNP, BETA, SE or P removed", "always"),
        *[("GWAS QC", f"Filter: {name}",
           ("applied" if f.get("apply") else "counted only") +
           "".join(f", {k} {v}" for k, v in f.items() if k != "apply"))
          for name, f in qc["filters"].items()],
        ("GWAS QC", "n_gwas summary of per-SNP N", qc["n_summary"]),
        ("GWAS QC", "MAF window for auto N", "-".join(str(x) for x in qc["auto_maf"])),
        ("Donors", "Ancestry kept from sample sheet", geno["euro_samples"]["ancestry"]),
        ("Donors", "Genotype PCA MAF", geno["ld_prune"]["maf"]),
        ("Donors", "LD pruning (window step r2)", geno["ld_prune"]["indep_pairwise"]),
        ("Donors", "Genotype PCs", geno["genotype_pcs"]["pcs"]),
        ("Scoring", "Genotypes", "imputation dosages (PLINK 2, autosomes)"),
        ("Scoring", "Score", "sum of dosage x posterior effect, standardised across donors"),
    ]
    return ";".join(f"{g}|{k}|{v}" for g, k, v in rows)

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
             # Autosomes only: the PRS-CS LD reference has no chrX, and importing chrX
             # would need donor sex
             plink2 --vcf {input} dosage=DS \
                    --autosome \
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

rule prs_report:
    input:   scores    = rules.combine_scores.output.scores,
             summary   = rules.combine_scores.output.summary,
             qc        = expand(config["prep_gwas"]["format_gwas"]["qc"], gwas = GWAS),
             qc_steps  = expand(config["prep_gwas"]["format_gwas"]["qc_steps"], gwas = GWAS),
             pcs       = config["genotypes"]["genotype_pcs"]["output"],
             ancestry  = PRS["prs_report"]["ancestry_inference"],
             rmd       = PRS["prs_report"]["rmd"]
    output:  prs_path("prs_report", "html")
    params:  gwas        = ",".join(GWAS),
             n_desc      = ";".join(f"{g}|{n_desc(g)}" for g in GWAS),
             qc          = lambda w, input: ",".join(input.qc),
             qc_steps    = lambda w, input: ",".join(input.qc_steps),
             eigenval    = config["genotypes"]["genotype_pcs"]["prefix"] + ".eigenval",
             id_strip    = PRS["prs_report"]["ancestry_id_strip"],
             flag_donors = ",".join(PRS["prs_report"]["flag_donors"]),
             settings    = report_settings(),
             tag         = PRSCS_TAG,
             bmark_dir   = PRS["prs_report"]["bmark_dir"]
    singularity: config["containers"]["prs"]
    resources: mem_mb = 10000, time = "00:30:00"
    message: "Generate the PGS report"
    benchmark: prs_path("prs_report", "benchmark")
    log:     prs_path("prs_report", "log")
    shell:
             """
             Rscript -e "wd <- getwd(); rmarkdown::render('{input.rmd}', \
                 output_file = file.path(wd, '{output}'), \
                 knit_root_dir = wd, \
                 params = list(scores      = '{input.scores}', \
                               summary     = '{input.summary}', \
                               qc          = '{params.qc}', \
                               qc_steps    = '{params.qc_steps}', \
                               pcs         = '{input.pcs}', \
                               eigenval    = '{params.eigenval}', \
                               ancestry    = '{input.ancestry}', \
                               id_strip    = '{params.id_strip}', \
                               flag_donors = '{params.flag_donors}', \
                               gwas        = '{params.gwas}', \
                               n_desc      = '{params.n_desc}', \
                               settings    = '{params.settings}', \
                               tag         = '{params.tag}', \
                               bmark_dir   = '{params.bmark_dir}'))" > {log} 2>&1
             """

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
