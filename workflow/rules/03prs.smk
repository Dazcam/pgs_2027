#--------------------------------------------------------------------------------------
#
#    03PRS: PRS-CS posterior SNP weights and per-donor polygenic scores
#
#--------------------------------------------------------------------------------------

PRS = config["prs"]

localrules: set_n_gwas   # writes one number; not worth a SLURM job
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

# N sets for --n_gwas: n_st1 (Supplementary Table 1, primary) and n_calc (from the file)
N_SETS    = PRS["n_sets"]
N_PRIMARY = PRS["n_primary"]

def prs_path(rule, key, nset = None):
    # Fill {tag}, and {nset} if given; leave the {gwas} and {chr} wildcards in place
    path = PRS[rule][key].replace("{tag}", PRSCS_TAG)
    return path if nset is None else path.replace("{nset}", nset)

def n_source(gwas, nset):
    # Config key that gives the N for a GWAS in an N set: n_st1 falls back to n_calc
    if nset == "n_st1" and "n_st1" in config["gwas"][gwas]:
        return "n_st1"
    return "n_calc"

def nset_gwas(nset):
    # GWAS run in an N set: all of them in the primary set; in the others, only those
    # whose N differs from the primary set (the rest would repeat it exactly)
    if nset == N_PRIMARY:
        return GWAS
    return [g for g in GWAS if n_source(g, nset) != n_source(g, N_PRIMARY)]

def n_value_desc(gwas, nset):
    # Where the N in an N set comes from, for the log and the report
    if n_source(gwas, nset) == "n_st1":
        return "Supplementary Table 1"
    return "file: " + n_desc(gwas)

def n_desc(gwas):
    # Short description of how per-SNP N is calculated from a GWAS file (n_calc)
    n = config["gwas"][gwas]["n_calc"]
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
        ("N for PRS-CS", "N sets run", ", ".join(N_SETS)),
        ("N for PRS-CS", "Primary N set (deliverable)", N_PRIMARY),
        ("GWAS QC", "n_gwas summary of per-SNP N (n_calc)", qc["n_summary"]),
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

rule set_n_gwas:
    input:   n_calc = config["prep_gwas"]["format_gwas"]["n_gwas"]
    output:  PRS["set_n_gwas"]["output"]
    wildcard_constraints: gwas = "|".join(GWAS), nset = "|".join(N_SETS)
    params:  source = lambda w: n_source(w.gwas, w.nset),
             n_st1  = lambda w: config["gwas"][w.gwas].get("n_st1", ""),
             desc   = lambda w: n_value_desc(w.gwas, w.nset)
    message: "Set n_gwas for {wildcards.gwas} in N set {wildcards.nset}"
    benchmark: PRS["set_n_gwas"]["benchmark"]
    log:     PRS["set_n_gwas"]["log"]
    shell:
             """
             if [ "{params.source}" = "n_st1" ]; then
                 echo {params.n_st1} > {output}
             else
                 cp {input.n_calc} {output}
             fi
             echo "{wildcards.gwas} {wildcards.nset}: n_gwas $(cat {output}) ({params.desc})" > {log}
             """

rule run_prscs:
    input:   sumstats = config["prep_gwas"]["format_gwas"]["sumstats"],
             n_gwas   = PRS["set_n_gwas"]["output"],
             bed      = config["genotypes"]["vcf_to_plink"]["output"],   # PRS-CS reads the .bim at this prefix
             prscs    = config["setup"]["get_prscs"]["output"],
             ld_ref   = config["setup"]["get_prscs_ref"]["output"]
    output:  prs_path("run_prscs", "output")
    wildcard_constraints: gwas = "|".join(GWAS), chr = "|".join(str(c) for c in CHROMS),
                          nset = "|".join(N_SETS)
    params:  script     = config["setup"]["get_prscs"]["dir"] + "/PRScs.py",
             ref_dir    = config["setup"]["get_prscs_ref"]["ld_dir"],
             bim_prefix = config["genotypes"]["vcf_to_plink"]["prefix"],
             out_prefix = lambda w: RUN["out_prefix"].format(nset = w.nset, gwas = w.gwas),
             a          = RUN["a"],
             b          = RUN["b"],
             phi_arg    = PHI_ARG,
             n_iter     = RUN["n_iter"],
             n_burnin   = RUN["n_burnin"],
             thin       = RUN["thin"],
             seed       = RUN["seed"]
    conda:   RUN["env"]
    threads: 2
    resources: mem_mb = 16000, time = "06:00:00", threads = 2  # 2 CPUs x 90 jobs (as > 300 jobs runs using maxjobs500_rome)
    message: "PRS-CS posterior effects for {wildcards.gwas} ({wildcards.nset}), chromosome {wildcards.chr}"
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
    input:   lambda w: expand(prs_path("run_prscs", "output", w.nset), gwas = w.gwas, chr = CHROMS)
    output:  prs_path("combine_weights", "output")
    wildcard_constraints: gwas = "|".join(GWAS), nset = "|".join(N_SETS)
    resources: time = "00:10:00"
    message: "Combine PRS-CS weights across chromosomes for {wildcards.gwas} ({wildcards.nset})"
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
    wildcard_constraints: gwas = "|".join(GWAS), nset = "|".join(N_SETS)
    params:  pfile  = PRS["vcf_to_pgen"]["prefix"],
             prefix = lambda w: prs_path("score", "prefix", w.nset).format(gwas = w.gwas)
    envmodules: "PLINK"
    resources: time = "00:30:00"
    message: "Score donors on {wildcards.gwas} PRS-CS weights ({wildcards.nset}, dosages)"
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
    input:   sscores = lambda w: expand(prs_path("score", "output", w.nset), gwas = nset_gwas(w.nset)),
             pcs     = config["genotypes"]["genotype_pcs"]["output"],
             samples = config["genotypes"]["euro_samples"]["output"],
             sheet   = config["genotypes"]["sample_sheet"]
    output:  scores  = prs_path("combine_scores", "scores"),
             summary = prs_path("combine_scores", "summary")
    wildcard_constraints: nset = "|".join(N_SETS)
    params:  gwas    = lambda w: nset_gwas(w.nset),
             sheet   = config["genotypes"]["euro_samples"]["sheet"],
             id_col  = config["genotypes"]["euro_samples"]["id_col"],
             sex_col = PRS["combine_scores"]["sex_col"],
             pcw_col = PRS["combine_scores"]["pcw_col"]
    singularity: config["containers"]["prs"]
    resources: time = "00:10:00"
    message: "Combine PGS across GWAS with genotype PCs, sex and PCW ({wildcards.nset})"
    benchmark: prs_path("combine_scores", "benchmark")
    log:     prs_path("combine_scores", "log")
    script:  PRS["combine_scores"]["script"]

rule prs_report:
    input:   scores    = prs_path("combine_scores", "scores", N_PRIMARY),
             summary   = prs_path("combine_scores", "summary", N_PRIMARY),
             set_scores = [prs_path("combine_scores", "scores", n) for n in N_SETS],
             n_gwas    = [PRS["set_n_gwas"]["output"].format(nset = n, gwas = g)
                          for n in N_SETS for g in nset_gwas(n)],
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
             n_primary   = N_PRIMARY,
             set_scores  = ";".join(f"{n}|{prs_path('combine_scores', 'scores', n)}" for n in N_SETS),
             n_gwas      = ";".join(f"{n}|{g}|{PRS['set_n_gwas']['output'].format(nset = n, gwas = g)}|{n_value_desc(g, n)}"
                                    for n in N_SETS for g in nset_gwas(n)),
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
                               n_primary   = '{params.n_primary}', \
                               set_scores  = '{params.set_scores}', \
                               n_gwas      = '{params.n_gwas}', \
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
