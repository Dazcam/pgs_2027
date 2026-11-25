#--------------------------------------------------------------------------------------
#
#    04DATA-SHARING: PGS, covariates and pseudobulk counts for the collaborators
#
#--------------------------------------------------------------------------------------

import os
import re

DS = config["data_sharing"]

def share_label(cell_type):
    # Shared name for a cell type: each cell_type_labels substitution, in order
    for old, new in DS["cell_type_labels"].items():
        cell_type = cell_type.replace(old, new)
    return cell_type

# Shared name -> pipeline name, in config order
SHARE_CELL_TYPES = {share_label(ct): ct for ct in DS["cell_types"]}
SHARE_COUNTS     = [DS["share_pseudobulk"]["output"].format(label = lab) for lab in SHARE_CELL_TYPES]

def allowed_missing(cell_type):
    # Donors allowed to be absent from a cell type's pseudobulk (pipeline name)
    donors = DS.get("allowed_missing") or {}
    return [str(d) for d, types in donors.items() if cell_type in types]

def readme_n_source(gwas):
    # Where the primary N set's n_gwas comes from, as worded in the README
    if n_source(gwas, N_PRIMARY) == "n_st1":
        return "Supplementary Table 1"
    if config["gwas"][gwas]["n_calc"]["method"] == "column":
        return "N column in file"
    return "calculated from file"

rule share_scores:
    input:   scores  = prs_path("combine_scores", "scores", N_PRIMARY),
             pcs     = config["genotypes"]["genotype_pcs"]["output"],
             samples = config["genotypes"]["euro_samples"]["output"],
             sheet   = config["genotypes"]["sample_sheet"]
    output:  scores  = DS["share_scores"]["output"]
    params:  sheet      = config["genotypes"]["euro_samples"]["sheet"],
             id_col     = config["genotypes"]["euro_samples"]["id_col"],
             sex_col    = PRS["combine_scores"]["sex_col"],
             pcw_col    = PRS["combine_scores"]["pcw_col"],
             sex_values = DS["share_scores"]["sex_values"]
    singularity: config["containers"]["prs"]
    resources: time = "00:10:00"
    message: "Check the PGS file against its sources (PCs, sex, PCW) and copy it for sharing"
    benchmark: DS["share_scores"]["benchmark"]
    log:     DS["share_scores"]["log"]
    script:  DS["share_scores"]["script"]

rule share_pseudobulk:
    input:   csv    = lambda w: DS["pseudobulk_csv"].format(cell_type = SHARE_CELL_TYPES[w.label]),
             scores = DS["share_scores"]["output"]
    output:  counts = DS["share_pseudobulk"]["output"]
    wildcard_constraints: label = "|".join(re.escape(lab) for lab in SHARE_CELL_TYPES)
    params:  cell_type       = lambda w: SHARE_CELL_TYPES[w.label],
             allowed_missing = lambda w: allowed_missing(SHARE_CELL_TYPES[w.label])
    singularity: config["containers"]["prs"]
    resources: mem_mb = 8000, time = "00:30:00"
    message: "Pseudobulk counts for {wildcards.label}: PGS donors, same order"
    benchmark: DS["share_pseudobulk"]["benchmark"]
    log:     DS["share_pseudobulk"]["log"]
    script:  DS["share_pseudobulk"]["script"]

rule share_readme:
    input:   scores  = DS["share_scores"]["output"],
             counts  = SHARE_COUNTS,
             summary = prs_path("combine_scores", "summary", N_PRIMARY),
             n_gwas  = [PRS["set_n_gwas"]["output"].format(nset = N_PRIMARY, gwas = g) for g in GWAS]
    output:  readme  = DS["share_readme"]["output"]
    params:  labels       = list(SHARE_CELL_TYPES),
             allowed_missing = {lab: allowed_missing(ct) for lab, ct in SHARE_CELL_TYPES.items()},
             gwas         = GWAS,
             gwas_files   = {g: config["gwas"][g]["file"] for g in GWAS},
             n_sources    = {g: readme_n_source(g) for g in GWAS},
             prscs        = {k: RUN[k] for k in ["a", "b", "phi", "n_iter", "n_burnin", "thin", "seed"]},
             ld_reference = DS["share_readme"]["ld_reference"],
             pca          = {"software":       DS["share_readme"]["pca_software"],
                             "maf":            config["genotypes"]["ld_prune"]["maf"],
                             "indep_pairwise": config["genotypes"]["ld_prune"]["indep_pairwise"]},
             sex_values   = DS["share_scores"]["sex_values"],
             title        = lambda w: DS["share_readme"]["title"],   # function: {n_donors} is filled in the script
             code         = DS["share_readme"]["code"]
    singularity: config["containers"]["prs"]
    resources: mem_mb = 8000, time = "00:20:00"
    message: "Check donors match across the shared files and write the README"
    benchmark: DS["share_readme"]["benchmark"]
    log:     DS["share_readme"]["log"]
    script:  DS["share_readme"]["script"]

rule share_tar:
    input:   readme = DS["share_readme"]["output"],
             scores = DS["share_scores"]["output"],
             counts = SHARE_COUNTS
    output:  tar     = DS["share_tar"]["output"],
             md5     = DS["share_tar"]["md5"],
             md5sums = DS["share_tar"]["md5sums"]
    params:  parent  = os.path.dirname(os.path.normpath(DS["dir"])),
             src     = os.path.basename(os.path.normpath(DS["dir"])),
             top_dir = DS["share_tar"]["top_dir"],
             n_files = len(SHARE_COUNTS) + 3   # pseudobulk files, PGS, README, md5sums
    resources: time = "00:30:00"
    message: "Pack the shared files into one tarball and check it"
    benchmark: DS["share_tar"]["benchmark"]
    log:     DS["share_tar"]["log"]
    shell:
             """
             share={params.parent}/{params.src}

             # Checksums of every shared file, packed with them
             (cd $share && find . -type f ! -name md5sums.txt | sed 's,^\./,,' | sort | \
                 xargs md5sum > md5sums.txt)

             tar -czf {output.tar} -C {params.parent} \
                 --transform 's,^{params.src},{params.top_dir},' {params.src} > {log} 2>&1

             # Round trip: unpack, compare checksums, test each gzip, count files
             tmp=$(mktemp -d -p {params.parent})
             trap 'rm -rf $tmp' EXIT
             tar -xzf {output.tar} -C $tmp
             (cd $tmp/{params.top_dir} && md5sum -c md5sums.txt) >> {log} 2>&1
             for f in $tmp/{params.top_dir}/pseudobulk/*.tsv.gz; do gzip -t $f; done
             n=$(find $tmp/{params.top_dir} -type f | wc -l)
             if [ "$n" -ne {params.n_files} ]; then
                 echo "Expected {params.n_files} files in the tarball, found $n" >> {log}
                 exit 1
             fi
             echo "Round trip OK: $n files, checksums match, all .tsv.gz files decompress" >> {log}

             # Checksum of the tarball itself, for checking the download
             (cd $(dirname {output.tar}) && md5sum $(basename {output.tar})) > {output.md5}
             cat {output.md5} >> {log}
             """

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
