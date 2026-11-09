#--------------------------------------------------------------------------------------
#
#    00SETUP: containers, PRS-CS, LD reference and GWAS downloads
#
#--------------------------------------------------------------------------------------

SETUP = config["setup"]

# Downloads need internet access, so run on the login node
localrules: get_containers, get_prscs, get_prscs_ref, get_gwas_sumstats

# GWAS with a public download link (collaborator files are copied by hand)
GWAS_DOWNLOAD = [g for g, v in config["gwas"].items() if v["url"]]

# Container <target.sif> <uri> pairs for setup_pull_containers.sh
CONTAINER_PAIRS = " ".join(f"{config['containers'][name]} {uri}"
                           for name, uri in SETUP["get_containers"]["uris"].items())

rule get_containers:
    output:  SETUP["get_containers"]["output"]
    params:  pairs  = CONTAINER_PAIRS,
             script = SETUP["get_containers"]["script"]
    message: "Pull Singularity / Apptainer containers from quay.io"
    benchmark: SETUP["get_containers"]["benchmark"]
    log:     SETUP["get_containers"]["log"]
    shell:
             """
             bash {params.script} {params.pairs} > {log} 2>&1
             touch {output}
             """

rule get_prscs:
    output:  SETUP["get_prscs"]["output"]
    params:  repo   = SETUP["get_prscs"]["repo"],
             commit = SETUP["get_prscs"]["commit"],
             dir    = SETUP["get_prscs"]["dir"]
    message: "Clone PRS-CS and check out the pinned commit"
    benchmark: SETUP["get_prscs"]["benchmark"]
    log:     SETUP["get_prscs"]["log"]
    shell:
             """
             (
             [ -d {params.dir}/.git ] || git clone {params.repo} {params.dir}
             git -C {params.dir} fetch origin
             git -C {params.dir} checkout {params.commit}
             git -C {params.dir} log -1 --format="Commit: %H (%cd)"
             ) > {log} 2>&1
             touch {output}
             """

rule get_prscs_ref:
    output:  SETUP["get_prscs_ref"]["output"]
    params:  url     = SETUP["get_prscs_ref"]["url"],
             tarball = SETUP["get_prscs_ref"]["tarball"],
             dir     = SETUP["get_prscs_ref"]["dir"]
    message: "Download and extract the PRS-CS LD reference panel"
    benchmark: SETUP["get_prscs_ref"]["benchmark"]
    log:     SETUP["get_prscs_ref"]["log"]
    shell:
             """
             (
             mkdir -p {params.dir}
             wget -O {params.tarball} "{params.url}"
             tar -zxvf {params.tarball} -C {params.dir}
             rm {params.tarball}
             ) > {log} 2>&1
             touch {output}
             """

rule get_gwas_sumstats:
    output:  SETUP["get_gwas_sumstats"]["output"]
    wildcard_constraints: gwas = "|".join(GWAS_DOWNLOAD)
    params:  url  = lambda w: config["gwas"][w.gwas]["url"],
             file = lambda w: config["gwas"][w.gwas]["file"],
             md5  = lambda w: str(config["gwas"][w.gwas]["md5"] or ""),
             dir  = SETUP["get_gwas_sumstats"]["dir"],
             user_agent = SETUP["get_gwas_sumstats"]["user_agent"]
    message: "Download {wildcards.gwas} GWAS summary statistics"
    benchmark: SETUP["get_gwas_sumstats"]["benchmark"]
    log:     SETUP["get_gwas_sumstats"]["log"]
    shell:
             """
             (
             mkdir -p {params.dir}
             wget -U "{params.user_agent}" -O {params.dir}/{params.file} "{params.url}"
             if [ -n "{params.md5}" ]; then
                 echo "{params.md5}  {params.dir}/{params.file}" | md5sum -c -
             fi
             ls -lh {params.dir}/{params.file}
             ) > {log} 2>&1
             touch {output}
             """

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

