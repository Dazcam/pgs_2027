#--------------------------------------------------------------------------------------
#
#    00SETUP: containers, PRS-CS, LD reference and GWAS downloads
#
#--------------------------------------------------------------------------------------

SETUP = config["setup"]

# Downloads need internet access, so run on the login node
localrules: get_containers, get_prscs, get_prscs_ref, get_gwas_sumstats

# GWAS with a file set in config (downloaded or copied in by hand)
GWAS = [g for g, v in config["gwas"].items() if v["file"]]

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
    # download: true  -> fetch from url (with the Figshare token if token: true)
    # download: false -> file copied in by hand; check it is present
    # Both log the file size and md5, and check md5 against config if set
    # The token is read inside the shell so it never appears in params or logs
    output:  SETUP["get_gwas_sumstats"]["output"]
    wildcard_constraints: gwas = "|".join(GWAS)
    params:  url        = lambda w: config["gwas"][w.gwas]["url"] or "none (local file)",
             download   = lambda w: str(config["gwas"][w.gwas]["download"]).lower(),
             token      = lambda w: str(config["gwas"][w.gwas]["token"]).lower(),
             token_file = SETUP["get_gwas_sumstats"]["figshare_token"],
             file       = lambda w: config["gwas"][w.gwas]["file"],
             md5        = lambda w: str(config["gwas"][w.gwas]["md5"] or ""),
             dir        = SETUP["get_gwas_sumstats"]["dir"]
    message: "Get {wildcards.gwas} GWAS summary statistics"
    benchmark: SETUP["get_gwas_sumstats"]["benchmark"]
    log:     SETUP["get_gwas_sumstats"]["log"]
    shell:
             """
             (
             mkdir -p {params.dir}
             if [ "{params.download}" = "true" ] && [ "{params.token}" = "true" ]; then
                 if [ ! -s {params.token_file} ]; then
                     echo "ERROR: Figshare token file {params.token_file} missing or empty"
                     exit 1
                 fi
                 curl -fsSL -H "Authorization: token $(cat {params.token_file})" \
                      -o {params.dir}/{params.file} "{params.url}"
             elif [ "{params.download}" = "true" ]; then
                 curl -fsSL -o {params.dir}/{params.file} "{params.url}"
             elif [ ! -s {params.dir}/{params.file} ]; then
                 echo "ERROR: {params.dir}/{params.file} missing or empty; copy it in by hand"
                 echo "Source: {params.url}"
                 exit 1
             else
                 echo "Using file copied in by hand. Source: {params.url}"
             fi
             ls -lh {params.dir}/{params.file}
             md5sum {params.dir}/{params.file}
             if [ -n "{params.md5}" ]; then
                 echo "{params.md5}  {params.dir}/{params.file}" | md5sum -c -
             fi
             ) > {log} 2>&1
             touch {output}
             """

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

