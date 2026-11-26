#--------------------------------------------------------------------------------------
#
#    Nucleus-level raw counts with cell type annotation, for sharing (cNMF)
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  04DATA-SHARING | Rule: share_nuclei_counts
#            Upstream:   eQTL_study_2025 02SCANPY L1 object and L2 subcluster objects;
#                        share_scores
#            Downstream: share_readme, share_tar
#
# Purpose:   - Nuclei from the donors in the shared PGS file, in L1 object order
#            - L1 cell type from the clustering (cluster_col via cluster_anns); L2 from
#              each subcluster object's subcluster column, joined by nucleus ID
#            - PCW and sex from the shared PGS file (one value per donor)
#            - A new AnnData with only: X = raw counts (counts layer, sparse); obs =
#              iid, cell_type_l1, cell_type_l2, PCW, sex (index cell_id, the source
#              nucleus IDs); var = gene IDs. Written gzip-compressed
#
# Checks:    Nucleus IDs unique; every kept nucleus has an L1 type; every L2 nucleus is
#            in the L1 object, in no other L2 object, and of the matching L1 type;
#            every donor has nuclei; sex in the object agrees with the PGS file;
#            counts are non-negative whole numbers
#
# Inputs:    l1       L1 h5ad (obs, var and the counts layer are read; X is not)
#            l2       L2 h5ad files (obs only), one per subclustered L1 type
#            scores   shared pgs_scores.tsv (donors, PCW, sex)
#
# Outputs:   h5ad     nuclei_counts.h5ad
#            summary  nuclei per L1 and L2 cell type, for the README
#
#--------------------------------------------------------------------------------------

import sys
import h5py
import numpy as np
import pandas as pd
import anndata as ad

try:
    from anndata.io import read_elem
except ImportError:
    from anndata.experimental import read_elem

##  Logging  --------------------------------------------------------------------------
log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log
print("\n\nPreparing nucleus-level counts for sharing ...")

##  Variables  ------------------------------------------------------------------------
l1_in        = snakemake.input["l1"]
l2_in        = list(snakemake.input["l2"])
scores_in    = snakemake.input["scores"]
h5ad_out     = snakemake.output["h5ad"]
summary_out  = snakemake.output["summary"]

l2_types     = list(snakemake.params["l2_types"])          # L1 type of each L2 file, same order
cluster_col  = snakemake.params["cluster_col"]
cluster_anns = {str(k): v for k, v in snakemake.params["cluster_anns"].items()}
sample_col   = snakemake.params["sample_col"]
sex_col      = snakemake.params["sex_col"]
counts_layer = snakemake.params["counts_layer"]
labels       = dict(snakemake.params["labels"])            # pipeline -> shared cell type names

print("\nVariables")
print("============================")
for name, value in [("l1_in", l1_in), ("l2_in", ", ".join(l2_in)), ("scores_in", scores_in),
                    ("h5ad_out", h5ad_out), ("summary_out", summary_out),
                    ("l2_types", ", ".join(l2_types)), ("cluster_col", cluster_col),
                    ("cluster_anns", cluster_anns), ("sample_col", sample_col),
                    ("sex_col", sex_col), ("counts_layer", counts_layer), ("labels", labels)]:
    print(f"{name:13s} {value}")
print("============================\n")


def relabel(x):
    for old, new in labels.items():
        x = x.replace(old, new)
    return x


def fail(msg):
    print(f"\nERROR: {msg}")
    log.flush()
    raise SystemExit(msg)


##  Donors  ---------------------------------------------------------------------------
scores = pd.read_csv(scores_in, sep="\t", dtype={"iid": str})
donors = scores.set_index("iid")[["PCW", "sex"]]
print(f"Donors in the PGS file: {len(donors)}")

##  L1 object: obs, var  --------------------------------------------------------------
with h5py.File(l1_in, "r") as f:
    obs = read_elem(f["obs"])
    var = read_elem(f["var"])
print(f"L1 object: {obs.shape[0]:,} nuclei, {var.shape[0]:,} genes")

if not obs.index.is_unique:
    fail(f"{obs.index.duplicated().sum()} duplicated nucleus IDs in the L1 object")

obs[sample_col] = obs[sample_col].astype(str)
keep = obs[sample_col].isin(donors.index).to_numpy()
print(f"Nuclei from the PGS donors: {keep.sum():,} of {len(obs):,}")

missing_donors = sorted(set(donors.index) - set(obs.loc[keep, sample_col]))
if missing_donors:
    fail(f"PGS donors with no nuclei: {', '.join(missing_donors)}")

##  Annotation  -----------------------------------------------------------------------
out = pd.DataFrame(index=obs.index[keep].astype(str))
out.index.name = "cell_id"
out["iid"] = obs.loc[keep, sample_col].to_numpy()

l1 = obs.loc[keep, cluster_col].astype(str).map(cluster_anns)
if l1.isna().any():
    bad = sorted(obs.loc[keep, cluster_col].astype(str)[l1.isna()].unique())
    fail(f"{cluster_col} clusters with no entry in cluster_anns: {', '.join(bad)}")
out["cell_type_l1"] = l1.to_numpy()
out["cell_type_l2"] = pd.NA

seen = pd.Series(dtype=str)
for path, l1_type in zip(l2_in, l2_types):
    with h5py.File(path, "r") as f:
        sub = read_elem(f["obs"])["subcluster"].astype(str)
    sub.index = sub.index.astype(str)
    not_in_l1 = ~sub.index.isin(obs.index.astype(str))
    if not_in_l1.any():
        fail(f"{path}: {not_in_l1.sum()} nuclei not in the L1 object")
    dup = sub.index[sub.index.isin(seen.index)]
    if len(dup) > 0:
        fail(f"{path}: {len(dup)} nuclei also in {seen[dup[0]]}")
    seen = pd.concat([seen, pd.Series(path, index=sub.index)])

    sub = sub[sub.index.isin(out.index)]
    wrong = out.loc[sub.index, "cell_type_l1"] != l1_type
    if wrong.any():
        fail(f"{path}: {wrong.sum()} nuclei are not {l1_type} at L1")
    out.loc[sub.index, "cell_type_l2"] = sub.to_numpy()
    print(f"{l1_type}: {len(sub):,} PGS-donor nuclei with an L2 label "
          f"({', '.join(sorted(sub.unique()))})")

no_l2 = out["cell_type_l1"].isin(l2_types) & out["cell_type_l2"].isna()
print(f"Nuclei of subclustered L1 types without an L2 label: {no_l2.sum():,}")

##  Donor covariates  -----------------------------------------------------------------
obj_sex = obs.loc[keep, sex_col].astype(str).to_numpy()
pgs_sex = donors.loc[out["iid"], "sex"].astype(str).to_numpy()
if (obj_sex != pgs_sex).any():
    bad = sorted(set(out["iid"][obj_sex != pgs_sex]))
    fail(f"Sex in the object differs from the PGS file for donors: {', '.join(bad)}")
out["PCW"] = donors.loc[out["iid"], "PCW"].to_numpy()
out["sex"] = pgs_sex

for col in ["cell_type_l1", "cell_type_l2"]:
    out[col] = out[col].astype("string").map(lambda x: relabel(x) if pd.notna(x) else x)
for col in ["iid", "cell_type_l1", "cell_type_l2", "sex"]:
    out[col] = out[col].astype("category")

##  Counts  ---------------------------------------------------------------------------
print(f"\nReading layers['{counts_layer}'] ...")
with h5py.File(l1_in, "r") as f:
    counts = read_elem(f["layers"][counts_layer])
counts = counts[np.flatnonzero(keep)]
counts.sort_indices()

if counts.data.size and (counts.data.min() < 0 or not np.all(np.mod(counts.data, 1) == 0)):
    fail("Counts layer has negative or non-whole values")
print(f"Counts: {counts.shape[0]:,} nuclei x {counts.shape[1]:,} genes, "
      f"{counts.nnz:,} non-zero values, {counts.dtype}")

##  Write  ----------------------------------------------------------------------------
adata = ad.AnnData(X=counts, obs=out, var=pd.DataFrame(index=var.index.astype(str)))
adata.var.index.name = "gene"
adata.write_h5ad(h5ad_out, compression="gzip")
print(f"\nWritten: {h5ad_out}")

summary = pd.concat([
    out.groupby("cell_type_l1", observed=True).size().rename("nuclei").reset_index()
       .rename(columns={"cell_type_l1": "cell_type"}).assign(level="L1"),
    out.groupby("cell_type_l2", observed=True).size().rename("nuclei").reset_index()
       .rename(columns={"cell_type_l2": "cell_type"}).assign(level="L2"),
    pd.DataFrame({"level": ["all"], "cell_type": ["nuclei"], "nuclei": [len(out)]}),
    pd.DataFrame({"level": ["all"], "cell_type": ["genes"], "nuclei": [counts.shape[1]]})])
summary[["level", "cell_type", "nuclei"]].to_csv(summary_out, sep="\t", index=False)
print(summary.to_string(index=False))
print("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
