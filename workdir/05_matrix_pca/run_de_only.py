#!/usr/bin/env python3
"""Ejecucion streamlined del bloque DE de deseq2_pipeline.py.
QC ya generado previamente. Solo corre los contrastes B + C."""
import sys, warnings, time, json
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds  import DeseqStats
from pydeseq2.default_inference import DefaultInference

warnings.filterwarnings("ignore")
T0 = time.time()
log = lambda m: print(f"[{time.time()-T0:6.1f}s] {m}", flush=True)

OUT = Path("deseq2_out"); OUT.mkdir(exist_ok=True)
ALPHA, LFC = 0.05, 1.0
SAMPLES_ALL = ["Ch14_1","Ch14_2","Ch14_3","EP08_1","EP08_2","EP08_3",
               "Se17_1","Se17_2","Se17_3","SES1","SES2","SES3"]
counts_all = pd.read_csv("gene_count_matrix.csv", index_col=0)[SAMPLES_ALL].astype(int)
meta_all   = pd.read_csv("sample_metadata.tsv", sep="\t",
                          index_col="sample").loc[SAMPLES_ALL]

# ============================================================
# B. Contrastes 2022 (Ch14, EP08, Se17)
# ============================================================
S22 = SAMPLES_ALL[:9]
c22 = counts_all[S22]
c22 = c22.loc[c22.sum(axis=1) >= 10]
log(f"B 2022 input: {c22.shape[0]} genes x 9 muestras")

dds = DeseqDataSet(counts=c22.T, metadata=meta_all.loc[S22],
                   design="~cohort", refit_cooks=True, min_replicates=3,
                   inference=DefaultInference(n_cpus=4), quiet=True)
dds.deseq2()
log("DESeq2 fitted (2022)")

# ranking 2022 - 3 contrastes
def run_contrast(name, contrast):
    s = DeseqStats(dds, contrast=contrast, alpha=ALPHA,
                   independent_filter=True, quiet=True)
    s.summary()
    df = s.results_df.copy()
    df["gene_id"] = df.index
    df = df.sort_values("padj", na_position="last")
    df.to_csv(OUT/f"DE_{name}_all.csv", index=False)
    sig = df[(df.padj < ALPHA) & (df.log2FoldChange.abs() > LFC)]
    sig.to_csv(OUT/f"DE_{name}_significant.csv", index=False)
    # Volcano
    nlogp = -np.log10(df["padj"].clip(lower=1e-300))
    is_sig = (df.padj < ALPHA) & (df.log2FoldChange.abs() > LFC)
    up = is_sig & (df.log2FoldChange >  LFC)
    dn = is_sig & (df.log2FoldChange < -LFC)
    fig, ax = plt.subplots(figsize=(8, 6.5))
    ax.scatter(df.loc[~is_sig,"log2FoldChange"], nlogp[~is_sig],
               s=5, color="grey", alpha=0.3, label="NS")
    ax.scatter(df.loc[up,"log2FoldChange"], nlogp[up],
               s=12, color="firebrick", alpha=0.8, label=f"up ({int(up.sum())})")
    ax.scatter(df.loc[dn,"log2FoldChange"], nlogp[dn],
               s=12, color="royalblue", alpha=0.8, label=f"down ({int(dn.sum())})")
    ax.axvline(+LFC, ls="--", color="black", lw=0.5)
    ax.axvline(-LFC, ls="--", color="black", lw=0.5)
    ax.axhline(-np.log10(ALPHA), ls="--", color="black", lw=0.5)
    ax.set_xlabel("log2 fold change"); ax.set_ylabel("-log10 padj")
    ax.set_title(f"{name} | padj<{ALPHA}, |LFC|>{LFC}", loc="left")
    ax.legend(frameon=False)
    for gid, row in df[is_sig].head(10).iterrows():
        ax.annotate(gid, (row.log2FoldChange, -np.log10(row.padj)),
                    fontsize=7, alpha=0.7)
    plt.tight_layout()
    plt.savefig(OUT/f"volcano_{name}.png", dpi=200, bbox_inches="tight")
    plt.savefig(OUT/f"volcano_{name}.pdf",            bbox_inches="tight")
    plt.close()
    # ranked .rnk for GSEA
    rk = df.dropna(subset=["pvalue","log2FoldChange"]).copy()
    rk["score"] = np.sign(rk.log2FoldChange) * -np.log10(rk["pvalue"].clip(lower=1e-300))
    rk = rk.sort_values("score", ascending=False)
    rk[["gene_id","score"]].to_csv(OUT/f"ranked_{name}.rnk",
                                   sep="\t", header=False, index=False)
    return df, int(up.sum()), int(dn.sum())

CONTR = [("EP08_vs_Ch14", ["cohort","EP08","Ch14"]),
         ("Se17_vs_Ch14", ["cohort","Se17","Ch14"]),
         ("Se17_vs_EP08", ["cohort","Se17","EP08"])]
rows = []
res_all = {}
for nm, c in CONTR:
    df, nu, nd = run_contrast(nm, c)
    res_all[nm] = df
    rows.append({"contrast":nm,"n_total":len(df),"up":nu,"down":nd,
                 "n_DEG":nu+nd,"notas":"clean 2022"})
    log(f"  {nm}: total={len(df)}  up={nu}  down={nd}")

# Heatmap union top-30
top = sorted({g for d in res_all.values() for g in d.dropna(subset=["padj"]).head(30).gene_id})
dds.vst(use_design=False)
vst22 = pd.DataFrame(dds.layers["vst_counts"],
                     index=dds.obs_names, columns=dds.var_names)
sub = vst22.loc[:, [g for g in top if g in vst22.columns]]
sub = sub.sub(sub.mean(axis=0), axis=1)
fig, ax = plt.subplots(figsize=(7, max(6, 0.07*sub.shape[1])))
im = ax.imshow(sub.T.values, aspect="auto", cmap="RdBu_r", vmin=-3, vmax=3)
ax.set_xticks(range(sub.shape[0])); ax.set_xticklabels(sub.index, rotation=45, ha="right")
ax.set_yticks([])
plt.colorbar(im, ax=ax, label="VST centered")
ax.set_title(f"Top DEGs (union top-30/contraste) | n={sub.shape[1]}", loc="left")
plt.tight_layout()
plt.savefig(OUT/"heatmap_top_DEGs_2022.png", dpi=200, bbox_inches="tight")
plt.close()
log("Heatmap 2022 saved")

# ============================================================
# C. SES vs Se17 (lote confundido)
# ============================================================
S_ses = ["Se17_1","Se17_2","Se17_3","SES1","SES2","SES3"]
c_ses = counts_all[S_ses]
c_ses = c_ses.loc[c_ses.sum(axis=1) >= 10]
dds_s = DeseqDataSet(counts=c_ses.T, metadata=meta_all.loc[S_ses],
                     design="~cohort", refit_cooks=True, min_replicates=3,
                     inference=DefaultInference(n_cpus=4), quiet=True)
dds_s.deseq2()
stats = DeseqStats(dds_s, contrast=["cohort","SES","Se17"],
                   alpha=ALPHA, independent_filter=True, quiet=True)
stats.summary()
df_s = stats.results_df.copy(); df_s["gene_id"] = df_s.index
df_s = df_s.sort_values("padj", na_position="last")
df_s.to_csv(OUT/"DE_SES_vs_Se17_BATCH_CONFOUNDED.csv", index=False)
sig_s = df_s[(df_s.padj < ALPHA) & (df_s.log2FoldChange.abs() > LFC)]
rows.append({"contrast":"SES_vs_Se17_BATCH","n_total":len(df_s),
             "up":int((sig_s.log2FoldChange>0).sum()),
             "down":int((sig_s.log2FoldChange<0).sum()),
             "n_DEG":len(sig_s),
             "notas":"batch confundido - exploratorio"})
log(f"  SES_vs_Se17: total={len(df_s)}  sig={len(sig_s)}  (batch confounded)")

summary = pd.DataFrame(rows)
summary.to_csv(OUT/"DE_summary.csv", index=False)
print("\n[RESUMEN]")
print(summary.to_string(index=False))
log("DONE")
