#!/usr/bin/env python3
"""
build_matrix_pca.py
-------------------
Equivalente Python del script R `build_matrix_pca.R`.
Se proporciona porque el sandbox de Cowork no tiene R instalado;
produce identicos archivos de salida (matriz gen-nivel + PCA).

Uso:
    python3 build_matrix_pca.py [STRINGTIE_DIR] [OUT_DIR]
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from adjustText import adjust_text

# ---- 1. Configuracion -------------------------------------------------------
STRINGTIE_DIR = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("../04_stringtie/per_sample")
OUT_DIR       = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(".")
OUT_DIR.mkdir(parents=True, exist_ok=True)

SAMPLES = ["Ch14_1","Ch14_2","Ch14_3",
           "EP08_1","EP08_2","EP08_3",
           "Se17_1","Se17_2","Se17_3",
           "SES1",  "SES2",  "SES3"]
COHORT  = {**{s:"Ch14" for s in SAMPLES[0:3]},
           **{s:"EP08" for s in SAMPLES[3:6]},
           **{s:"Se17" for s in SAMPLES[6:9]},
           **{s:"SES"  for s in SAMPLES[9:12]}}
ORIGIN  = {"Ch14":"Chapala (alta, invierno 2022)",
           "EP08":"Punta Banda (alta, 2022)",
           "Se17":"Sesma (baja, 2022)",
           "SES" :"Sesma (lote 2024, sin clasif.)"}
BATCH   = {**{s:"2022" for s in SAMPLES[0:9]},
           **{s:"2024" for s in SAMPLES[9:12]}}
PALETTE = {"Ch14":"#1f77b4","EP08":"#d62728","Se17":"#2ca02c","SES":"#9467bd"}

meta = pd.DataFrame({
    "sample": SAMPLES,
    "cohort": [COHORT[s] for s in SAMPLES],
    "origin": [ORIGIN[COHORT[s]] for s in SAMPLES],
    "batch":  [BATCH[s]  for s in SAMPLES],
})
meta.to_csv(OUT_DIR/"sample_metadata.tsv", sep="\t", index=False)
print(f"[INFO] {len(SAMPLES)} librerias | {meta['cohort'].nunique()} cohortes")

# ---- 2. Lectura modular ------------------------------------------------------
def read_gene_abund(sample: str) -> pd.DataFrame:
    f = STRINGTIE_DIR / f"{sample}.gene_abund.tab"
    if not f.exists():
        raise FileNotFoundError(f)
    df = pd.read_csv(f, sep="\t")
    # Colapsar duplicados de Gene ID (transcritos multilocus)
    agg = (df.groupby("Gene ID", as_index=False)
             .agg({"Gene Name":"first",
                   "Reference": lambda x: ",".join(sorted(set(x))),
                   "Coverage":"sum","FPKM":"sum","TPM":"sum"}))
    agg = agg.rename(columns={"Coverage": f"{sample}__Coverage",
                              "FPKM":     f"{sample}__FPKM",
                              "TPM":      f"{sample}__TPM"})
    return agg

frames = [read_gene_abund(s) for s in SAMPLES]

# ---- 3. Union outer por Gene ID ---------------------------------------------
union = frames[0]
for fr in frames[1:]:
    union = union.merge(fr, on="Gene ID", how="outer",
                        suffixes=("","_y"))
    # consolidar Gene Name / Reference (priorizar lo existente)
    for col in ("Gene Name","Reference"):
        if f"{col}_y" in union.columns:
            union[col] = union[col].combine_first(union[f"{col}_y"])
            union = union.drop(columns=[f"{col}_y"])

num_cols = [c for c in union.columns if c.endswith(("__Coverage","__FPKM","__TPM"))]
union[num_cols] = union[num_cols].fillna(0)

# ---- 4. Matrices finales -----------------------------------------------------
def make_matrix(metric: str) -> pd.DataFrame:
    cols = [f"{s}__{metric}" for s in SAMPLES]
    m = union[["Gene ID","Gene Name","Reference"] + cols].copy()
    m.columns = ["Gene ID","Gene Name","Reference"] + SAMPLES
    return m

mat_tpm  = make_matrix("TPM")
mat_fpkm = make_matrix("FPKM")
mat_cov  = make_matrix("Coverage")
mat_tpm .to_csv(OUT_DIR/"gene_matrix_TPM.tsv",      sep="\t", index=False)
mat_fpkm.to_csv(OUT_DIR/"gene_matrix_FPKM.tsv",     sep="\t", index=False)
mat_cov .to_csv(OUT_DIR/"gene_matrix_Coverage.tsv", sep="\t", index=False)
print(f"[INFO] Matriz TPM: {mat_tpm.shape[0]} genes x {len(SAMPLES)} muestras")

# ---- 5. Pre-procesamiento PCA ------------------------------------------------
expr = mat_tpm.set_index("Gene ID")[SAMPLES].astype(float)
keep = (expr >= 1).sum(axis=1) >= 3
expr_f = np.log2(expr.loc[keep] + 1)
print(f"[INFO] Genes tras filtro (TPM>=1 en >=3 muestras): {keep.sum()}")

top_n = min(5000, expr_f.shape[0])
vars_ = expr_f.var(axis=1)
expr_f = expr_f.loc[vars_.sort_values(ascending=False).index[:top_n]]
print(f"[INFO] Top variables retenidos: {expr_f.shape[0]}")

# ---- 6. PCA ------------------------------------------------------------------
X = StandardScaler().fit_transform(expr_f.T.values)
pca = PCA(n_components=min(8, X.shape[0]))
scores = pca.fit_transform(X)
var_exp = pca.explained_variance_ratio_ * 100

scores_df = pd.DataFrame(scores[:, :4],
                         columns=[f"PC{i+1}" for i in range(min(4, scores.shape[1]))],
                         index=SAMPLES).reset_index().rename(columns={"index":"sample"})
scores_df = scores_df.merge(meta, on="sample")
scores_df.to_csv(OUT_DIR/"pca_scores.tsv", sep="\t", index=False)

loadings = pd.DataFrame(pca.components_[:4].T,
                        columns=[f"PC{i+1}" for i in range(min(4, pca.components_.shape[0]))],
                        index=expr_f.index).reset_index()
loadings.to_csv(OUT_DIR/"pca_loadings.tsv", sep="\t", index=False)

# ---- 7. Visualizacion --------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 6.5))
ax.axhline(0, ls="--", color="grey", lw=0.7)
ax.axvline(0, ls="--", color="grey", lw=0.7)
markers = {"2022":"o", "2024":"^"}
texts = []
for _, r in scores_df.iterrows():
    ax.scatter(r["PC1"], r["PC2"],
               color=PALETTE[r["cohort"]],
               marker=markers[r["batch"]],
               s=120, edgecolor="black", linewidth=0.6, alpha=0.9,
               zorder=3)
    texts.append(ax.text(r["PC1"], r["PC2"], r["sample"], fontsize=9))
try:
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle="-", color="grey", lw=0.5))
except Exception:
    pass

# Leyenda de cohorte
from matplotlib.lines import Line2D
cohort_handles = [Line2D([0],[0], marker="o", color="w",
                         markerfacecolor=PALETTE[c],
                         markeredgecolor="black", markersize=10,
                         label=f"{c} - {ORIGIN[c].split(' (')[0]}")
                  for c in ["Ch14","EP08","Se17","SES"]]
batch_handles  = [Line2D([0],[0], marker="o", color="w",
                         markerfacecolor="grey", markeredgecolor="black",
                         markersize=10, label="Lote 2022"),
                  Line2D([0],[0], marker="^", color="w",
                         markerfacecolor="grey", markeredgecolor="black",
                         markersize=10, label="Lote 2024")]
leg1 = ax.legend(handles=cohort_handles, title="Cohorte",
                 loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
ax.add_artist(leg1)
ax.legend(handles=batch_handles, title="Lote",
          loc="upper left", bbox_to_anchor=(1.01, 0.55), frameon=False)

ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)")
ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)")
ax.set_title("PCA - expresion genica (StringTie TPM, log2)\n"
             f"Top {expr_f.shape[0]} genes mas variables | n = {len(SAMPLES)} librerias",
             loc="left", fontsize=12)
plt.tight_layout()
plt.savefig(OUT_DIR/"PCA_PC1_PC2.png", dpi=300, bbox_inches="tight")
plt.savefig(OUT_DIR/"PCA_PC1_PC2.pdf",            bbox_inches="tight")
plt.close()

# Scree
fig, ax = plt.subplots(figsize=(7, 4.5))
n_show = min(10, len(var_exp))
ax.bar(range(1, n_show+1), var_exp[:n_show], color="#4477AA")
for i, v in enumerate(var_exp[:n_show]):
    ax.text(i+1, v + max(var_exp)*0.02, f"{v:.1f}%",
            ha="center", fontsize=9)
ax.set_xticks(range(1, n_show+1))
ax.set_xticklabels([f"PC{i}" for i in range(1, n_show+1)])
ax.set_ylabel("% varianza")
ax.set_title("Varianza explicada por componente", loc="left")
ax.set_ylim(0, max(var_exp[:n_show]) * 1.18)
plt.tight_layout()
plt.savefig(OUT_DIR/"PCA_scree.png", dpi=300, bbox_inches="tight")
plt.close()

# ---- 8. Resumen -------------------------------------------------------------
print(f"\n[RESUMEN PCA]")
print(f"  PC1: {var_exp[0]:.2f}%   PC2: {var_exp[1]:.2f}%   PC3: {var_exp[2]:.2f}%")
print(f"\n[OUTPUTS] -> {OUT_DIR.resolve()}")
for f in sorted(OUT_DIR.iterdir()):
    print(f"  {f.name}")
