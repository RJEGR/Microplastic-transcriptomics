#!/usr/bin/env python3
"""
prepDE_counts.py
----------------
Genera gene_count_matrix.csv y transcript_count_matrix.csv compatibles con
DESeq2/edgeR a partir de los GTF per-sample de StringTie.

Equivalente a prepDE.py3 (Pertea lab) PERO sin requerir el flag `-e` de
StringTie (los GTF de este proyecto se generaron en modo assembly).
La conversion cov -> count usa la misma formula oficial:

        count_transcrito = ceil( cov * longitud_efectiva / read_length )

donde longitud_efectiva = suma de exones del transcrito.
El conteo gen-nivel es la suma de sus transcritos.

Uso:
    python3 prepDE_counts.py <sample_list.tsv> <out_dir> [--read-length 150]

sample_list.tsv: dos columnas TAB-separadas
    sample_id   /path/al/sample.gtf

Autor: R. Gomez (UABC) - 2026-06-24
"""
from __future__ import annotations
import argparse, csv, re, sys
from collections import defaultdict
from math import ceil
from pathlib import Path

RE_GENE_ID       = re.compile(r'gene_id "([^"]+)"')
RE_GENE_NAME     = re.compile(r'gene_name "([^"]+)"')
RE_TRANSCRIPT_ID = re.compile(r'transcript_id "([^"]+)"')
RE_COV           = re.compile(r'cov "([\-\+\d\.]+)"')

def get_gene_id(attrs: str, tid: str) -> str:
    g  = RE_GENE_ID.search(attrs)
    gn = RE_GENE_NAME.search(attrs)
    if g:
        return f"{g.group(1)}|{gn.group(1)}" if gn else g.group(1)
    return tid

def parse_gtf(path: Path, read_len: int):
    """Devuelve dict transcript_id -> (gene_id, count)."""
    t2g, t2cnt = {}, {}
    cur_tid = cur_gid = None
    cur_cov = 0.0
    cur_len = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            v = line.rstrip("\n").split("\t")
            if len(v) < 9:
                continue
            ftype, attrs = v[2], v[8]
            if ftype == "transcript":
                # cerrar transcrito previo
                if cur_tid is not None and cur_len > 0:
                    t2cnt[cur_tid] = int(ceil(cur_cov * cur_len / read_len))
                cur_tid = RE_TRANSCRIPT_ID.search(attrs).group(1)
                cur_gid = get_gene_id(attrs, cur_tid)
                m = RE_COV.search(attrs)
                cur_cov = max(0.0, float(m.group(1))) if m else 0.0
                cur_len = 0
                t2g[cur_tid] = cur_gid
            elif ftype == "exon":
                cur_len += int(v[4]) - int(v[3]) + 1
    if cur_tid is not None and cur_len > 0:
        t2cnt[cur_tid] = int(ceil(cur_cov * cur_len / read_len))
    return t2g, t2cnt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sample_list", type=Path,
                    help="TSV con columnas: sample_id <TAB> gtf_path")
    ap.add_argument("out_dir", type=Path)
    ap.add_argument("--read-length", type=int, default=150,
                    help="Longitud promedio de lectura PE (default 150)")
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    samples = []
    for line in args.sample_list.read_text().splitlines():
        if not line.strip() or line.startswith("#"): continue
        sid, gtf = line.split("\t")
        samples.append((sid, Path(gtf)))
    samples.sort()
    print(f"[INFO] {len(samples)} muestras | read_len = {args.read_length}",
          file=sys.stderr)

    # transcript_dict[tid][sid] = count
    t_dict = defaultdict(lambda: defaultdict(int))
    t2g_global = {}

    for sid, gtf in samples:
        print(f"  -> parsing {sid}", file=sys.stderr)
        t2g, t2cnt = parse_gtf(gtf, args.read_length)
        t2g_global.update(t2g)
        for tid, cnt in t2cnt.items():
            t_dict[tid][sid] = cnt

    # ---- transcript matrix ---------------------------------------------------
    tx_csv = args.out_dir / "transcript_count_matrix.csv"
    with tx_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["transcript_id"] + [s for s,_ in samples])
        for tid in sorted(t_dict):
            w.writerow([tid] + [t_dict[tid].get(s, 0) for s,_ in samples])
    print(f"[OK] {tx_csv} ({len(t_dict)} transcritos)", file=sys.stderr)

    # ---- gene matrix (suma de transcritos por gen) ---------------------------
    g_dict = defaultdict(lambda: defaultdict(int))
    for tid, sample_counts in t_dict.items():
        gid = t2g_global[tid]
        for s, c in sample_counts.items():
            g_dict[gid][s] += c

    gene_csv = args.out_dir / "gene_count_matrix.csv"
    with gene_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene_id"] + [s for s,_ in samples])
        for gid in sorted(g_dict):
            w.writerow([gid] + [g_dict[gid].get(s, 0) for s,_ in samples])
    print(f"[OK] {gene_csv} ({len(g_dict)} genes)", file=sys.stderr)

if __name__ == "__main__":
    main()
