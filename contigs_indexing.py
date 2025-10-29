#!/usr/bin/env python3
import argparse
import csv
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

BASE_HEADER = [
    "query_name", "subject_name", "identity", "length",
    "mismatch", "gap", "query_start", "query_end",
    "subject_start", "subject_end", "E_value",
    "bit_score", "gene_length",
]
FINAL_HEADER = [
    "query_name","subject_name","identity","length","mismatch","gap",
    "query_start","query_end","subject_start","subject_end",
    "genomic_start","genomic_end","E_value","bit_score","gene_length",
    "ref_scaff","scaff_start","scaff_stop","orientation","contig_length"
]
IDX_LENGTH = 3
IDX_EVALUE = 10
IDX_GENELEN = 12

def read_tsv_13(path: Path):
    rows = []
    with path.open("r", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        for rec in r:
            if len(rec) != 13:
                raise ValueError(f"{path}: expected 13 columns, got {len(rec)}")
            rows.append(rec)
    return rows

def load_parameters(params_path: Path):
    # parameters_dScaff.txt: header + 8 lines like the ones you sent
    with params_path.open("r", encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    body = lines[1:]
    if len(body) < 8:
        raise ValueError("parameters file must have at least 8 rows after header")

    def third_token(line):
        parts = line.split()
        if len(parts) < 3:
            raise ValueError(f"Bad parameter line: {line}")
        return parts[2]

    parm1 = float(third_token(body[0]))
    parm2 = float(third_token(body[1]))
    parm6 = float(third_token(body[5]))
    parm7 = float(third_token(body[6]))
    parm8 = float(third_token(body[7]))

    def special(line):
        tok = third_token(line)        # e.g. "0,5"
        parts = tok.split(",")
        if len(parts) < 2:
            raise ValueError(f"Expected a comma in 3rd token: {line}")
        return float(parts[1]) / 10.0  # matches your R: /10

    parm3 = special(body[2])
    parm4 = special(body[3])
    parm5 = special(body[4])

    return parm1, parm2, parm3, parm4, parm5, parm6, parm7, parm8

def load_contig_lengths(path: Path):
    contig_len = {}
    with path.open("r", newline="") as f:
        r = csv.reader(f, delimiter=",")
        for rec in r:
            if not rec:
                continue
            name = rec[0]
            try:
                val = float(rec[1])
            except ValueError:
                val = float("nan")
            contig_len[name] = val
    return contig_len

def discover_subdirs(main_dir: Path):
    return [p for p in main_dir.iterdir() if p.is_dir()]

def list_query_csvs(subdir: Path):
    q = subdir / "queries"
    if not q.is_dir():
        return []
    # Non-recursive, like your R `list.files()` in that folder
    return sorted([p for p in q.iterdir() if p.is_file() and p.suffix.lower() == ".csv"])

def filter_one_file(file_path: str, parm2: float, parm3: float, parm4: float, parm5: float):
    p = Path(file_path)
    try:
        rows = read_tsv_13(p)
        if not rows:
            return []
        try:
            gl = float(rows[0][IDX_GENELEN])
        except ValueError:
            return []
        if gl <= parm2:
            return []

        def pick(mult: float):
            thr = gl * mult
            out = []
            for r in rows:
                try:
                    lv = float(r[IDX_LENGTH])
                except ValueError:
                    continue
                if lv >= thr:
                    out.append(tuple(r))
            return out

        for m in (parm3, parm4, parm5):
            out = pick(m)
            if out:
                return out
        return []
    except Exception as e:
        print(f"[warn] {p}: {e}")
        return []

def enrich_and_write(subdir: Path, base_rows, contig_len_map, genes_path: Path):
    # Load query_filtered.csv (per chromosome)
    genes = {}
    if genes_path.is_file():
        with genes_path.open("r", newline="") as f:
            r = csv.DictReader(f)
            for rec in r:
                fid = rec.get("full_id") or rec.get("fullid") or rec.get("id")
                if not fid:
                    continue
                try:
                    start = float(rec.get("start", "nan"))
                except ValueError:
                    start = float("nan")
                try:
                    stop = float(rec.get("stop", "nan"))
                except ValueError:
                    stop = float("nan")
                scaff = rec.get("scaff", "")
                genes[fid] = (start, stop, scaff)
    else:
        print(f"[warn] {genes_path} not found; enrichment will be partial.")

    # Deduplicate like R's !duplicated
    seen = set()
    dedup = []
    for r in base_rows:
        if r not in seen:
            seen.add(r)
            dedup.append(list(r))

    out_rows = []
    for base in dedup:
        try:
            subj_start = float(base[8])
            subj_end   = float(base[9])
            orientation = "plus" if subj_start < subj_end else "minus"
        except ValueError:
            orientation = "minus"

        qname = base[0]
        start, stop, scaff = genes.get(qname, (float("nan"), float("nan"), ""))

        subj = base[1]
        clen = contig_len_map.get(subj, float("nan"))

        # Insert genomic_start/genomic_end BEFORE E_value (index 10)
        base.insert(IDX_EVALUE, str(start))   # now at index 10
        base.insert(IDX_EVALUE + 1, str(stop))# now at index 11

        # Append the rest (after gene_length)
        base += [str(scaff), str(start), str(stop), orientation, str(clen)]

        row_map = {
            "query_name": base[0],
            "subject_name": base[1],
            "identity": base[2],
            "length": base[3],
            "mismatch": base[4],
            "gap": base[5],
            "query_start": base[6],
            "query_end": base[7],
            "subject_start": base[8],
            "subject_end": base[9],
            "genomic_start": base[10],
            "genomic_end": base[11],
            "E_value": base[12],
            "bit_score": base[13],
            "gene_length": base[14],
            "ref_scaff": base[15],
            "scaff_start": base[16],
            "scaff_stop": base[17],
            "orientation": base[18],
            "contig_length": base[19],
        }
        out_rows.append([row_map[k] for k in FINAL_HEADER])

    # Arrange by genomic_start asc
    gs_idx = FINAL_HEADER.index("genomic_start")
    def srt(row):
        try:
            return float(row[gs_idx])
        except ValueError:
            return math.inf
    out_rows.sort(key=srt)

    out_path = subdir / "contigs_of_interest.csv"
    with out_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(FINAL_HEADER)
        w.writerows(out_rows)
    print(f"[ok] {subdir.name}: wrote {len(out_rows)} rows -> {out_path}")

def main():
    ap = argparse.ArgumentParser(description="Parallel dScaff filter (per-chromosome outputs).")
    ap.add_argument("--mainDir", default=".", help="Root directory (contains 2L,2R,3L,3R,4,X, etc.)")
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1, help="Processes to use")
    ap.add_argument("--params", default="parameters_dScaff.txt", help="Parameters file (in mainDir)")
    ap.add_argument("--contig_lengths", default="lungimi_contiguri.csv", help="Contig lengths CSV (in mainDir)")
    args = ap.parse_args()

    main_dir = Path(args.mainDir).resolve()
    parm1, parm2, parm3, parm4, parm5, parm6, parm7, parm8 = load_parameters(main_dir / args.params)
    contig_len_map = load_contig_lengths(main_dir / args.contig_lengths)

    # Find chromosome subdirs and their queries/*.csv
    subdirs = discover_subdirs(main_dir)
    per_dir_files = {sd: list_query_csvs(sd) for sd in subdirs}
    # Keep only dirs that actually have files to process
    per_dir_files = {sd: fs for sd, fs in per_dir_files.items() if fs}

    if not per_dir_files:
        print("No queries/*.csv found under any chromosome folder.")
        return

    # Submit every file as a separate future, but **track which dir it belongs to**
    by_dir_rows = {sd: [] for sd in per_dir_files}
    fut_to_dir = {}
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for sd, files in per_dir_files.items():
            for fp in files:
                fut = ex.submit(filter_one_file, str(fp), parm2, parm3, parm4, parm5)
                fut_to_dir[fut] = sd

        # Collect into the correct chromosome bucket
        for fut in as_completed(fut_to_dir):
            sd = fut_to_dir[fut]
            rows = fut.result()
            if rows:
                by_dir_rows[sd].extend(rows)

    # Now write one output per chromosome directory
    for sd, rows in by_dir_rows.items():
        genes_path = sd / "query_filtered.csv"
        enrich_and_write(sd, rows, contig_len_map, genes_path)

if __name__ == "__main__":
    main()

