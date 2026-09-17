#!/usr/bin/env python3
"""
01_extract_libraries.py -- turn the raw experiment table into one clean
peptide library per family.

===========================================================================
TL;DR
===========================================================================
    in   PCA_result.csv          (EXPERIMENT_CSV in paths.env)
         pep-pwmax-<family>.fa      (REFERENCE_FASTA, one per family)
    out  results/01_libraries/fam<family>.csv   one library per family
         results/01_libraries/summary.csv        counts and ranges

    run  python 01_extract_libraries.py
         python 01_extract_libraries.py --families 363
         python 01_extract_libraries.py --dry-run

No external tools. Pure pandas filtering plus a sequence diff against the
family reference peptide.

===========================================================================
WHAT IT DOES
===========================================================================
For every family listed in FAMILIES:

  1. select the rows of the experiment table belonging to that family
  2. drop anything unusable and say why:
       - sequences containing a stop codon (*)
       - sequences with non-standard residue letters (X, B, Z, gaps)
       - sequences whose length differs from the family reference
       - the reference peptide itself (zero substitutions: ddG is 0 by
         definition, so it is not something to predict)
  3. diff each surviving sequence against the reference and record the
     substitutions and how many there are
  4. average the measured value across repeated measurements of the same
     peptide, and keep the spread -- that spread is the assay noise floor
     and is the honest yardstick for judging any prediction later
  5. assign a stable id of the form <family>_<0001..>, numbered over the
     CLEAN sequences in file order

Nothing is silently dropped: every filter reports its count, and
summary.csv records the totals per family.

===========================================================================
WHY THE ID IS NOT THE JOIN KEY
===========================================================================
Ids are positional. If the input table ever gains or loses a row, every id
after that point shifts, and two files generated from different versions of
the table will join to different peptides without raising an error. The id
is written for convenience and for readable filenames; downstream scripts
join on `sequence`. Both columns are present in every output file.

===========================================================================
OUTPUT COLUMNS
===========================================================================
    family        family id, as given in FAMILIES
    id            <family>_0001 ... , positional, see caveat above
    sequence      the peptide, uppercase, one letter per residue
    reference     the family reference peptide, repeated on every row
    n_mut         number of substitutions from the reference
    mutations     e.g. S1A,G2D -- wild type, position (1-based), mutant
    positions     e.g. 1,2 -- just the positions, for grouping later
    ppi_exp       measured interaction score, averaged over replicates
    ppi_sd        spread across replicates, blank when measured once
    n_rep         how many rows of the input table gave this sequence

summary.csv holds one row per family: counts kept and dropped by reason,
the reference sequence, and the range of the measured value.

===========================================================================
REQUIRES
===========================================================================
pandas.  paths.env must define EXPERIMENT_CSV, REFERENCE_FASTA, FAMILIES,
SEQ_COL, FAMILY_COL, PPI_COL and RESULTS_DIR.
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: pandas not found.  pip install pandas")

STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


# --------------------------------------------------------------------------
def load_env(path="paths.env"):
    """Read KEY=value lines. The file is also valid bash, so the same
    settings can be picked up with `source paths.env` in a shell."""
    cfg = {}
    p = Path(path)
    if not p.is_file():
        sys.exit(f"ERROR: {p} not found.\n"
                 f"  cp paths.template.env paths.env   and edit it,\n"
                 f"  or pass --paths /path/to/paths.env")
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        if line.startswith("export "):
            line = line[7:].strip()
        k, _, v = line.partition("=")
        cfg[k.strip()] = v.strip().strip('"').strip("'")
    return cfg


def read_fasta(path):
    """First sequence in a FASTA, or the whole file if it has no header."""
    text = Path(path).read_text().strip()
    if not text:
        return ""
    if not text.startswith(">"):
        return "".join(text.split()).upper()
    out = []
    for line in text.splitlines()[1:]:
        if line.startswith(">"):
            break
        out.append(line.strip())
    return "".join(out).upper()


def to_number(series):
    """Numbers, tolerating a decimal comma and surrounding quotes.

    The experiment table stores values like "8,4284e-01" -- a European
    decimal comma inside a quoted field. Left to pandas these become
    strings or NaN, so the conversion is explicit.
    """
    return pd.to_numeric(
        series.astype(str).str.strip().str.strip('"')
              .str.replace(",", ".", regex=False)
              .replace({"": None, "nan": None, "NA": None, "None": None}),
        errors="coerce")


def diff_to_reference(seq, ref):
    """-> (n_mut, 'S1A,G2D', '1,2'). Assumes equal lengths."""
    diffs = [(ref[i], i + 1, seq[i])
             for i in range(len(ref)) if seq[i] != ref[i]]
    return (len(diffs),
            ",".join(f"{w}{p}{m}" for w, p, m in diffs),
            ",".join(str(p) for _w, p, _m in diffs))


# --------------------------------------------------------------------------
def process_family(family, raw, ref, cfg):
    """-> (DataFrame of the clean library, dict of counts)."""
    seq_col, ppi_col = cfg["SEQ_COL"], cfg["PPI_COL"]
    counts = Counter()

    df = raw.copy()
    counts["input rows"] = len(df)

    df[seq_col] = df[seq_col].astype(str).str.strip().str.upper()

    keep = df[seq_col].str.len() > 0
    counts["blank sequence"] = int((~keep).sum())
    df = df[keep]

    has_stop = df[seq_col].str.contains(r"\*", regex=True)
    counts["stop codon (*)"] = int(has_stop.sum())
    df = df[~has_stop]

    clean = df[seq_col].apply(lambda s: set(s) <= STANDARD)
    counts["non-standard residues"] = int((~clean).sum())
    df = df[clean]

    right_len = df[seq_col].str.len() == len(ref)
    counts[f"length != {len(ref)}"] = int((~right_len).sum())
    df = df[right_len]

    is_ref = df[seq_col] == ref
    counts["identical to reference"] = int(is_ref.sum())
    df = df[~is_ref]

    if df.empty:
        return pd.DataFrame(), counts

    df[ppi_col] = to_number(df[ppi_col])
    missing = df[ppi_col].isna()
    counts["no measured value"] = int(missing.sum())
    df = df[~missing]

    # collapse repeated measurements of the same peptide; the spread across
    # them is the assay noise floor
    grouped = (df.groupby(seq_col, sort=False)[ppi_col]
                 .agg(ppi_exp="mean", ppi_sd="std", n_rep="count")
                 .reset_index())
    counts["unique sequences"] = len(grouped)
    counts["repeated measurements"] = int((grouped["n_rep"] > 1).sum())

    rows = []
    for i, r in enumerate(grouped.itertuples(index=False), 1):
        seq = getattr(r, seq_col) if hasattr(r, seq_col) else r[0]
        n_mut, muts, pos = diff_to_reference(seq, ref)
        rows.append({
            "family": family,
            "id": f"{family}_{i:04d}",
            "sequence": seq,
            "reference": ref,
            "n_mut": n_mut,
            "mutations": muts,
            "positions": pos,
            "ppi_exp": round(float(r.ppi_exp), 6),
            "ppi_sd": ("" if pd.isna(r.ppi_sd)
                       else round(float(r.ppi_sd), 6)),
            "n_rep": int(r.n_rep),
        })

    out = pd.DataFrame(rows)
    for k, v in Counter(out["n_mut"]).items():
        counts[f"{k}-point"] = v
    return out, counts


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paths", default="paths.env",
                    help="settings file (default: paths.env)")
    ap.add_argument("--families", default=None,
                    help="comma-separated subset; default is FAMILIES "
                         "from paths.env")
    ap.add_argument("--experiment-csv", default=None,
                    help="override EXPERIMENT_CSV")
    ap.add_argument("--out-dir", default=None,
                    help="where to write (default: "
                         "RESULTS_DIR/01_libraries from paths.env)")
    ap.add_argument("--dry-run", action="store_true",
                    help="report the counts, write nothing")
    args = ap.parse_args()

    cfg = load_env(args.paths)
    csv_path = Path(args.experiment_csv or cfg["EXPERIMENT_CSV"])
    if not csv_path.is_file():
        sys.exit(f"ERROR: experiment table not found: {csv_path}")

    families = [f.strip() for f in
                (args.families or cfg["FAMILIES"]).split(",") if f.strip()]
    out_dir = Path(args.out_dir) if args.out_dir else \
        Path(cfg.get("RESULTS_DIR", "results")) / "01_libraries"

    seq_col, fam_col = cfg["SEQ_COL"], cfg["FAMILY_COL"]
    ppi_col = cfg["PPI_COL"]

    table = pd.read_csv(csv_path, dtype=str, skipinitialspace=True)
    table.columns = [c.strip() for c in table.columns]
    for col in (seq_col, fam_col, ppi_col):
        if col not in table.columns:
            sys.exit(f"ERROR: column '{col}' not in {csv_path}\n"
                     f"  columns present: {list(table.columns)}")
    table[fam_col] = table[fam_col].astype(str).str.strip()

    print(f"experiment table : {csv_path}  ({len(table)} rows)")
    print(f"families in file : "
          f"{', '.join(sorted(set(table[fam_col].dropna()))[:10])}")
    print(f"processing       : {', '.join(families)}")
    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"writing to       : {out_dir.resolve()}")

    summary = []
    for family in families:
        ref_path = Path(cfg["REFERENCE_FASTA"].format(family))
        if not ref_path.is_file():
            print(f"\n{family}: reference FASTA missing ({ref_path}), "
                  f"skipping")
            continue
        ref = read_fasta(ref_path)
        raw = table[table[fam_col] == family]
        if raw.empty:
            print(f"\n{family}: no rows in the experiment table, skipping")
            continue

        lib, counts = process_family(family, raw, ref, cfg)

        print(f"\n\033[1mfamily {family}\033[0m")
        print(f"  reference {ref}  ({len(ref)} aa, {ref_path.name})")
        print(f"  {counts['input rows']} row(s) in the table")
        for reason in ("blank sequence", "stop codon (*)",
                       "non-standard residues", f"length != {len(ref)}",
                       "identical to reference", "no measured value"):
            if counts.get(reason):
                print(f"    dropped {counts[reason]:>5}  {reason}")
        if lib.empty:
            print("  nothing left after filtering")
            continue

        by_n = {k: v for k, v in sorted(counts.items())
                if k.endswith("-point")}
        print(f"  kept {len(lib)} peptide(s): "
              + ", ".join(f"{k} {v}" for k, v in by_n.items()))
        if counts.get("repeated measurements"):
            print(f"    {counts['repeated measurements']} peptide(s) "
                  f"measured more than once; their spread is in ppi_sd")
        print(f"  {ppi_col} range "
              f"{lib['ppi_exp'].min():.3f} .. {lib['ppi_exp'].max():.3f}")

        if not args.dry_run:
            dest = out_dir / f"fam{family}.csv"
            lib.to_csv(dest, index=False)
            print(f"  -> {dest}")

        summary.append({
            "family": family,
            "reference": ref,
            "length": len(ref),
            "input_rows": counts["input rows"],
            "kept": len(lib),
            **{f"n_mut_{k.split('-')[0]}": v for k, v in by_n.items()},
            "dropped_stop": counts.get("stop codon (*)", 0),
            "dropped_nonstandard": counts.get("non-standard residues", 0),
            "dropped_length": counts.get(f"length != {len(ref)}", 0),
            "dropped_reference": counts.get("identical to reference", 0),
            "repeated_measurements": counts.get("repeated measurements", 0),
            "ppi_min": round(float(lib["ppi_exp"].min()), 6),
            "ppi_max": round(float(lib["ppi_exp"].max()), 6),
        })

    if not summary:
        sys.exit("\nnothing written")

    if not args.dry_run:
        s = pd.DataFrame(summary)
        dest = out_dir / "summary.csv"
        s.to_csv(dest, index=False)
        print(f"\nwrote {dest}")
        print("\nNext: python 02_foldx_ddg.py")
    else:
        print("\n(dry run: nothing written)")


if __name__ == "__main__":
    main()
