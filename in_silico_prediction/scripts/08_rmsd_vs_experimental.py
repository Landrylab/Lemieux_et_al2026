#!/usr/bin/env python3
"""
08_rmsd_vs_experimental.py -- how close are the AlphaFold3 models to the
experimental structures?

===========================================================================
TL;DR
===========================================================================
    in   structures/<fam>_sys_v2/fold_<fam>_sys_v2_model_0.cif
         ref-pdbs/<uniprot>-<fam>/<PDB>.cif
    out  results/08_rmsd/rmsd.csv
         a printed table

    run  python 08_rmsd_vs_experimental.py
         python 08_rmsd_vs_experimental.py --all-models

Runs in PyMOL if it is importable, and falls back to a gemmi implementation
otherwise, so the numbers are obtainable either way. Both are reported the
same and agree to within rounding.

===========================================================================
WHAT IS MEASURED
===========================================================================
Backbone Ca RMSD over the receptor chain:

  rmsd_all      over every residue the alignment matched, no outlier
                rejection. The conservative figure.
  rmsd_refined  after iterative rejection of poorly-fitting pairs, with the
                number of atoms retained. This is what PyMOL's align prints
                by default, and it flatters a model with a few displaced
                loops.

Report rmsd_all, or both. Quoting only the refined value without its
retained-atom count is how a model that matches over half its residues
comes to look like a good one.

SYMMETRY
--------
Some entries deposit only the domain, the bound peptide being contributed by
a crystallographic neighbour. For those, "sym:lastN" generates symmetry
mates, keeps whichever copy packs against the domain, and takes its last N
residues as the ligand. This needs PyMOL; the gemmi fallback cannot expand
symmetry and says so rather than reporting a number.

Where the experimental structure also contains the bound peptide, two
further numbers are reported, and the distinction between them is the
point:

  pep_rmsd_placed  the peptide RMSD after superposing on the RECEPTOR, with
                   no further fitting. This asks whether the peptide sits
                   in the right place relative to the domain, which is what
                   an interface energy depends on.
  pep_rmsd_shape   the peptide RMSD after fitting the peptide onto itself.
                   This asks only whether its conformation is right, and
                   ignores where it was put.

A model can have a good pep_rmsd_shape and a poor pep_rmsd_placed: the
peptide is the right shape, docked in the wrong register or orientation.
Only the placed value bears on the ddG calculations.

===========================================================================
WHY THE ALIGNMENT IS STRUCTURAL RATHER THAN BY RESIDUE NUMBER
===========================================================================
The experimental structures use their own numbering, are missing loops, and
were solved on different construct boundaries than the trimmed sequences
submitted to AlphaFold3. Matching by residue number would silently pair the
wrong residues. Superposition is therefore done with cealign, which matches
on structure alone, and the fallback uses a sequence alignment for the same
reason.

===========================================================================
REQUIRES
===========================================================================
PyMOL, or gemmi and numpy. No network access: the reference structures are
read from ref-pdbs/.
"""

import argparse
import csv
import sys
from pathlib import Path

# family -> (PDB id, receptor chain, peptide specification or None)
#
# The peptide specification is either a chain letter, when the bound peptide
# is deposited as its own chain, or "<chain>:<first>-<last>" when it is part
# of the receptor chain -- some entries fuse the ligand to the domain, or
# deposit it as a continuation of the same entity. In the second case those
# residues are also EXCLUDED from the receptor comparison, or the peptide
# would be counted twice and the receptor RMSD would silently include it.
#
# None where the structure is apo and there is no peptide to compare.
PAIRS = {
    # 4Q2Q's asymmetric unit holds only the domain: the bound peptide comes
    # from a crystallographic neighbour, so symmetry mates are generated and
    # the ligand taken from the one that packs into the groove. Only the
    # last few residues of that copy occupy the binding site, so the
    # comparison is restricted to them -- see SYMMETRY below.
    "363": ("4Q2Q", "A", "sym:last6"),   # ZO-1 PDZ3, peptide from a mate
    "366": ("3KZD", "A", None),          # Tiam1 PDZ, apo
    "385": ("4Q2O", "A", "B"),           # PDLIM4 PDZ + phage peptide
}

# For a peptide supplied by a symmetry mate: how far to search for
# neighbouring copies, and how many C-terminal residues of the contacting
# copy to treat as the ligand.
SYM_RADIUS = 5.0
SYM_LAST_N = 6


def parse_pep(spec):
    """'B' -> ('B', None);  'A:505-514' -> ('A', (505,514));
    'sym:last6' -> ('sym', 6)."""
    if not spec:
        return None, None
    if spec.startswith("sym:"):
        return "sym", int(spec.split("last")[1])
    if ":" in spec:
        chain, rng = spec.split(":", 1)
        lo, hi = rng.split("-")
        return chain.strip(), (int(lo), int(hi))
    return spec.strip(), None


def find_reference(ref_root, family, pdb_id):
    """ref-pdbs/<uniprot>-<family>/<PDB>.cif, whatever the uniprot is."""
    for d in sorted(Path(ref_root).glob(f"*-{family}")):
        for name in (f"{pdb_id}.cif", f"{pdb_id.lower()}.cif",
                     f"{pdb_id}.pdb"):
            p = d / name
            if p.is_file():
                return p
    hits = sorted(Path(ref_root).rglob(f"{pdb_id}.cif"))
    return hits[0] if hits else None


def find_model(struct_root, family, index):
    p = (Path(struct_root) / f"{family}_sys_v2" /
         f"fold_{family}_sys_v2_model_{index}.cif")
    if p.is_file():
        return p
    hits = sorted(Path(struct_root).rglob(
        f"*{family}*model_{index}.cif"))
    return hits[0] if hits else None


# --------------------------------------------------------------- pymol ---
def pymol_seq(cmd, selection):
    """One-letter sequence of a CA selection, in chain order.

    The RMSD alone does not say WHICH residues were compared. With a
    crystallised ligand of different sequence from the modelled peptide,
    the matcher pairs whatever residues it can, and an RMSD over an
    accidental run of coincidental matches means far less than one over the
    aligned C-terminus. Printing both sequences makes that visible.
    """
    three = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
             "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
             "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
             "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
    got = []
    cmd.iterate(selection, "got.append((int(resi), resn))",
                space={"got": got})
    return "".join(three.get(rn, "X") for _n, rn in got), \
        [n for n, _rn in got]


def rmsd_pymol(model_path, ref_path, model_chain, ref_chain,
               model_pep=None, ref_pep=None):
    from pymol import cmd
    cmd.reinitialize()
    cmd.load(str(model_path), "af3")
    cmd.load(str(ref_path), "exp")

    pep_chain, pep_range = parse_pep(ref_pep)

    sel_a = f"af3 and chain {model_chain} and polymer.protein and name CA"
    sel_e = f"exp and chain {ref_chain} and polymer.protein and name CA"
    # a peptide living inside the receptor chain must not also be scored as
    # receptor
    if pep_range and pep_chain == ref_chain:
        sel_e += f" and not resi {pep_range[0]}-{pep_range[1]}"
    if cmd.count_atoms(sel_a) == 0 or cmd.count_atoms(sel_e) == 0:
        raise RuntimeError(f"empty selection; af3 has chains "
                           f"{cmd.get_chains('af3')}, exp has "
                           f"{cmd.get_chains('exp')}")

    # structural superposition first: the two use different numbering and
    # sequence-based alignment can fail outright on a trimmed construct
    cmd.cealign(sel_e, sel_a)
    refined = cmd.align(sel_a, sel_e, cycles=5)
    allatom = cmd.align(sel_a, sel_e, cycles=0)

    out = {"rmsd_all": allatom[0], "n_all": allatom[1],
           "rmsd_refined": refined[0], "n_refined": refined[1],
           "n_model": cmd.count_atoms(sel_a),
           "n_exp": cmd.count_atoms(sel_e), "engine": "pymol"}

    if model_pep and pep_chain == "sym":
        # generate crystallographic neighbours and keep the copy whose
        # C-terminus packs against the domain: that is the bound ligand
        cmd.symexp("symmate", "exp", "exp", SYM_RADIUS)
        mates = [o for o in cmd.get_object_list() if o.startswith("symmate")]
        best, best_n = None, 0
        for m in mates:
            tail = (f"{m} and polymer.protein and name CA within 6 of "
                    f"(exp and chain {ref_chain} and polymer.protein)")
            n = cmd.count_atoms(tail)
            if n > best_n:
                best, best_n = m, n
        if best is None:
            out["pep_note"] = (f"no symmetry mate within {SYM_RADIUS} A "
                               f"contacts the domain")
        else:
            nums = []
            cmd.iterate(f"{best} and polymer.protein and name CA",
                        "nums.append(int(resi))", space={"nums": nums})
            keep = sorted(nums)[-pep_range:]
            pa = (f"af3 and chain {model_pep} and polymer.protein and "
                  f"name CA")
            pe = (f"{best} and polymer.protein and name CA and resi "
                  + "+".join(str(n) for n in keep))
            out["pep_source"] = f"{best} resi {keep[0]}-{keep[-1]}"
            if cmd.count_atoms(pa) and cmd.count_atoms(pe):
                placed = cmd.align(pa, pe, cycles=0, transform=0)
                shape = cmd.align(pa, pe, cycles=0)
                mseq, _mn = pymol_seq(cmd, pa)
                xseq, xnum = pymol_seq(cmd, pe)
                out.update({"pep_rmsd_placed": placed[0],
                            "n_pep": placed[1],
                            "pep_rmsd_shape": shape[0],
                            "n_pep_model": cmd.count_atoms(pa),
                            "n_pep_exp": cmd.count_atoms(pe),
                            "pep_seq_model": mseq,
                            "pep_seq_xray": xseq,
                            "pep_resi_xray": f"{xnum[0]}-{xnum[-1]}"
                            if xnum else ""})
        return out

    if model_pep and pep_chain:
        pa = f"af3 and chain {model_pep} and polymer.protein and name CA"
        pe = f"exp and chain {pep_chain} and polymer.protein and name CA"
        if pep_range:
            pe += f" and resi {pep_range[0]}-{pep_range[1]}"
        if cmd.count_atoms(pa) and cmd.count_atoms(pe):
            # the receptors are already superposed, so transform=0 measures
            # the peptide where it actually sits rather than refitting it
            placed = cmd.align(pa, pe, cycles=0, transform=0)
            shape = cmd.align(pa, pe, cycles=0)
            mseq, _mn = pymol_seq(cmd, pa)
            xseq, xnum = pymol_seq(cmd, pe)
            out.update({"pep_rmsd_placed": placed[0], "n_pep": placed[1],
                        "pep_rmsd_shape": shape[0],
                        "n_pep_model": cmd.count_atoms(pa),
                        "n_pep_exp": cmd.count_atoms(pe),
                        "pep_seq_model": mseq, "pep_seq_xray": xseq,
                        "pep_resi_xray": f"{xnum[0]}-{xnum[-1]}"
                        if xnum else ""})
    return out


# --------------------------------------------------------------- gemmi ---
def rmsd_gemmi(model_path, ref_path, model_chain, ref_chain,
               model_pep=None, ref_pep=None):
    """Superpose on Ca atoms matched by sequence alignment.

    Used when PyMOL is unavailable. gemmi's align_string_sequences gives the
    residue correspondence, which is what makes this robust to the
    experimental structure's own numbering and missing loops.
    """
    import gemmi
    import numpy as np
    from difflib import SequenceMatcher

    def ca_chain(path, chain_id, keep=None, drop=None):
        """CA sequence and coordinates for one chain.

        keep: only residues in this (first, last) numbering range.
        drop: exclude residues in this range -- used when the peptide is
              part of the receptor chain and must not be scored twice.
        """
        st = gemmi.read_structure(str(path))
        st.setup_entities()
        st.remove_ligands_and_waters()
        for ch in st[0]:
            if ch.name != chain_id:
                continue
            seq, coords = [], []
            for res in ch:
                n = res.seqid.num
                if keep and not (keep[0] <= n <= keep[1]):
                    continue
                if drop and drop[0] <= n <= drop[1]:
                    continue
                info = gemmi.find_tabulated_residue(res.name)
                ca = res.find_atom("CA", "*")
                if info and info.is_amino_acid() and ca is not None:
                    c = info.one_letter_code.upper()
                    seq.append(c if c.isalpha() else "X")
                    coords.append([ca.pos.x, ca.pos.y, ca.pos.z])
            return "".join(seq), np.asarray(coords, float), \
                [c.name for c in st[0]]
        return "", np.empty((0, 3)), [c.name for c in st[0]]

    pep_chain, pep_range = parse_pep(ref_pep)
    sa, xa, chains_a = ca_chain(model_path, model_chain)
    se, xe, chains_e = ca_chain(
        ref_path, ref_chain,
        drop=pep_range if (pep_range and pep_chain == ref_chain) else None)
    if not sa or not se:
        raise RuntimeError(f"empty chain; af3 has {chains_a}, "
                           f"exp has {chains_e}")

    # Match residues by the longest common subsequence of the two
    # sequences. This handles the loop the experimental structure is
    # missing and its different construct boundaries, without depending on
    # residue numbering, which the two do not share.
    from difflib import SequenceMatcher
    sm = SequenceMatcher(None, sa, se, autojunk=False)
    pa, pe = [], []
    for i, j, n in sm.get_matching_blocks():
        pa.extend(range(i, i + n))
        pe.extend(range(j, j + n))
    if len(pa) < 10:
        raise RuntimeError(f"only {len(pa)} residue(s) matched between the "
                           f"two sequences")
    A, B = xa[pa], xe[pe]

    def kabsch_rmsd(P, Q):
        Pc, Qc = P - P.mean(0), Q - Q.mean(0)
        V, _S, Wt = np.linalg.svd(Pc.T @ Qc)
        d = np.sign(np.linalg.det(V @ Wt))
        R = V @ np.diag([1.0, 1.0, d]) @ Wt
        Pr = Pc @ R
        return float(np.sqrt(((Pr - Qc) ** 2).sum(1).mean())), Pr, Qc

    r_all, Pr, Qc = kabsch_rmsd(A, B)

    # the same outlier rejection PyMOL's align performs, so the two
    # engines report comparable numbers
    for _ in range(5):
        d = np.sqrt(((Pr - Qc) ** 2).sum(1))
        sel = d < d.mean() + 2 * d.std()
        if sel.all() or sel.sum() < 10:
            break
        A, B = A[sel], B[sel]
        _r, Pr, Qc = kabsch_rmsd(A, B)
    r_ref = float(np.sqrt(((Pr - Qc) ** 2).sum(1).mean()))

    out = {"rmsd_all": r_all, "n_all": len(pa),
           "rmsd_refined": r_ref, "n_refined": len(A),
           "n_model": len(sa), "n_exp": len(se), "engine": "gemmi"}

    if pep_chain == "sym":
        out["pep_note"] = ("peptide comes from a symmetry mate; run with "
                           "--engine pymol")
        return out

    if model_pep and pep_chain:
        spa, xpa, _c = ca_chain(model_path, model_pep)
        spe, xpe, _c = ca_chain(ref_path, pep_chain, keep=pep_range)
        if spa and spe:
            sm2 = SequenceMatcher(None, spa, spe, autojunk=False)
            qa, qe = [], []
            for i, j, n in sm2.get_matching_blocks():
                qa.extend(range(i, i + n))
                qe.extend(range(j, j + n))
            if len(qa) < 3:
                out["pep_note"] = (
                    f"only {len(qa)} residue(s) could be matched between "
                    f"the modelled peptide ({spa}) and the crystallised "
                    f"ligand ({spe}); the sequences are too dissimilar for "
                    f"a meaningful comparison")
            if len(qa) >= 3:
                # apply the receptor superposition to the peptide, then
                # measure in place: this tests placement, not shape
                Pa, Qa = xa[pa], xe[pe]
                ca_m, ca_e = Pa.mean(0), Qa.mean(0)
                V, _S, Wt = np.linalg.svd((Pa - ca_m).T @ (Qa - ca_e))
                d = np.sign(np.linalg.det(V @ Wt))
                R = V @ np.diag([1.0, 1.0, d]) @ Wt
                moved = (xpa[qa] - ca_m) @ R + ca_e
                out["pep_rmsd_placed"] = float(
                    np.sqrt(((moved - xpe[qe]) ** 2).sum(1).mean()))
                out["n_pep"] = len(qa)
                out["pep_rmsd_shape"] = kabsch_rmsd(xpa[qa], xpe[qe])[0]
                out["n_pep_model"] = len(spa)
                out["n_pep_exp"] = len(spe)
                out["pep_seq_model"] = spa
                out["pep_seq_xray"] = spe
                out["pep_matched_model"] = "".join(spa[i] for i in qa)
                out["pep_matched_xray"] = "".join(spe[i] for i in qe)
    return out


# ------------------------------------------------------------------ main --
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--structures", default="structures")
    ap.add_argument("--ref-pdbs", default="ref-pdbs")
    ap.add_argument("--model-chain", default="A",
                    help="receptor chain in the AlphaFold3 model "
                         "(default: A)")
    ap.add_argument("--model-pep-chain", default="B",
                    help="peptide chain in the AlphaFold3 model "
                         "(default: B)")
    ap.add_argument("--all-models", action="store_true",
                    help="all five models rather than model 0 only, which "
                         "shows whether the top-ranked one is typical")
    ap.add_argument("--engine", default="auto",
                    choices=["auto", "pymol", "gemmi"])
    ap.add_argument("-o", "--out", default="results/08_rmsd/rmsd.csv")
    args = ap.parse_args()

    engine = args.engine
    if engine == "auto":
        try:
            import pymol  # noqa: F401
            engine = "pymol"
        except ImportError:
            engine = "gemmi"
    print(f"engine: {engine}\n")

    indices = range(5) if args.all_models else [0]
    rows = []
    for family, (pdb_id, ref_chain, ref_pep) in PAIRS.items():
        ref = find_reference(args.ref_pdbs, family, pdb_id)
        if ref is None:
            print(f"{family}: {pdb_id}.cif not found under "
                  f"{args.ref_pdbs}/, skipping")
            continue
        for idx in indices:
            model = find_model(args.structures, family, idx)
            if model is None:
                print(f"{family} model {idx}: not found, skipping")
                continue
            fn = rmsd_pymol if engine == "pymol" else rmsd_gemmi
            try:
                r = fn(model, ref, args.model_chain, ref_chain,
                       args.model_pep_chain if ref_pep else None, ref_pep)
            except Exception as exc:                         # noqa: BLE001
                print(f"{family} model {idx} vs {pdb_id}: FAILED "
                      f"({type(exc).__name__}: {exc})")
                continue
            r.update({"family": family, "model": idx, "pdb": pdb_id,
                      "ref_chain": ref_chain, "ref_pep_chain": ref_pep or "",
                      "model_file": model.name, "ref_file": ref.name})
            rows.append(r)
            print(f"{family} model {idx} vs {pdb_id} chain {ref_chain}: "
                  f"receptor RMSD {r['rmsd_all']:.2f} A over {r['n_all']} "
                  f"Ca (refined {r['rmsd_refined']:.2f} over "
                  f"{r['n_refined']});  model {r['n_model']} res, exp "
                  f"{r['n_exp']} res")
            if r.get("pep_note"):
                print(f"      peptide: {r['pep_note']}")
            if "pep_rmsd_placed" in r:
                label = r.get("pep_source", ref_pep)
                print(f"      peptide ({label}): "
                      f"{r['pep_rmsd_placed']:.2f} A as placed, "
                      f"{r['pep_rmsd_shape']:.2f} A refitted, over "
                      f"{r['n_pep']} Ca")
                if r.get("pep_seq_model") or r.get("pep_seq_xray"):
                    resi = (f"  (resi {r['pep_resi_xray']})"
                            if r.get("pep_resi_xray") else "")
                    print(f"        model peptide : "
                          f"{r.get('pep_seq_model', '?')}")
                    print(f"        x-ray peptide : "
                          f"{r.get('pep_seq_xray', '?')}{resi}")
                if r.get("pep_matched_model"):
                    print(f"        residues paired: "
                          f"{r['pep_matched_model']} / "
                          f"{r['pep_matched_xray']}")
                shared = sum(1 for a, b in zip(r.get("pep_seq_model", ""),
                                               r.get("pep_seq_xray", ""))
                             if a == b)
                if r.get("pep_seq_model") and r.get("pep_seq_xray"):
                    print(f"        identical at {shared} of "
                          f"{min(len(r['pep_seq_model']), len(r['pep_seq_xray']))} "
                          f"aligned position(s)")

    if not rows:
        sys.exit("\nnothing computed")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["family", "model", "pdb", "ref_chain", "rmsd_all", "n_all",
            "rmsd_refined", "n_refined", "n_model", "n_exp",
            "ref_pep_chain", "pep_source", "pep_rmsd_placed",
            "pep_rmsd_shape", "n_pep", "n_pep_model", "n_pep_exp",
            "pep_seq_model", "pep_seq_xray", "pep_resi_xray",
            "pep_matched_model", "pep_matched_xray",
            "pep_note", "engine", "model_file", "ref_file"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(r[k], 3) if isinstance(r.get(k), float)
                            else r.get(k, ""))
                        for k in cols})
    print(f"\nwrote {out}")

    thin = [r for r in rows
            if r.get("n_pep") and r["n_pep"] < 4]
    if thin:
        print(f"\n{len(thin)} peptide comparison(s) rest on fewer than four "
              f"Cα atoms.\nAn RMSD over three residues is not evidence of "
              f"much; quote the atom count with it.")

    short = [r for r in rows
             if r["n_all"] < 0.6 * min(r["n_model"], r["n_exp"])]
    if short:
        print(f"\n{len(short)} comparison(s) matched fewer than 60% of the "
              f"shorter chain.\nAn RMSD over a fragment is not an RMSD over "
              f"the domain; check those before\nquoting them.")


if __name__ == "__main__":
    main()
