#!/usr/bin/env python3
"""
03_flexddg_ddg.py -- binding ddG for every peptide variant, by flex ddG.

===========================================================================
TL;DR
===========================================================================
    in   results/01_libraries/fam<N>.csv        (from 01)
         structures/**/fold_<N>*model_<i>.cif   (i = MODEL_INDEX)
    out  results/03_flexddg/fam<N>.csv          one row per variant
         results/03_flexddg/fam<N>_perjob.csv   one row per ensemble member
         results/03_flexddg/work_fam<N>/        the db3 audit trail

    run  python 03_flexddg_ddg.py                  every family
         python 03_flexddg_ddg.py --families 363   one family
         python 03_flexddg_ddg.py --max-mut 1      single mutants only
         python 03_flexddg_ddg.py --dry-run        plan and cost, no compute
         python 03_flexddg_ddg.py --resume         after an interruption

Positive ddG means destabilising, i.e. weaker binding. Same convention as
02_foldx_ddg.py.

===========================================================================
WHAT IT RUNS
===========================================================================
flex ddG (Barlow et al., J Phys Chem B 2018) through PyRosetta, executing
the tutorial's own RosettaScripts protocol, ddG-backrub.xml, unmodified.

One job = one (variant, ensemble member). Inside each job the protocol:

    BackrubMover        FLEXDDG_BACKRUB_TRIALS Monte Carlo backbone moves,
                        giving one ensemble member
    PackRotamersMover   repacks side chains twice on that same member: once
                        with the resfile applied (mutant), once without (wt)
    MinMover            constrained minimisation of both
    InterfaceDdGMover   separates the peptide chain and scores four states
    ReportToDB          writes them to ddG.db3

    ddG = (bound_mut - unbound_mut) - (bound_wt - unbound_wt)

Building wild type and mutant on the SAME backbone is the point: sampling
noise largely cancels in the subtraction, and backbone accommodation is
modelled rather than assumed. That is the difference from FoldX BuildModel
in 02, which repacks on a frozen backbone.

Score function: talaris2014, via -restore_talaris_behavior, which is what
the protocol was parameterised against.

Multi-point variants are applied simultaneously, one PIKAA line per position
in a single resfile, so coupling between positions is sampled rather than
assumed additive.

===========================================================================
THE ONE THING THAT SILENTLY GIVES WRONG NUMBERS
===========================================================================
Run through PyRosetta's XmlObjects there is no JD2 job, so ReportToDB tags
all four output structures "EMPTY_JOB_use_jd2_0000" and the
bound/unbound x wt/mutant assignment must be positional. STATE_ORDER below
holds that order. Getting it backwards flips the sign of every ddG while the
numbers still look reasonable.

Two guards: 00_verify_setup.py check 6 tests the order on a known complex,
and every job here checks that its bound states score below its unbound
states, discarding the job if they do not.

===========================================================================
SAMPLING: TWO PUBLISHED PROTOCOLS
===========================================================================
--protocol selects a preset. Both are published and citable.

  barlow   nstruct 35, 35000 backrub trials, stride 35000
           The protocol as published: Barlow et al., J Phys Chem B 2018,
           122:5389. Reported Pearson 0.63 against experiment on a curated
           1240-mutant benchmark.

  hummer   nstruct 1, 3500 backrub trials, stride 35000        [default]
           A reduced protocol used to generate a 20,829-mutation synthetic
           dataset. From the Methods of Hummer, Schneider, Chinery & Deane,
           Nat Comput Sci 2025, 5:635-647 (doi 10.1038/s43588-025-00823-8),
           section "Synthetic ddG data preparation / Flex ddG":

               "to reduce the runtime of Flex ddG, we set the 'nstruct'
               parameter (number of structures to model) to 1, instead of
               the default value of 35. The 'backrub_trajectory_stride' was
               set to 35,000, and all other parameters were set to default
               values (max_minimization_iter = 5,000,
               abs_score_convergence_thresh = 1.0,
               number_backrub_trials = 3,500)."

           Benchmarked by those authors on 608 single-point antibody-antigen
           mutations from SKEMPI 2.0: Pearson 0.46 with GAM reweighting,
           0.42 without, against 0.20 for FoldX on the same set.

  prelim   nstruct 5, 3500 backrub trials
           Not published. A first pass: cheap enough to run all families
           in an afternoon, and unlike nstruct 1 it yields a standard
           deviation across the ensemble, so a variant whose ddG is
           unstable is visible. Anything reported in a paper should be
           rerun at 'barlow'; use --tag to keep the two side by side.

  custom   whatever FLEXDDG_NSTRUCT and FLEXDDG_BACKRUB_TRIALS say in
           paths.env, or the matching command-line flags.

Which to use is a budget decision. Measured on this hardware, a job costs
about 20 s of fixed work plus roughly 0.0034 s per backrub trial, so the
reduced protocol is about 37 s per job against 144 s, and it needs 35 times
fewer jobs. Against ~32 h per family for barlow, hummer is ~15 min. The
accuracy difference those authors report, 0.46 against Barlow's 0.63, is
the price. --dry-run prints the estimate without running anything.

A note on sign: Hummer et al. define ddG = dG_WT - dG_Mutant, so in their
data a NEGATIVE value is destabilising. This repository uses the opposite
convention throughout, positive = destabilising, matching 02_foldx_ddg.py.
Their numbers are not directly comparable in sign to these.

===========================================================================
OUTPUT
===========================================================================
fam<N>.csv, one row per variant, same columns as 02 so the two concatenate:

    family, id, sequence, n_mut, mutations, ppi_exp, method, structure,
    ddg_binding, ddg_binding_sd, ddg_folding, ddg_folding_sd, n_obs, status

    ddg_binding_sd   spread across the ensemble members
    ddg_folding      blank; InterfaceDdGMover reports binding only
    status           ok | failed -- failed variants keep their row with a
                     blank ddG, so the row count always matches the library

fam<N>_perjob.csv keeps every individual member, for convergence analysis.

===========================================================================
REQUIRES
===========================================================================
PyRosetta      free for academic/non-profit use, not bundled.
               Chaudhury et al., Bioinformatics 2010
ddG-backrub.xml  from Kortemme-Lab/flex_ddG_tutorial (MIT).
               Barlow et al., J Phys Chem B 2018, 122:5389
gemmi          mmCIF to PDB. Wojdyr, JOSS 2022
"""

import argparse
import csv
import os
import re
import shutil
import sqlite3
import sys
import tempfile
import time
from pathlib import Path

# Silence PyRosetta's banner. It prints once per worker process, which for a
# few thousand jobs buries the progress bar. The licence still applies and is
# reproduced in THIRD_PARTY_NOTICES.md.
os.environ.setdefault("PYROSETTA_INIT_SILENT", "1")

try:
    import gemmi
except ImportError:
    sys.exit("ERROR: gemmi not found.  pip install gemmi")

from concurrent.futures import ProcessPoolExecutor, as_completed

# Order in which InterfaceDdGMover writes its four structures.
#
# Determined empirically, not from documentation. Dumping the four totals
# for real jobs gives e.g. [-159.96, -121.86, -159.51, -121.42]: they
# alternate bound, unbound, bound, unbound with a ~38 kcal/mol gap, so
# structures 0 and 2 are the bound states and 1 and 3 the unbound ones.
#
# Which pair is wild type is not visible in those numbers, since swapping
# wt and mut only flips the sign. It is fixed here by the earlier
# observation that the opposite assignment made ddG correlate positively
# with the measured interaction score, i.e. predicted destabilisation
# tracking tighter binding, which is backwards.
#
# Getting this wrong is silent: the magnitudes stay plausible and only the
# sign changes. read_ddg() below rejects any job whose bound states do not
# score below its unbound states, which catches a wrong bound/unbound
# grouping; nothing can catch a wt/mut swap automatically.
STATE_ORDER = [("bound", "wt"), ("unbound", "wt"),
               ("bound", "mut"), ("unbound", "mut")]

# Published parameter sets. See the docstring for citations. Values are
# (nstruct, number_backrub_trials, backrub_trajectory_stride,
#  max_minimization_iter, abs_score_convergence_thresh).
PROTOCOLS = {
    # The protocol as published: Barlow et al., J Phys Chem B 2018.
    "barlow": (35, 35000, 35000, 5000, 1.0),
    # Reduced protocol from Hummer et al., Nat Comput Sci 2025, 5:635.
    "hummer": (1, 3500, 35000, 5000, 1.0),
    # A middle setting for a first pass. nstruct 5 is enough for an
    # ensemble average with a usable standard deviation, which nstruct 1
    # cannot give at all; 3500 trials is the Hummer trial count. This is
    # NOT a published parameter set -- use it to get numbers quickly and
    # to check the pipeline end to end, then rerun at 'barlow' for
    # anything that goes in a paper.
    "prelim": (5, 3500, 3500, 5000, 1.0),
}

OUT_COLUMNS = ["family", "id", "sequence", "n_mut", "mutations", "ppi_exp",
               "method", "structure", "ddg_binding", "ddg_binding_sd",
               "ddg_folding", "ddg_folding_sd", "n_obs", "status"]


# ---------------------------------------------------------------- config --
def load_env(path="paths.env"):
    """Read KEY=value lines. The file is also valid bash."""
    cfg, p = {}, Path(path)
    if not p.is_file():
        sys.exit(f"ERROR: {p} not found. cp paths.template.env paths.env")
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        if line.startswith("export "):
            line = line[7:].strip()
        k, _, v = line.partition("=")
        cfg[k.strip()] = v.strip().strip('"').strip("'")
    return cfg


def hms(seconds):
    if seconds < 90:
        return f"{seconds:.0f} s"
    if seconds < 5400:
        return f"{seconds / 60:.0f} min"
    if seconds < 172800:
        return f"{seconds / 3600:.1f} h"
    return f"{seconds / 86400:.1f} days"


# ------------------------------------------------------------- structures --
def cif_to_pdb(src, dst):
    """Convert to PDB, which is what PyRosetta reads. Returns chain ids."""
    st = gemmi.read_structure(str(src))
    st.setup_entities()
    st.remove_ligands_and_waters()
    st.write_pdb(str(dst))
    return [ch.name for ch in st[0]]


def fill_xml(text, chain, trials, stride=None, max_min=5000,
             conv_thresh=1.0):
    """Substitute the %%script_vars%% the tutorial XML expects.

    XmlObjects does not read these from the options system the way the
    command-line application does, so they are substituted textually and the
    result is checked for anything left over.
    """
    variables = {
        "chainstomove": chain,
        "max_minimization_iter": str(max_min),
        "abs_score_convergence_thresh": str(conv_thresh),
        "number_backrub_trials": str(trials),
        # stride equal to the trial count is the fastest: the final
        # minimisation and packing then run once instead of at every
        # checkpoint. A smaller stride yields intermediate ddG values along
        # the trajectory, which is useful for a convergence curve.
        "backrub_trajectory_stride": str(stride if stride else trials),
        "resfile_relpath": "mutation.resfile",
        "mutate_resfile_relpath": "mutation.resfile",
    }
    for k, v in variables.items():
        text = text.replace(f"%%{k}%%", v)
    missing = sorted(set(re.findall(r"%%(\w+)%%", text)))
    if missing:
        sys.exit(f"ERROR: unsubstituted script_vars in the XML: {missing}")
    return text


# -------------------------------------------------------------- one job ---
def read_ddg(db_path):
    """Binding ddG from one ddG.db3, or None if the states look wrong.

    Per-structure totals come from the weights table in the same file rather
    than being assumed, so the result follows whatever score function the
    run actually used.
    """
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        have = {r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        if "structures" not in have:
            # ReportToDB created the file but the protocol never reached a
            # checkpoint; almost always a stride larger than the trial count
            return None
        rows = sorted(con.execute("SELECT struct_id, tag FROM structures"))
        if len(rows) < 4:
            return None
        names = dict(con.execute(
            "SELECT score_type_id, score_type_name FROM score_types"))
        weights = dict(con.execute(
            "SELECT score_type_id, weight FROM score_function_weights"))
        scores = {}
        for sid, tid, val in con.execute(
                "SELECT struct_id, score_type_id, score_value "
                "FROM structure_scores"):
            scores.setdefault(sid, {})[tid] = val
    finally:
        con.close()

    def total(sid):
        terms = scores.get(sid, {})
        for tid, val in terms.items():
            if names.get(tid, "").lower() == "total_score":
                return val
        return sum(v * weights.get(t, 0.0) for t, v in terms.items())

    state = {STATE_ORDER[i]: total(sid) for i, (sid, _tag)
             in enumerate(rows[:4])}

    # Separating the chains removes favourable contacts, so bound must score
    # below unbound. If it does not, STATE_ORDER is wrong for this build and
    # the ddG would be meaningless.
    bound = (state[("bound", "mut")] + state[("bound", "wt")]) / 2
    unbound = (state[("unbound", "mut")] + state[("unbound", "wt")]) / 2
    if bound >= unbound:
        return None

    return ((state[("bound", "mut")] - state[("unbound", "mut")])
            - (state[("bound", "wt")] - state[("unbound", "wt")]))


def run_job(task):
    """One (variant, ensemble member). Runs in a worker process."""
    vid, muts, member, pdb, xml_text, chain, workdir, init_flags = task
    # absolute: the worker chdirs into wd below, so a relative
    # -out:path:all would be resolved against wd itself and Rosetta
    # would fail to create the database
    wd = (Path(workdir) / vid / f"n{member}").resolve()
    started = time.time()
    cwd = os.getcwd()
    try:
        import pyrosetta
        from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects

        wd.mkdir(parents=True, exist_ok=True)
        shutil.copy(Path(pdb).resolve(), wd / "input.pdb")

        # one PIKAA line per position; several lines means the substitutions
        # are introduced together
        lines = ["NATAA", "start"]
        lines += [f"{pos} {chain} PIKAA {mut}" for _wt, pos, mut in muts]
        (wd / "mutation.resfile").write_text("\n".join(lines) + "\n")

        os.chdir(wd)
        pyrosetta.init(f"{init_flags} -out:path:all {wd}", silent=True)
        pose = pyrosetta.pose_from_pdb("input.pdb")
        XmlObjects.create_from_string(xml_text)\
                  .get_mover("ParsedProtocol").apply(pose)
        os.chdir(cwd)

        db = wd / "ddG.db3"
        if not db.is_file():
            found = list(wd.glob("*.db3"))
            if not found:
                return vid, member, None, "no db3 written", 0.0
            db = found[0]
        ddg = read_ddg(db)
        if ddg is None:
            return vid, member, None, "state check failed", 0.0
        return vid, member, ddg, None, time.time() - started

    except Exception as exc:                                 # noqa: BLE001
        return vid, member, None, f"{type(exc).__name__}: {exc}", 0.0
    finally:
        os.chdir(cwd)


# ------------------------------------------------------------- per family --
def load_library(path):
    """Rows of fam<N>.csv, with the mutations parsed into tuples."""
    out = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            muts = []
            for tok in (row.get("mutations") or "").split(","):
                m = re.match(r"^([A-Z])(\d+)([A-Z])$", tok.strip())
                if m:
                    muts.append((m.group(1), int(m.group(2)), m.group(3)))
            if muts:
                row["_muts"] = muts
                out.append(row)
    return out


def find_structure(cfg, family):
    root = Path(cfg.get("STRUCTURES_DIR", "structures"))
    idx = cfg.get("MODEL_INDEX", "0")
    hits = sorted(root.rglob(f"*{family}*model_{idx}.cif"))
    if not hits:
        sys.exit(f"ERROR: no *{family}*model_{idx}.cif under {root}")
    return hits[0]


def progress(done, total, started, label):
    """Single-line progress bar with an ETA from measured throughput."""
    frac = done / total
    bar = ("#" * int(28 * frac)).ljust(28, ".")
    elapsed = time.time() - started
    eta = elapsed / done * (total - done) if done else 0
    sys.stdout.write(f"\r  [{bar}] {done}/{total} {frac * 100:5.1f}%  "
                     f"{label}  elapsed {hms(elapsed)}  ETA {hms(eta)}   ")
    sys.stdout.flush()


# Cost model fitted to two measurements on the reference workstation:
# 42 s per job at 5000 backrub trials and 144 s at 35000. The fixed term is
# pose loading, packing and minimisation, which do not scale with the number
# of backrub moves; that is why sampling 7x harder costs only 3.4x the time.
# Timings are hardware specific -- rerun bench_flexddg_timing.py and pass
# --per-job-seconds if this machine differs.
JOB_FIXED_SECONDS = 25.0
JOB_SECONDS_PER_TRIAL = 0.0034


def estimate_job_seconds(trials):
    return JOB_FIXED_SECONDS + JOB_SECONDS_PER_TRIAL * trials


def resolve_protocol(args, cfg):
    """-> (nstruct, trials, stride, max_min, conv_thresh).

    A preset supplies all five; an explicit flag overrides its value so a
    single parameter can be varied without leaving the preset behind.
    """
    if args.protocol in PROTOCOLS:
        nstruct, trials, stride, max_min, conv = PROTOCOLS[args.protocol]
    else:
        nstruct = int(cfg.get("FLEXDDG_NSTRUCT", 35))
        trials = int(cfg.get("FLEXDDG_BACKRUB_TRIALS", 35000))
        stride, max_min, conv = trials, 5000, 1.0
    if args.nstruct:
        nstruct = args.nstruct
    if args.backrub_trials:
        trials = args.backrub_trials
        if args.protocol not in PROTOCOLS:
            stride = trials
    if args.stride:
        stride = args.stride

    # A stride larger than the trial count means the trajectory finishes
    # before the first checkpoint, so InterfaceDdGMover never runs and
    # ReportToDB writes a database with no tables in it. The published
    # Hummer settings pair stride 35000 with 3500 trials, which evidently
    # still emitted output on Rosetta 2020.08 but does not on current
    # builds. Clamping to the trial count preserves the intent -- a single
    # checkpoint at the end of the trajectory, which is also the fastest
    # option -- and produces identical sampling.
    if stride > trials:
        stride = trials
    return nstruct, trials, stride, max_min, conv


def process_family(family, cfg, args, xml_raw):
    lib_path = (Path(args.library) if args.library else
                Path(cfg.get("RESULTS_DIR", "results")) / "01_libraries" /
                f"fam{family}.csv")
    if not lib_path.is_file():
        print(f"\n{family}: {lib_path} not found; run 01 first")
        return None

    variants = load_library(lib_path)
    if args.max_mut:
        variants = [v for v in variants if len(v["_muts"]) <= args.max_mut]
    if args.min_mut:
        variants = [v for v in variants if len(v["_muts"]) >= args.min_mut]
    if args.limit:
        variants = variants[:args.limit]
    if not variants:
        print(f"\n{family}: no variants selected")
        return None

    chain = cfg.get("PEPTIDE_CHAIN", "B")
    nstruct, trials, stride, max_min, conv = resolve_protocol(args, cfg)
    jobs = int(args.jobs or cfg.get("JOBS", 8))

    out_dir = Path(args.out_dir or
                   Path(cfg.get("RESULTS_DIR", "results")) / "03_flexddg")
    suffix = f"_{args.tag}" if args.tag else ""
    work = out_dir / f"work_fam{family}{suffix}"
    out_dir.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    cif = Path(args.structure) if args.structure else find_structure(cfg,
                                                                     family)
    pdb = work / "model.pdb"
    chains = cif_to_pdb(cif, pdb)
    if chain not in chains:
        sys.exit(f"ERROR: chain {chain} not in {cif.name} (has {chains})")

    # the reference in the library must match the structure's peptide chain,
    # or every mutation would be numbered against the wrong sequence
    ref = variants[0].get("reference", "")
    st = gemmi.read_structure(str(pdb))
    st.setup_entities()
    pep = "".join(
        gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
        for ch in st[0] if ch.name == chain for r in ch
        if gemmi.find_tabulated_residue(r.name)
        and gemmi.find_tabulated_residue(r.name).is_amino_acid())
    if ref and pep != ref:
        sys.exit(f"ERROR: {family}: library reference {ref} does not match "
                 f"chain {chain} of {cif.name} ({pep})")

    by_n = {}
    for v in variants:
        by_n[len(v["_muts"])] = by_n.get(len(v["_muts"]), 0) + 1
    total_jobs = len(variants) * nstruct

    print(f"\n\033[1mfamily {family}\033[0m")
    print(f"  library    : {lib_path.name}   {len(variants)} variant(s)  ("
          + ", ".join(f"{k}-point {v}" for k, v in sorted(by_n.items())) + ")")
    print(f"  structure  : {cif.name}   chains {chains}   peptide {pep}")
    requested_stride = (args.stride or
                        (PROTOCOLS[args.protocol][2]
                         if args.protocol in PROTOCOLS else trials))
    print(f"  sampling   : nstruct {nstruct} x {trials} backrub trials, "
          f"stride {stride}   [{args.protocol}]")
    if requested_stride > trials:
        print(f"               (stride clamped from {requested_stride} to "
              f"the trial count: a larger stride never reaches a checkpoint "
              f"and\n                writes an empty database)")
    print(f"  jobs       : {total_jobs}  on {jobs} worker(s)")
    per_job = args.per_job_seconds or estimate_job_seconds(trials)
    est = total_jobs * per_job / jobs
    print(f"  estimate   : ~{hms(est)}   (at {per_job:.0f} s/job"
          + ("" if args.per_job_seconds else
             ", modelled from measured timings")
          + "; check with bench_flexddg_timing.py)")
    print(f"  work dir   : {work}")

    if args.dry_run:
        return None

    xml_text = fill_xml(xml_raw, chain, trials, stride, max_min, conv)
    init_flags = ("-restore_talaris_behavior -ex1 -ex2 -use_input_sc "
                  "-ignore_unrecognized_res -ignore_zero_occupancy false "
                  "-mute all")

    tasks = []
    for v in variants:
        for member in range(nstruct):
            if args.resume and (work / v["id"] / f"n{member}" /
                                "ddG.db3").is_file():
                continue
            tasks.append((v["id"], v["_muts"], member, str(pdb), xml_text,
                          chain, str(work), init_flags))
    skipped = total_jobs - len(tasks)
    if skipped:
        print(f"  resume     : {skipped} job(s) already done")

    values, failures = {}, {}
    perjob = []
    if tasks:
        started = time.time()
        done = n_failed = 0
        with ProcessPoolExecutor(max_workers=min(jobs, len(tasks))) as pool:
            futures = [pool.submit(run_job, t) for t in tasks]
            for fut in as_completed(futures):
                vid, member, ddg, err, dt = fut.result()
                done += 1
                if err:
                    failures.setdefault(vid, []).append(err)
                    n_failed += 1
                else:
                    values.setdefault(vid, []).append(ddg)
                    perjob.append((vid, member, round(ddg, 5), round(dt, 1)))
                progress(done, len(tasks), started, f"{n_failed} failed")
                if done == 20 and n_failed == done:
                    print("\n  every one of the first 20 jobs failed; "
                          "stopping.\n  last error: "
                          + list(failures.values())[-1][-1][:200])
                    pool.shutdown(cancel_futures=True)
                    break
        print()

    # jobs skipped by --resume still have their db3 on disk; read them back
    if args.resume:
        for v in variants:
            for member in range(nstruct):
                db = work / v["id"] / f"n{member}" / "ddG.db3"
                if db.is_file() and not any(
                        p[0] == v["id"] and p[1] == member for p in perjob):
                    ddg = read_ddg(db)
                    if ddg is not None:
                        values.setdefault(v["id"], []).append(ddg)
                        perjob.append((v["id"], member, round(ddg, 5), 0.0))

    # ---- aggregate -------------------------------------------------------
    rows = []
    for v in variants:
        vals = values.get(v["id"], [])
        if vals:
            mean = sum(vals) / len(vals)
            sd = (sum((x - mean) ** 2 for x in vals) / (len(vals) - 1)) ** 0.5 \
                if len(vals) > 1 else 0.0
            rows.append({
                "family": family, "id": v["id"], "sequence": v["sequence"],
                "n_mut": len(v["_muts"]), "mutations": v["mutations"],
                "ppi_exp": v.get("ppi_exp", ""), "method": "flexddg",
                "structure": cif.name,
                "ddg_binding": round(mean, 5),
                "ddg_binding_sd": round(sd, 5),
                "ddg_folding": "", "ddg_folding_sd": "",
                "n_obs": len(vals), "status": "ok"})
        else:
            rows.append({
                "family": family, "id": v["id"], "sequence": v["sequence"],
                "n_mut": len(v["_muts"]), "mutations": v["mutations"],
                "ppi_exp": v.get("ppi_exp", ""), "method": "flexddg",
                "structure": cif.name,
                "ddg_binding": "", "ddg_binding_sd": "",
                "ddg_folding": "", "ddg_folding_sd": "",
                "n_obs": 0, "status": "failed"})

    dest = out_dir / f"fam{family}{suffix}.csv"
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS)
        w.writeheader()
        w.writerows(rows)

    pj = out_dir / f"fam{family}{suffix}_perjob.csv"
    with open(pj, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "member", "ddg_binding", "seconds"])
        w.writerows(sorted(perjob))

    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"  wrote      : {dest}   {ok}/{len(rows)} variant(s) ok")
    print(f"               {pj}   {len(perjob)} ensemble member(s)")
    if failures:
        reasons = {}
        for errs in failures.values():
            for e in errs:
                reasons[e.split(":")[0]] = reasons.get(e.split(":")[0], 0) + 1
        print(f"  failures   : {sum(len(v) for v in failures.values())} "
              f"job(s) -- " + ", ".join(f"{k} x{v}"
                                        for k, v in reasons.items()))
    return dest


# ------------------------------------------------------------------ main --
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paths", default="paths.env")
    ap.add_argument("--families", default=None,
                    help="comma-separated subset of FAMILIES")
    ap.add_argument("--library", default=None,
                    help="one library CSV, implies a single family")
    ap.add_argument("--structure", default=None,
                    help="override the .cif chosen from STRUCTURES_DIR")
    ap.add_argument("--out-dir", default=None,
                    help="default: RESULTS_DIR/03_flexddg")
    ap.add_argument("--tag", default=None,
                    help="suffix for the output files and work directory, "
                         "so runs at different settings do not overwrite "
                         "each other (e.g. --tag prelim gives "
                         "fam363_prelim.csv)")
    ap.add_argument("--protocol", default="hummer",
                    choices=list(PROTOCOLS) + ["custom"],
                    help="published parameter set. 'barlow' is the original "
                         "35 x 35000; 'hummer' is the reduced 1 x 3500 from "
                         "Nat Comput Sci 2025 5:635 (default); 'custom' "
                         "takes the values from paths.env. Individual flags "
                         "below override whichever is chosen")
    ap.add_argument("--stride", type=int, default=None,
                    help="backrub_trajectory_stride. Equal to the trial "
                         "count is fastest; smaller gives intermediate ddG "
                         "checkpoints along each trajectory")
    ap.add_argument("--nstruct", type=int, default=None,
                    help="ensemble members per variant "
                         "(default: FLEXDDG_NSTRUCT)")
    ap.add_argument("--backrub-trials", type=int, default=None,
                    help="default: FLEXDDG_BACKRUB_TRIALS")
    ap.add_argument("-j", "--jobs", type=int, default=None,
                    help="parallel workers (default: JOBS)")
    ap.add_argument("--max-mut", type=int, default=None,
                    help="only variants with at most this many substitutions")
    ap.add_argument("--min-mut", type=int, default=None)
    ap.add_argument("--limit", type=int, default=None,
                    help="only the first N variants, for testing")
    ap.add_argument("--resume", action="store_true",
                    help="skip jobs whose ddG.db3 already exists")
    ap.add_argument("--per-job-seconds", type=float, default=None,
                    help="override the header time estimate. By default it "
                         "is modelled from the backrub trial count using "
                         "timings measured with bench_flexddg_timing.py")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the plan and the estimate, run nothing")
    args = ap.parse_args()

    cfg = load_env(args.paths)
    xml_path = Path(cfg["FLEXDDG_XML"])
    if not xml_path.is_file():
        sys.exit(f"ERROR: FLEXDDG_XML not found: {xml_path}")
    xml_raw = xml_path.read_text()

    families = [f.strip() for f in
                (args.families or cfg.get("FAMILIES", "")).split(",")
                if f.strip()]
    if not families:
        sys.exit("ERROR: no families; set FAMILIES or pass --families")

    print(f"protocol : flex ddG, {xml_path.name}   preset: {args.protocol}")
    if args.protocol == "hummer":
        print("           Hummer et al., Nat Comput Sci 2025, 5:635-647")
    elif args.protocol == "barlow":
        print("           Barlow et al., J Phys Chem B 2018, 122:5389-5399")
    elif args.protocol == "prelim":
        print("           preliminary settings, not a published parameter "
              "set.\n           Rerun at --protocol barlow before "
              "reporting these numbers.")
    print(f"families : {', '.join(families)}")

    started = time.time()
    written = []
    for family in families:
        dest = process_family(family, cfg, args, xml_raw)
        if dest:
            written.append(dest)

    if written:
        print(f"\ntotal {hms(time.time() - started)}")
        for d in written:
            print(f"  {d}")
        print("\nNext: python 04_merge_and_figures.py")


if __name__ == "__main__":
    main()
