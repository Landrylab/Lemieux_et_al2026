#!/usr/bin/env python3
"""
02_foldx_ddg.py -- binding and folding ddG for every peptide variant, by
FoldX.

===========================================================================
TL;DR
===========================================================================
    in   results/01_libraries/fam<N>.csv        (from 01)
         structures/**/fold_<N>*model_<i>.cif   (i = MODEL_INDEX)
    out  results/02_foldx/fam<N>.csv            one row per variant
         results/02_foldx/work_fam<N>/          the FoldX audit trail

    run  python 02_foldx_ddg.py                  every family
         python 02_foldx_ddg.py --families 363   one family
         python 02_foldx_ddg.py --max-mut 1      single mutants only
         python 02_foldx_ddg.py --dry-run        plan and cost, no compute
         python 02_foldx_ddg.py --resume         after an interruption

Positive ddG means destabilising, i.e. weaker binding. Same convention and
the same output columns as 03_flexddg_ddg.py, so the two concatenate.

===========================================================================
WHAT IT RUNS
===========================================================================
Three FoldX 5 commands, in order:

  RepairPDB        FOLDX_REPAIR_ROUNDS times on the input structure, each
                   round taking the previous output. FoldX energies are
                   meaningless on an unrepaired structure, so this step is
                   not optional.

                   NOTE ON THE ROUND COUNT. The FoldX manual specifies one
                   round, and MutateX, the most widely used FoldX wrapper,
                   also runs exactly one. Repeating it is a community
                   practice rather than a documented protocol. Five was
                   chosen here from the convergence curve measured on these
                   structures, where the total energy plateaus by round
                   four; every round's total is written to
                   repair_energies.csv so the choice can be checked rather
                   than taken on trust. Set FOLDX_REPAIR_ROUNDS=1 to follow
                   the manual.

  BuildModel       --vdwDesign controls how harshly a new side chain is
                   penalised for not fitting. The FoldX manual recommends
                   level 2, the softest, when introducing mutations, which
                   is what BuildModel does; the binary defaults to 0.

                   This matters for the comparison this repository makes.
                   FoldX repacks on a fixed backbone, so a bulky
                   substitution has nowhere to go and level 0 punishes it
                   twice: once for the real clash and once for a clash a
                   flexible backbone would have relieved. Running FoldX at
                   level 0 while comparing it against a backbone-sampling
                   method would stack the comparison. Level 2 is therefore
                   the default here, and --vdw-design 0 reproduces the
                   binary's own default if you want to see the difference.

                   The call itself is one per chunk of variants, driven by an
                   individual_list.txt in which each line is one variant.
                   Substitutions within a line are comma separated and are
                   introduced together. FoldX rebuilds BOTH the mutant and a
                   matched wild type for every run, so the pairing cancels
                   most of the repacking noise. Folding ddG comes from
                   Dif_*.fxout.

  AnalyseComplex   on each matched mutant/wild-type pair, separating
                   RECEPTOR_CHAIN from PEPTIDE_CHAIN. Binding ddG is the
                   difference of the two interaction energies. Skipped
                   unless --binding, since it roughly doubles the runtime.

Everything is repacked on a FIXED BACKBONE. That is the defining limitation
of this method and the reason 03 exists: flex ddG samples backbone
accommodation, FoldX does not. Expect the gap to show up on multi-point
variants, where a frozen backbone has the least room to absorb several
substitutions at once.

===========================================================================
DIFFERENCES FROM 03 THAT ARE REAL, NOT BOOKKEEPING
===========================================================================
  * FoldX requires RepairPDB; flex ddG does not and minimises internally.
    The two arms therefore start from slightly different coordinates. This
    is stated in the paper rather than hidden.
  * FoldX gives folding AND binding ddG; InterfaceDdGMover gives binding
    only, so ddg_folding is blank in 03 and populated here.
  * FoldX is close to deterministic per structure, so FOLDX_NRUNS is small
    (3) and ddg_binding_sd is correspondingly tiny. It measures repacking
    scatter, not structural uncertainty.

===========================================================================
ROBUSTNESS
===========================================================================
FoldX occasionally hangs on a structure it cannot resolve. Each call runs in
its own process group and is killed by --timeout, group and all: a plain
subprocess timeout leaves a descendant holding the pipe open and the worker
waits forever. Work is split into chunks of --chunk-size variants so a hang
costs one chunk rather than the whole family, and --resume reuses any chunk
whose Dif file is already complete.

===========================================================================
OUTPUT
===========================================================================
    family, id, sequence, n_mut, mutations, ppi_exp, method, structure,
    ddg_binding, ddg_binding_sd, ddg_folding, ddg_folding_sd, n_obs, status

    status   ok | failed -- failed variants keep their row with blank ddG,
             so the row count always matches the library

===========================================================================
REQUIRES
===========================================================================
FoldX 5.x    free academic licence from https://foldxsuite.crg.eu/, not
             bundled. Schymkowitz et al., NAR 2005; Delgado et al.,
             Bioinformatics 2019
gemmi        mmCIF to PDB. Wojdyr, JOSS 2022
"""

import argparse
import csv
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    import gemmi
except ImportError:
    sys.exit("ERROR: gemmi not found.  pip install gemmi")

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


# -------------------------------------------------------------- foldx io --
def run_foldx(binary, args_list, cwd, timeout):
    """Run one FoldX command, killing the whole process group if it hangs.

    subprocess.run's timeout kills only the direct child. FoldX can leave a
    descendant holding the pipe open, so the parent blocks forever waiting
    on output and the worker never returns. Its own session plus a group
    kill ends that.
    """
    proc = subprocess.Popen([str(binary)] + args_list, cwd=cwd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT,
                            text=True, start_new_session=True)
    try:
        out, _ = proc.communicate(timeout=timeout)
        return proc.returncode == 0, out or ""
    except subprocess.TimeoutExpired:
        try:
            os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            proc.kill()
        try:
            out, _ = proc.communicate(timeout=30)
        except subprocess.TimeoutExpired:
            out = ""
        return False, (out or "") + f"\n[killed after {timeout}s]"


def parse_dif(path):
    """Dif_*.fxout -> [(pdb_name, total_ddg)] in file order."""
    rows, started = [], False
    for line in Path(path).read_text(errors="ignore").splitlines():
        parts = line.rstrip("\n").split("\t")
        if not started:
            if parts and parts[0].strip().lower().startswith("pdb"):
                started = True
            continue
        if len(parts) < 2 or not parts[0].strip():
            continue
        try:
            rows.append((parts[0].strip(), float(parts[1])))
        except ValueError:
            continue
    return rows


def parse_interaction(path, chains):
    """Interaction energy for one chain pair out of Interaction_*.fxout."""
    want = set(chains)
    for line in Path(path).read_text(errors="ignore").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 6 and {parts[1].strip(), parts[2].strip()} == want:
            try:
                return float(parts[5])
            except ValueError:
                return None
    return None


def variant_index(pdb_name):
    """model_Repair_12_1.pdb -> 12, the line number in individual_list."""
    m = re.search(r"_(\d+)(?:_(\d+))?\.pdb$", pdb_name)
    return int(m.group(1)) if m else None


# --------------------------------------------------------------- workers --
def parse_total_energy(stdout, wd, stem):
    """Total energy after a repair round.

    FoldX 5 prints the repaired structure's energy breakdown to stdout, not
    to the .fxout, which holds the per-residue table. stdout also carries a
    running "Total = <n>" from every optimisation cycle, so the value that
    matters is the last one before "End of Repair" -- taking the first, or
    the last line of the file, gives an intermediate.

    The .fxout is still checked as a fallback for builds that write the
    total there instead.
    """
    if stdout:
        head = stdout.split("End of Repair")[0]
        totals = re.findall(r"^\s*Total\s*=\s*(-?[\d.]+)", head, re.M)
        if totals:
            try:
                return float(totals[-1])
            except ValueError:
                pass

    for name in (f"{stem}_Repair.fxout", "Raw_" + f"{stem}_Repair.fxout"):
        f = wd / name
        if not f.is_file():
            continue
        text = f.read_text(errors="ignore")
        m = re.search(r"^\s*Total\s*=?\s*(-?[\d.]+)", text, re.M)
        if m:
            try:
                return float(m.group(1))
            except ValueError:
                pass
        lines = [ln for ln in text.splitlines() if ln.strip()]
        for i, ln in enumerate(lines):
            if ln.split("\t")[0].strip().lower().startswith("pdb"):
                if i + 1 < len(lines):
                    cells = lines[i + 1].split("\t")
                    if len(cells) > 1:
                        try:
                            return float(cells[1])
                        except ValueError:
                            pass
                break
    return None


def repair_structure(task):
    """RepairPDB N times, feeding each output back in. One worker.

    Every round's structure and .fxout are kept, in repair/round<N>/. FoldX
    appends _Repair to the filename each time, so a naive loop would produce
    model_Repair_Repair_Repair...; instead each round's products are moved
    into their own folder and the input is renamed back, which keeps the
    filenames stable and leaves a complete per-round record.

    The record is what justifies the round count: the total energies go into
    repair_energies.csv for the convergence figure, and the structures are
    there if anyone wants to check that the coordinates settled too, not
    just the energy.
    """
    workdir, binary, rounds, timeout, vdw = task
    wd = Path(workdir).resolve()
    keep = wd / "repair"
    keep.mkdir(exist_ok=True)
    energies = []
    try:
        for rnd in range(1, rounds + 1):
            if rnd == 1:
                src = "model.pdb"
            else:
                (wd / "model_Repair.pdb").rename(wd / "model_prev.pdb")
                src = "model_prev.pdb"

            ok, msg = run_foldx(binary, ["--command=RepairPDB",
                                         f"--pdb={src}",
                                         f"--vdwDesign={vdw}"], wd, timeout)
            produced = wd / f"{Path(src).stem}_Repair.pdb"
            if not ok or not produced.is_file():
                return None, f"RepairPDB round {rnd} failed: {msg[-300:]}"

            total = parse_total_energy(msg, wd, Path(src).stem)
            energies.append(total)

            # archive this round before the next one overwrites anything
            rd = keep / f"round{rnd}"
            rd.mkdir(exist_ok=True)
            shutil.copy(produced, rd / "model_Repair.pdb")
            for pattern in ("*.fxout", "*.txt"):
                for f in wd.glob(pattern):
                    if f.is_file():
                        shutil.copy(f, rd / f.name)
            (rd / "foldx.log").write_text(msg)

            if rnd > 1:
                produced.rename(wd / "model_Repair.pdb")
                (wd / "model_prev.pdb").unlink(missing_ok=True)

        return energies, None
    except Exception as exc:                                 # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"


def build_chunk(task):
    """BuildModel (+ AnalyseComplex) on one slice of the library."""
    (chunk_id, offset, lines, struct_dir, binary, chains, nruns,
     do_binding, timeout, resume, labels, vdw) = task
    wd = (Path(struct_dir) / f"chunk_{chunk_id:03d}").resolve()
    try:
        expected = len(lines) * nruns
        dif = wd / "Dif_model_Repair.fxout"
        if resume and dif.is_file() and len(parse_dif(dif)) >= expected:
            folding = [(offset + variant_index(n) - 1, d)
                       for n, d in parse_dif(dif)
                       if variant_index(n) is not None]
            binding = {}
            if do_binding:
                binding = read_binding(wd, offset, chains)
            return chunk_id, folding, binding, None, True

        if wd.exists():
            shutil.rmtree(wd)
        wd.mkdir(parents=True)
        shutil.copy(Path(struct_dir).resolve() / "model_Repair.pdb",
                    wd / "model_Repair.pdb")
        (wd / "individual_list.txt").write_text("\n".join(lines) + "\n")

        # FoldX names its outputs by line number, so record what each line
        # is; without this the artifacts are uninterpretable later
        with open(wd / "variant_map.tsv", "w") as fh:
            fh.write("foldx_index\tid\tsequence\tmutations\n")
            for j, (vid, seq, muts) in enumerate(labels, 1):
                fh.write(f"{j}\t{vid}\t{seq}\t{muts}\n")

        ok, msg = run_foldx(binary, [
            "--command=BuildModel", "--pdb=model_Repair.pdb",
            "--mutant-file=individual_list.txt",
            f"--numberOfRuns={nruns}",
            f"--out-pdb={'true' if do_binding else 'false'}",
            f"--vdwDesign={vdw}"], wd, timeout)
        (wd / "foldx.log").write_text(msg)

        if not dif.is_file():
            found = list(wd.glob("Dif_*.fxout"))
            if not found:
                return chunk_id, [], {}, f"no Dif file: {msg[-200:]}", False
            dif = found[0]

        folding = [(offset + variant_index(n) - 1, d)
                   for n, d in parse_dif(dif)
                   if variant_index(n) is not None]
        binding = {}
        if do_binding:
            for name, _ in parse_dif(dif):
                idx = variant_index(name)
                mut, wt = wd / name, wd / f"WT_{name}"
                if idx is None or not (mut.is_file() and wt.is_file()):
                    continue
                pair = []
                for pdb in (wt, mut):
                    run_foldx(binary, [
                        "--command=AnalyseComplex", f"--pdb={pdb.name}",
                        f"--analyseComplexChains={','.join(chains)}"],
                        wd, timeout)
                    ix = wd / f"Interaction_{pdb.stem}_AC.fxout"
                    pair.append(parse_interaction(ix, chains)
                                if ix.is_file() else None)
                if None not in pair:
                    binding.setdefault(offset + idx - 1, []).append(
                        pair[1] - pair[0])
        return chunk_id, folding, binding, None, False

    except Exception as exc:                                 # noqa: BLE001
        return chunk_id, [], {}, f"{type(exc).__name__}: {exc}", False


def read_binding(wd, offset, chains):
    """Recover binding values from a chunk that was already computed."""
    out = {}
    for ix in wd.glob("Interaction_WT_model_Repair_*_AC.fxout"):
        stem = ix.name[len("Interaction_WT_"):-len("_AC.fxout")]
        idx = variant_index(stem + ".pdb")
        mut_ix = wd / f"Interaction_{stem}_AC.fxout"
        if idx is None or not mut_ix.is_file():
            continue
        wt = parse_interaction(ix, chains)
        mu = parse_interaction(mut_ix, chains)
        if wt is not None and mu is not None:
            out.setdefault(offset + idx - 1, []).append(mu - wt)
    return out


# ------------------------------------------------------------- per family --
def load_library(path):
    rows = []
    with open(path) as fh:
        for row in csv.DictReader(fh):
            muts = [t.strip() for t in (row.get("mutations") or "").split(",")
                    if t.strip()]
            if muts:
                row["_muts"] = muts
                rows.append(row)
    return rows


def find_structure(cfg, family):
    root = Path(cfg.get("STRUCTURES_DIR", "structures"))
    idx = cfg.get("MODEL_INDEX", "0")
    hits = sorted(root.rglob(f"*{family}*model_{idx}.cif"))
    if not hits:
        sys.exit(f"ERROR: no *{family}*model_{idx}.cif under {root}")
    return hits[0]


def peptide_sequence(pdb, chain):
    st = gemmi.read_structure(str(pdb))
    st.setup_entities()
    letters = []
    for ch in st[0]:
        if ch.name != chain:
            continue
        for r in ch:
            info = gemmi.find_tabulated_residue(r.name)
            if info and info.is_amino_acid():
                c = info.one_letter_code.upper()
                if c.isalpha() and c != "X":
                    letters.append(c)
    return "".join(letters)


class Progress:
    """Live progress bar.

    A chunk of variants takes many minutes, so a bar that only redraws when
    a chunk finishes looks frozen. A background thread redraws every second,
    and the fraction is estimated from the .fxout files FoldX writes as it
    goes rather than from completed chunks alone -- that gives movement
    inside a chunk, which is where nearly all the time is spent.
    """

    def __init__(self, total_variants, workdir, jobs):
        self.total = max(1, total_variants)
        self.workdir = Path(workdir)
        self.jobs = jobs
        self.done_chunks = 0
        self.total_chunks = 0
        self.failed = 0
        self.started = time.time()
        self._stop = threading.Event()
        self._thread = None

    def _count_started(self):
        """Variants FoldX has produced output for, across every chunk."""
        n = 0
        for dif in self.workdir.glob("chunk_*/Dif_*.fxout"):
            try:
                n += sum(1 for line in dif.open()
                         if line.startswith("model_Repair_"))
            except OSError:
                pass
        return n

    def _draw(self):
        seen = self._count_started()
        frac = min(1.0, seen / self.total) if self.total else 0.0
        bar = ("#" * int(28 * frac)).ljust(28, ".")
        elapsed = time.time() - self.started
        eta = elapsed / frac - elapsed if frac > 0.01 else 0
        tail = f"{self.failed} failed  " if self.failed else ""
        sys.stdout.write(
            f"\r  [{bar}] {frac * 100:5.1f}%  "
            f"chunk {self.done_chunks}/{self.total_chunks}  {tail}"
            f"elapsed {hms(elapsed)}  ETA "
            f"{hms(eta) if eta else '--'}      ")
        sys.stdout.flush()

    def start(self, total_chunks):
        self.total_chunks = total_chunks
        self._thread = threading.Thread(target=self._loop, daemon=True)
        self._thread.start()

    def _loop(self):
        while not self._stop.wait(1.0):
            self._draw()

    def finish(self):
        self._stop.set()
        if self._thread:
            self._thread.join(timeout=2)
        self._draw()
        print()


def mutation_line(muts, chain):
    """['S1A','G2D'] -> 'SB1A,GB2D;' -- FoldX wants the chain inside."""
    out = []
    for tok in muts:
        m = re.match(r"^([A-Z])(\d+)([A-Z])$", tok)
        if not m:
            return None
        out.append(f"{m.group(1)}{chain}{m.group(2)}{m.group(3)}")
    return ",".join(out) + ";"


def process_family(family, cfg, args):
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

    binary = Path(cfg["FOLDX_BINARY"]).expanduser()
    if not binary.is_file():
        sys.exit(f"ERROR: FOLDX_BINARY not found: {binary}")
    binary = binary.resolve()

    chain = cfg.get("PEPTIDE_CHAIN", "B")
    pair = [cfg.get("RECEPTOR_CHAIN", "A"), chain]
    rounds = int(args.repair_rounds or cfg.get("FOLDX_REPAIR_ROUNDS", 5))
    nruns = int(args.nruns or cfg.get("FOLDX_NRUNS", 3))
    jobs = int(args.jobs or cfg.get("JOBS", 8))

    out_dir = Path(args.out_dir or
                   Path(cfg.get("RESULTS_DIR", "results")) / "02_foldx")
    work = out_dir / f"work_fam{family}"
    out_dir.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    cif = Path(args.structure) if args.structure else find_structure(cfg,
                                                                     family)
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    st.remove_ligands_and_waters()
    st.write_pdb(str(work / "model.pdb"))
    chains = [c.name for c in st[0]]
    for c in pair:
        if c not in chains:
            sys.exit(f"ERROR: chain {c} not in {cif.name} (has {chains})")

    pep = peptide_sequence(work / "model.pdb", chain)
    ref = variants[0].get("reference", "")
    if ref and pep != ref:
        sys.exit(f"ERROR: {family}: library reference {ref} does not match "
                 f"chain {chain} of {cif.name} ({pep})")

    lines, labels, bad = [], [], []
    for v in variants:
        line = mutation_line(v["_muts"], chain)
        if line is None:
            bad.append(v["id"])
            continue
        lines.append(line)
        labels.append((v["id"], v["sequence"], v["mutations"]))
    variants = [v for v in variants if v["id"] not in set(bad)]

    by_n = {}
    for v in variants:
        by_n[len(v["_muts"])] = by_n.get(len(v["_muts"]), 0) + 1
    # aim for about four chunks per worker: enough waves that the bar and
    # the ETA mean something, few enough that FoldX startup stays negligible
    chunk = args.chunk_size or max(5, -(-len(variants) // max(1, jobs * 4)))
    n_chunks = -(-len(variants) // chunk)

    print(f"\n\033[1mfamily {family}\033[0m")
    print(f"  library    : {lib_path.name}   {len(variants)} variant(s)  ("
          + ", ".join(f"{k}-point {v}" for k, v in sorted(by_n.items())) + ")")
    print(f"  structure  : {cif.name}   chains {chains}   peptide {pep}")
    print(f"  protocol   : RepairPDB x{rounds}, BuildModel x{nruns} run(s), "
          f"vdwDesign {args.vdw_design}"
          + (", AnalyseComplex" if args.binding else ", no binding"))
    if rounds != 1:
        print(f"               (the FoldX manual and MutateX both use one "
              f"repair round; see repair_energies.csv)")
    print(f"  chunks     : {n_chunks} of <={chunk} variant(s) "
          f"on {jobs} worker(s)")
    print(f"  work dir   : {work}")
    if bad:
        print(f"  skipped    : {len(bad)} variant(s) with unparsable "
              f"mutations, e.g. {bad[:3]}")
    if args.dry_run:
        return None

    # ---- repair ---------------------------------------------------------
    repaired = work / "model_Repair.pdb"
    if args.resume and repaired.is_file():
        print(f"  repair     : reusing {repaired.name}")
        energies = []
    else:
        t0 = time.time()
        sys.stdout.write(f"  repair     : {rounds} round(s) ...")
        sys.stdout.flush()
        energies, err = repair_structure((str(work), str(binary), rounds,
                                          args.timeout, args.vdw_design))
        if err:
            print(f"\n  FAILED during repair: {err}")
            return None
        print(f"\r  repair     : {rounds} round(s) in {hms(time.time() - t0)}"
              f"   total energy "
              + " -> ".join("?" if e is None else f"{e:.1f}"
                            for e in energies))
        print(f"               every round kept in "
              f"{work / 'repair'}/round<N>/")
        with open(work / "repair_energies.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["round", "total_energy", "delta_from_previous"])
            for i, e in enumerate(energies, 1):
                prev = energies[i - 2] if i > 1 else None
                delta = ("" if e is None or prev is None
                         else round(e - prev, 4))
                w.writerow([i, "" if e is None else round(e, 4), delta])

    # ---- build ----------------------------------------------------------
    tasks = []
    for cid, start in enumerate(range(0, len(lines), chunk)):
        tasks.append((cid, start, lines[start:start + chunk], str(work),
                      str(binary), pair, nruns, args.binding, args.timeout,
                      args.resume, labels[start:start + chunk],
                      args.vdw_design))

    fold_vals, bind_vals, failed = {}, {}, 0
    bar = Progress(len(variants) * nruns, work, jobs)
    bar.start(len(tasks))
    with ProcessPoolExecutor(max_workers=min(jobs, len(tasks))) as pool:
        futures = [pool.submit(build_chunk, t) for t in tasks]
        for fut in as_completed(futures):
            cid, folding, binding, err, reused = fut.result()
            bar.done_chunks += 1
            if err:
                failed += 1
                bar.failed = failed
            for gidx, val in folding:
                fold_vals.setdefault(gidx, []).append(val)
            for gidx, vals in binding.items():
                bind_vals.setdefault(gidx, []).extend(vals)
    bar.finish()

    # ---- aggregate -------------------------------------------------------
    def stats(vals):
        if not vals:
            return "", ""
        mean = sum(vals) / len(vals)
        if len(vals) == 1:
            return round(mean, 5), 0.0
        sd = (sum((x - mean) ** 2 for x in vals) / (len(vals) - 1)) ** 0.5
        return round(mean, 5), round(sd, 5)

    rows = []
    for i, v in enumerate(variants):
        f_mean, f_sd = stats(fold_vals.get(i, []))
        b_mean, b_sd = stats(bind_vals.get(i, []))
        got = bool(fold_vals.get(i))
        rows.append({
            "family": family, "id": v["id"], "sequence": v["sequence"],
            "n_mut": len(v["_muts"]), "mutations": v["mutations"],
            "ppi_exp": v.get("ppi_exp", ""), "method": "foldx",
            "structure": cif.name,
            "ddg_binding": b_mean, "ddg_binding_sd": b_sd,
            "ddg_folding": f_mean, "ddg_folding_sd": f_sd,
            "n_obs": len(fold_vals.get(i, [])),
            "status": "ok" if got else "failed"})

    dest = out_dir / f"fam{family}.csv"
    with open(dest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLUMNS)
        w.writeheader()
        w.writerows(rows)

    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"  wrote      : {dest}   {ok}/{len(rows)} variant(s) ok")
    if failed:
        print(f"  failures   : {failed} chunk(s); rerun with --resume")
    return dest


# ------------------------------------------------------------------ main --
def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--paths", default="paths.env")
    ap.add_argument("--families", default=None)
    ap.add_argument("--library", default=None,
                    help="one library CSV, implies a single family")
    ap.add_argument("--structure", default=None)
    ap.add_argument("--out-dir", default=None,
                    help="default: RESULTS_DIR/02_foldx")
    ap.add_argument("--repair-rounds", type=int, default=None,
                    help="default: FOLDX_REPAIR_ROUNDS")
    ap.add_argument("--nruns", type=int, default=None,
                    help="BuildModel runs per variant (default: FOLDX_NRUNS)")
    ap.add_argument("-j", "--jobs", type=int, default=None)
    ap.add_argument("--vdw-design", type=int, default=2, choices=[0, 1, 2],
                    help="FoldX van der Waals design level. 2 is the "
                         "manual's recommendation for mutation work and the "
                         "default here; 0 is the binary's own default and "
                         "penalises new side chains hardest (default: 2)")
    ap.add_argument("--chunk-size", type=int, default=None,
                    help="variants per BuildModel call. Default is about "
                         "four chunks per worker: smaller chunks give finer "
                         "progress and lose less to a hang, at the cost of "
                         "repeating FoldX startup more often")
    ap.add_argument("--binding", action="store_true", default=True,
                    help="run AnalyseComplex for binding ddG (default on)")
    ap.add_argument("--no-binding", dest="binding", action="store_false",
                    help="folding ddG only, roughly halves the runtime")
    ap.add_argument("--max-mut", type=int, default=None)
    ap.add_argument("--min-mut", type=int, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--timeout", type=int, default=3600,
                    help="seconds before a FoldX call is killed (default: "
                         "3600)")
    ap.add_argument("--resume", action="store_true",
                    help="reuse the repaired structure and any complete chunk")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cfg = load_env(args.paths)
    families = [f.strip() for f in
                (args.families or cfg.get("FAMILIES", "")).split(",")
                if f.strip()]
    if not families:
        sys.exit("ERROR: no families; set FAMILIES or pass --families")

    print(f"method   : FoldX  {cfg.get('FOLDX_BINARY')}")
    print(f"families : {', '.join(families)}")

    started = time.time()
    written = []
    for family in families:
        dest = process_family(family, cfg, args)
        if dest:
            written.append(dest)

    if written:
        print(f"\ntotal {hms(time.time() - started)}")
        for d in written:
            print(f"  {d}")
        print("\nNext: python 03_flexddg_ddg.py")


if __name__ == "__main__":
    main()
