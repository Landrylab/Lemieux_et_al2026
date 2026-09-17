# release

Computational analysis accompanying [citation].

## Layout

    scripts/     the pipeline, numbered in the order it is run
    data/        experiment tables and reference peptide sequences
    structures/  AlphaFold3 models
    ref-pdbs/    experimental structures used for validation
    results/     per-variant predictions and the analysis

## Running it

    cp scripts/paths.template.env paths.env   # then edit the paths
    python scripts/01_extract_libraries.py
    python scripts/02_foldx_ddg.py
    python scripts/03_flexddg_ddg.py
    python scripts/04_merge_and_figures.py

FoldX and PyRosetta are not included and are not redistributable. Both are
free for academic use under their own licences and must be obtained
separately.

## What is not here

The FoldX and PyRosetta working directories are omitted: they hold every
intermediate structure and score file and run to several gigabytes. Every
value they contain that bears on a conclusion is already aggregated into
the CSVs under results/. Re-deriving an individual prediction from its raw
output requires rerunning the corresponding step.

## The file to start from

`results/04_analysis/all_ddg.csv` -- one row per variant per method, with
the measured interaction score alongside. No licensed software is needed to
read or replot it. See the README in that folder for the columns.
