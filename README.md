# k-MDLRDP

Research code for the **k-multiplex double Roman domination problem**
(`k-MLDRDP`).  The folder name follows the requested `k-MDLRDP-main`
spelling; the mathematical acronym used in the paper is `k-MLDRDP`.

The implementation follows the organization and command-line conventions of
the sibling `k-MLRDP-main` release, but all feasibility logic is rebuilt for
double Roman labels `0, 1, 2, 3`:

- label `0`: in every valid layer it has a label-3 neighbor or two distinct
  label-2 neighbors;
- label `1`: in every valid layer it has a label-2 or label-3 neighbor;
- one common labeling must be valid on at least `k` complete layers.

No source code was available in the reference release (it contains binaries,
instances, and results), so this is an independent, auditable implementation.

## Included methods

- `kmdrdp_greedy.exe`: randomized greedy construction, feasibility repair,
  and redundant-label reduction;
- `kmdrdp_ga.exe`: tournament-selection GA with elitism, one-point crossover,
  four-label mutation, and first/best-improvement local search;
- `kmdrdp_vns.exe`: VNS with Hamming shaking, first/best-improvement search,
  greedy or random initialization, and adaptive plateau acceptance;
- `python/exact_ilp.py`: exact DR-D (direct label) and DR-C (cumulative
  label) binary models solved with SciPy/HiGHS.
- `python/ga_ilp_hybrid.py`: support-aware GA--ILP V2. A complete GA incumbent
  is followed by restricted DR-D neighborhoods with strict improvement
  cuts, incumbent reinjection after success, and growth after repeated failure.

The C++ evaluator incrementally maintains, for every vertex-layer pair, the
numbers of label-2 and label-3 neighbors, the violation flag, and the
protection deficit.  A one-vertex move updates only the changed vertex and its
neighbors.

## Instance format

An instance is a directory containing `layer_1.txt`, `layer_2.txt`, and so on.
Each file has the same format as the reference repository:

```text
number_of_vertices number_of_edges
u1 v1
u2 v2
...
```

Vertex indices are one-based.  The reconstructed PPI instances and a small
two-layer smoke instance are included.  For example:

```powershell
./bin/kmdrdp_greedy.exe Instances/PPI/PPI-5Layers 5 13579 10
```

The PPI instances are under `Instances/PPI`; their identifier mappings and
edge-set validation artifacts are under `data_reconstruction/output`.  The
smoke instance is under `Instances/sample-2layer`.

## Build and test on Windows

Requirements: GCC with C++17 support, PowerShell, Python, NumPy, and SciPy.

```powershell
./build.ps1
./test.ps1
./run_examples.ps1
```

Alternatively, with MinGW Make:

```powershell
mingw32-make
mingw32-make test
```

## Commands

Greedy:

```text
kmdrdp_greedy INSTANCE_DIR K [SEED=13579] [CONSTRUCTIONS=10]
```

GA (same positional order as the reference executable):

```text
kmdrdp_ga INSTANCE_DIR K POP GENERATIONS ELITES MUTATION_PCT \
  TOURNAMENT CROSSOVER_RATE SEED LOCAL_SEARCH TIME_LIMIT \
  [--init greedy|random]
```

`LOCAL_SEARCH=1` selects first improvement; `2` selects best improvement.

VNS (same positional order as the reference executable):

```text
kmdrdp_vns vns INSTANCE_DIR K SEED KAPPA_MIN KAPPA_MAX TIE_PROB \
  LOCAL_SEARCH ITER_MAX NO_IMPROVE_MAX TIME_LIMIT \
  [--init greedy|random]
```

Exact ILP:

```powershell
python python/exact_ilp.py INSTANCE_DIR K `
  --model direct --time-limit 3600 --output result.json

python python/exact_ilp.py INSTANCE_DIR K `
  --model cumulative --time-limit 3600 --output result.json
```

Adaptive GA--ILP hybrid:

```powershell
python python/ga_ilp_hybrid.py INSTANCE_DIR K --init greedy `
  --time-limit 60 --ga-time 10 --subproblem-time 5 `
  --seed 13579 --output Results/hybrid.json
```

Batch experiments (rows are appended immediately and successful rows are
skipped on restart):

```powershell
python python/run_experiments.py PATH_TO_REFERENCE_INSTANCES `
  --methods greedy,ga-r-fi,ga-r-hybrid,vns-r-fi,vns-r-bi,ga-g-fi,ga-g-hybrid,vns-g-fi,vns-g-bi,dr-new1,dr-new2 `
  --k protocol --seeds 13579,24680 --time-limit 3600 `
  --output Results/batch.csv
```

Every solver prints or saves the complete label vector.  The exact solver
independently replays the returned solution before reporting it.

## Mapping to the manuscript

| Manuscript component | Implementation |
| --- | --- |
| DR-D | `python/exact_ilp.py --model direct` |
| DR-C | `python/exact_ilp.py --model cumulative` |
| Layer violation and deficit functions | `src/evaluator.cpp` |
| Incremental move evaluation | `EvaluationState::apply_move` |
| Greedy construction/repair/reduction | `greedy_construct`, `repair`, `reduce_feasible` |
| GA-R / GA-G | `kmdrdp_ga.exe --init random|greedy` |
| VNS-R / VNS-G, FI / BI | `kmdrdp_vns.exe --init ...`, `LOCAL_SEARCH=1|2` |
| Adaptive GA--ILP | `python/ga_ilp_hybrid.py` |

The hybrid keeps the high-throughput population search in compiled C++ and
uses SciPy/HiGHS only for restricted exact neighborhoods.  Its JSON report
records every neighborhood size, subproblem objective, acceptance decision,
and final certificate so that the adaptive rule is reproducible.

## Experiments and reproducibility

The preliminary 11-method study is retained as
`formal_experiments/screening_protocol.json`; it is not the final protocol.
The final synthetic experiment is frozen in `formal_experiments/protocol.json`:
three selected methods (Greedy, VNS-G-FI, and GA-G-Hybrid), 24 `n=1000`
instance--`k` configurations, and ten runs per configuration (the gate seed
plus nine new seeds), giving 720 runs.  The complete
run ledger, per-configuration summary, and validation report are:

- `Results/formal/n1000_formal_10seed.csv` (720 rows);
- `Results/formal/n1000_formal_by_config.csv`;
- `Results/formal/n1000_formal_summary.json` (720 label vectors replayed,
  zero validation errors);
- `Results/formal/n1000_sensitivity_audit.json` (gate-seed-excluded and
  nine-underlying-graph blocked sensitivity analyses).

The post-experiment environment audit is in
`formal_experiments/environment_audit.json`. It records the workstation,
compiler, Python/HiGHS versions, and the fact that the archived ILP launchers
did not explicitly pin solver thread counts. The synthetic edge lists are
fixed inputs from the reference release; their original graph-generation seeds
are not available.

The `n=30` method screening and the `n=100` exact-model experiment remain
available in `Results/formal/n30_runs.csv`, `Results/formal/n100_runs.csv`,
`Results/formal/n100_exact_dr_new1.csv`, and
`Results/formal/n100_exact_dr_new2.csv`.  These files support the algorithm
selection and the DR-D/DR-C comparison; they are not part of the final
720-run heuristic table.

`python/analyze_n30_optimal_gaps.py` derives gaps to the certified `n=30`
optima from the archived ledger. Its output,
`Results/formal/n30_optimal_gap_summary.json`, includes the source-ledger
checksum, 720 independently replayed heuristic certificates, 48 exact
certificates, and configuration-level gap summaries. Historical raw filenames
and solver identifiers retain `dr-new1`/`dr-new2`; these correspond to DR-D/DR-C.
To replay the certificates, supply the `Instances` directory from the
reference `k-MLRDP-main` release, for example:

```powershell
python python/analyze_n30_optimal_gaps.py --reference-instances ../k-MLRDP-main/Instances
```

### Real PPI optimization experiments

The reconstructed three-, four-, and five-layer PPI instances are under
`Instances/PPI`.  Their layer order and identifier-resolution requirements
are recorded in `formal_experiments/ppi_gate_protocol.json`.  The optimization
evidence consists of:

- the 36-run one-seed gate:
  `Results/formal/ppi_gate_merged.csv`;
- the 30-run, ten-run GA-G-Hybrid extension at `k=1` (one reused gate seed
  plus nine new seeds per instance):
  `Results/formal/ppi_formal_ga_g_hybrid_k1_10seeds.csv`;
- the deduplicated union of these experiments:
  `Results/formal/ppi_optimization_63run_ledger.csv`.

The two source tables contain three identical overlapping rows (the seed
13579 GA-G-Hybrid run at `k=1` for each PPI network), so the union contains
`36 + 30 - 3 = 63` distinct runs.  The merge manifest records the keys,
provenance, path normalization, and validation checks in
`Results/formal/ppi_optimization_63run_ledger_manifest.json`.  Rebuild it with:

```powershell
python python/merge_ppi_optimization_ledger.py
```

The later paired `k=1` PPI experiment matches Greedy, VNS-G-FI, and
GA-G-Hybrid on ten seeds for each of the three PPI networks: 90 method runs.
The frozen protocol and runner are
`formal_experiments/ppi_paired_k1_protocol.json` and
`python/run_ppi_paired_k1.py`. The directory
`Results/paired_ppi/k1_10seeds_20261013/` contains the combined 90-run
ledger, run-level solution certificates, validation status, and paired summary.
It reuses 36 archived runs and adds 54 Greedy/VNS runs. Together with the
36-run gate, this is 117 distinct runs; the separate nine new all-layer runs
bring the double-Roman PPI optimization total across these stages to 126.

### Exploratory VNS ablation

The small mechanism ablation is separate from the frozen formal comparison.
Its six selected configurations, three variants, three search seeds, and
60-second search budget are specified in
`formal_experiments/vns_ablation_protocol.json`. Build and test the variants
with `formal_experiments/build_vns_ablation.ps1`; the runner is
`python/run_vns_ablation.py`. The 54 validated run certificates, initial
labelings, evaluator benchmark, and summary are under
`Results/ablation/vns_small_20261013/`. The variants isolate repair-deficit
ranking and adaptive tie probability. These results are exploratory rather
than a new formal method comparison.

The all-layer biological analysis uses a ten-run VNS-G-FI ledger comprising
the gate seed and nine new seeds on the five-layer network with `k=M=5`.
Its frozen protocol and
validated outputs are `formal_experiments/ppi_5layer_k5_vns_protocol.json`,
`Results/formal/ppi_5layer_k5_vns_10seeds.csv`, and
`Results/formal/ppi_5layer_k5_vns_summary.json`.

### Biological and ordinary-Roman analyses

Machine-readable biological outputs are collected in `Results/biology`:

- `ppi5_k5_protein_consensus.tsv` and its manifest define strict (10/10) and
  stable (at least 8/10) double-Roman label sets;
- `ppi5_functional_enrichment.csv` and `ppi5_reactome_enrichment.csv` contain
  GO/g:Profiler and Reactome results, with request/response and protocol files;
- `gprofiler_root_excluded_comparison.json` records the ordinary/double-Roman
  functional overlap after excluding the nonbiological KEGG root term;
- `ppi5_essential_protein_tests.csv` contains essential-protein tests;
- `ppi5_degree_matched_centrality.csv` and
  `ppi5_degree_matched_diagnostics.json` contain the 10,000-replicate
  degree-matched centrality randomization;
- `ppi5_degree_matched_go.csv` and `ppi5_degree_matched_go_protocol.json`
  contain the degree-matched GO analysis and its full annotation/matching
  specification;
- `ppi5_degree_matched_go_sparse_summary.json` reports a sensitivity check
  requiring at least five annotated candidate proteins per reported term,
  while retaining the original Benjamini--Hochberg correction universe;
- files prefixed `ppi5_roman_` contain the ordinary-Roman control and its
  direct comparison with the double-Roman roles.

The final ordinary-Roman control uses the same compiled VNS framework in
`--mode roman`, which restricts labels to `{0,1,2}` and enforces the ordinary
Roman rule on every selected layer.  Its frozen ten-seed protocol is
`formal_experiments/certified_roman_vns_protocol.json`.  Complete label
certificates and the independent replay summary are stored in
`Results/formal/ppi_5layer_k5_roman_vns_certified_10seeds.csv` and
`Results/formal/ppi_5layer_k5_roman_vns_certified_summary.json`.  All ten
certificates pass replay on all five layers.  The resulting ordinary label-2
consensus contains 394 proteins at 10/10 and 459 proteins at 8/10.

The corresponding analysis programs are under `python/`, including
`prepare_ppi_biology.py`, `analyze_ppi_centrality_degree_matched.py`,
`analyze_ppi_essential_proteins.py`, `analyze_degree_matched_go.py`,
`prepare_roman_ppi_consensus.py`, `analyze_roman_double_comparison.py`, and
`summarize_degree_matched_go_sparse_terms.py`.
Derived networks, identifier mappings, validation reports, and source
provenance are retained in `data_reconstruction`; network construction is
therefore separated from the optimization and downstream statistical
analyses.  Large third-party downloads are intentionally not redistributed.
Their versions, commit, checksums, and acquisition limitations are documented
in `data_reconstruction/README.md`.

All CSV ledgers store the complete label vector.  Summary JSON files document
row counts and independent feasibility replay; protocol JSON files freeze the
method, seed, time limit, and network/layer choices used for each reported
experiment.

Machine-specific absolute paths in the newly released ablation and paired-PPI
JSON records have been converted to relative paths or `REFERENCE_INSTANCES`
placeholders. Their objective values, label vectors, seeds, and timings are
unchanged. SHA-256 fields in the frozen manifests identify the original
experimental inputs and build artifacts, so they should not be interpreted as
checksums of these path-normalized release JSON files. Temporary logs,
process IDs, compiled executables, and third-party source downloads are not
redistributed.

