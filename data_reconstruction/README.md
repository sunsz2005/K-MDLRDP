# PPI data reconstruction

This directory reconstructs the PPI instances used in *Solving Roman
domination in multiplex networks*.  The reconstruction follows the paper's
sequence:

1. normalise all yeast identifiers to systematic ORF names;
2. take the common vertex set of the selected layers;
3. preserve/recover the paper's nested integer numbering;
4. emit one named and one numeric edge list per layer;
5. compare every recoverable edge against the published instances.

## Version evidence

The article reports graph sizes but not database release numbers, download
dates, source filenames, score thresholds, or identifier conversion tables.
The following sources were therefore selected by exact graph fingerprints,
not merely by database name:

| Layer | Recovered source | Original simple graph fingerprint | Status |
|---|---|---:|---|
| STRING | `4932.protein.links.v9.1.txt` | 6,400 vertices; 777,589 undirected edges | exact |
| WI-PHI | author repository `wiphi.txt` | 5,955 vertices; 49,607 non-loop edges | exact |
| BioGRID | author repository `BioGRID.txt` | 5,640 vertices; 59,748 edges | exact |
| iRefIndex | paper gives no release; release 13 is a candidate, not proven | 6,090 vertices; 111,569 simple edges | raw named file still missing |
| MINT | paper gives no release | 4,466 vertices; 18,305 simple edges | raw named file still missing |

The three exact source files come from the paper author's public
`SupportnessOfPPIComplexes` repository (commit
`a30eb7e54490315e00e6ef941f5a7cefccd0e9fd`). Clone that repository at the
stated commit and place its contents under
`data_reconstruction/sources/SupportnessOfPPIComplexes` before rerunning the
reconstruction. The source data are not redistributed here. Their SHA-256
digests are:

```text
STRING   98BF27A0DE7F2F3AE8311079E62E869E8FC9373ED46AF57F6B8C8DDEB33649E9
WI-PHI   AFE529B3295C383C7A0D2B9BDEF156D1A1A6EE6218BECBA839C540B6BC74CE5F
BioGRID  A0931D01FEFF8A4569586D3437E96B70542B818428BAAD7828C37E9D1633520F
```

## Recovered mappings

The 4-layer permutation is recovered by exact multilayer graph matching.  A
1,298-vertex mutually consistent seed is expanded by adjacency fingerprints to
all 5,237 vertices.  The paper's identifiers are nested prefixes, which is
verified on every shared layer for 2→3→4→5 layers.  This recovers all mappings
for 3, 4, and 5 layers.

| Instance | Recovered identifiers | Status |
|---|---:|---|
| 2 layers | 6,044 / 6,052 | eight STRING-indistinguishable IDs remain |
| 3 layers | 5,708 / 5,708 | complete |
| 4 layers | 5,237 / 5,237 | complete |
| 5 layers | 4,223 / 4,223 | complete |

The eight unresolved 2-layer identifiers form five automorphism groups in the
STRING induced graph.  Their candidate ORF names are recorded in
`output/2Layers/unresolved_candidates.tsv`.  A named historical iRefIndex file
is required to select the correct permutation within each group.  No arbitrary
assignment is made.

## Exact validation

For every layer having an independent named source, the reconstructed induced
edge set equals the published edge set after applying the recovered mapping:

| Instance/layer | Compared edges | Missing | Extra |
|---|---:|---:|---:|
| 3L STRING | 718,942 | 0 | 0 |
| 3L WI-PHI | 47,869 | 0 | 0 |
| 4L STRING | 690,049 | 0 | 0 |
| 4L WI-PHI | 45,867 | 0 | 0 |
| 4L BioGRID | 56,899 | 0 | 0 |
| 5L STRING | 524,070 | 0 | 0 |
| 5L WI-PHI | 38,935 | 0 | 0 |
| 5L BioGRID | 51,440 | 0 | 0 |

`output/reconstruction_report.json` contains the machine-readable validation
record.  For 3--5 layers, `named_edges/` contains ORF-name edge lists for all
layers, including the iRefIndex and MINT restricted graphs recovered from the
published numeric instances.  `rebuilt_numeric/` contains canonical sorted
numeric edge lists.  The 2-layer named files are explicitly labelled
`named_edges_partial/`.

## Reproduce

From this directory, run:

```powershell
python recover_4layer_mapping.py
python recover_all_mappings.py
```

The scripts stop on any failed exact edge-set check.  Use
`diagnose_2layer_ambiguity.py` to display the unresolved 2-layer equivalence
classes.

## Download limitation

The Dryad archive `ILPMinPPI_data.zip` is publicly listed with MD5
`c8f7c94e042f3b963502ae00ba481418`, but its file endpoint is currently guarded
by AWS WAF and returned HTTP 403 in this environment.  The historical
iRefIndex server and its old FTP mirror are also unavailable.  Consequently,
the exact original named iRefIndex and MINT downloads have not been claimed as
recovered; their complete *restricted* named layers are nevertheless available
for the fully mapped 3--5 layer instances.
