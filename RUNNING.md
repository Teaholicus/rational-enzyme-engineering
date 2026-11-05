# Running the pipeline

## Environment

Everything runs in WSL/Ubuntu **except the notebooks**, which run under Windows Jupyter reading the same files (`N:\...` from Windows is `/mnt/n/...` from WSL — one drive, two paths).

| Step | Needs | Reason |
|---|---|---|
| `01` | either | pure Python |
| `02` | **WSL** | `mmseqs` |
| `03` | **WSL** | `mafft` |
| `04` | either | `foldx5.exe` is a Windows binary; the script detects WSL and converts paths with `wslpath`, so WSL interop launches it fine |
| `05` `06` | either | Biopython / FreeSASA, both pip-installed |
| `07` | **WSL** | `mkdssp` |
| `08`, `validation/` | **Windows Jupyter** | no external binaries needed |
| `util_*` | either | pure Python |

Simplest habit: run every numbered script from WSL, open notebooks in Windows Jupyter.

```bash
cd /mnt/n/Jupiter_Projects/rational-enzyme-engineering
which mmseqs mafft mkdssp      # confirm the WSL tools are on PATH
```

---

## One-time setup

```bash
bash tools/setup_uniref90.sh /absolute/path/to/db_dir
#   ~20 GB download, 1-3 h indexing. Point it at a fast drive.
#   Skips createindex by default -- not needed for occasional searches.
```

FoldX needs an academic licence (register, download, place in `tools/foldx5/`). Not redistributable, so it is not in this repo.

---

## Selecting a target

Edit one line in `config/active.yaml`:

```yaml
active_protein: "2JIE_A"       # names config/2JIE_A.yaml
```

This is deliberately a separate file from `base.yaml`: `base.yaml` holds shared
pipeline parameters, `active.yaml` holds which target you are currently working
on. Updating one never overwrites the other.

Every script and notebook picks this up and prints `[config] target: ...` on startup. To run a one-off against a different target without editing the file, pass `--config config/<name>.yaml` to any script.

Results are namespaced under `results/<PDBID>_<CHAIN>/`, so targets never overwrite each other.

---

## Main sequence

```bash
python workflow/01_acquire_protein.py
#   Downloads structure (RCSB, or AlphaFold DB if structure_source: alphafold)
#   and the UniProt sequence; verifies structure numbering against UniProt.
#   --skip-numbering-check     skip the verification
#
#   READ THE NUMBERING CHECK. 100% = clean. A handful of mismatches = documented
#   UniProt sequence conflicts (fine, handled downstream). Below 90% = an offset;
#   stop and resolve it before spending time on anything else.

python workflow/02_search_homologs.py
#   MMseqs2 vs UniRef90, filter by identity/coverage, cluster at 90%.
#   Watch the printed counts. If raw hits hit max_seqs, the set is truncated by
#   score and biased toward close homologs -- raise max_seqs and re-run.

python workflow/03_conservation_analysis.py
#   MAFFT alignment, Henikoff-weighted Shannon entropy, per-AA frequencies
#   (weighted and unweighted).

python workflow/04_foldx_ddg.py
#   RepairPDB + PositionScan: ddG for all 20 substitutions at every position.
#   SLOW -- 12 h for a 450-residue protein. Start it when you can leave it.
#   Skips RepairPDB and the scan if their outputs already exist.
#   --target-positions-override 73,226,351   scan only these positions
#   --output-suffix _patch                   don't overwrite the main results
#   --artifacts delete|archive|keep          per-mutation file cleanup (default delete)

python workflow/05_bfactor_and_geometry.py
#   B-factor (or pLDDT for predicted models), distances to catalytic and
#   binding-site residues, contact density, disulfide detection.

python workflow/06_solvent_accessibility.py
#   FreeSASA: absolute SASA, relative SASA, burial class. Run on Win

python workflow/07_secondary_structure.py
#   DSSP: 3-state secondary structure, phi/psi, glycine positive-phi flags.
#   --dssp-exe /path/to/mkdssp    if Biopython cannot find the binary

# then, in Windows Jupyter:
#   workflow/08_merge_and_analysis.ipynb
```

`05`, `06` and `07` depend only on the structure, so they can run while `04` is still going.

---

## Preparing a config for a new target

These run once, before `02`, and write the values you paste into the target's config file. They are utilities rather than pipeline steps — they prepare configuration, they do not produce pipeline output.

**When the target's own structure has a ligand:**

```bash
python workflow/util_find_binding_site.py --list-ligands
python workflow/util_find_binding_site.py --ligand <CODE>
#   --cutoff 4.5        contact distance (default 4.5 A)
#   --structure <path>  use a different structure file
```

**When it does not, but a homologous co-crystal exists:**

```bash
python workflow/util_transfer_ligand.py --reference 1WD4 --list-ligands
python workflow/util_transfer_ligand.py --reference 1WD4 --ligand AHR --ligand-resnum 601
#   --ligand-resnum     which copy, when a structure holds several
#   --local-radius 12   radius for the local refit around the ligand
#   --reference-chain   if the reference has several chains

python workflow/util_find_binding_site.py --structure results/<TARGET>/00_binding_site/<file>.pdb
```

Read both RMSDs that `util_transfer_ligand` prints. A high *global* RMSD at high sequence identity usually means the domains sit at different relative angles, not that the folds differ — the *local* RMSD around the ligand is what determines whether the transfer is trustworthy. Under ~1 Å is good.

Fill in `functional_residues.catalytic` before running the transfer: it uses the catalytic residues to pick the right ligand copy, and without them falls back to the first one.

---

## Maintenance

```bash
python workflow/util_clean_foldx_artifacts.py            # dry run, lists only
python workflow/util_clean_foldx_artifacts.py --execute
#   --archive-to results/foldx_artifacts    move instead of delete
```

A full scan leaves one PDB and one energies file per mutation (~9,000 files for a 450-residue protein) in the FoldX directory. `04` now cleans up after itself; this clears older leftovers.

---

## Validation notebooks

```
validation/bglb_kinetic_validation.ipynb   # needs active_protein: "2JIE_A"
validation/consurf_validation.ipynb        # needs a ConSurf run under results/<TARGET>/consurf/
```

Both assert on the active target and fail loudly if it is wrong.

---

## Notes on predicted structures

With `structure_source: alphafold` in the target config:

- The B-factor column holds **pLDDT**, where high means confident — the opposite direction to crystallographic B-factor, where high means mobile. `05` detects this, writes `struct_plddt` and `struct_plddt_band`, and drops the B-factor percentile columns, which would otherwise be read backwards.
- The model spans the whole UniProt sequence including signal peptide and disordered termini. Positions below `plddt_min` (default 70, AlphaFold's own "confident" boundary) are flagged `struct_confidence_ok = False` and excluded from candidate sets in `08`.
- Set `plddt_min` per target in its config if the default is too strict or too loose.
