"""
07_secondary_structure.py

Per-residue secondary structure via DSSP, which assigns structure from
backbone hydrogen-bond geometry (Kabsch & Sander 1983), not just phi/psi
lookup. Also pulls phi/psi directly and flags glycines in a positive-phi
backbone conformation -- a region of conformational space essentially
inaccessible to any other amino acid, so mutating such a position away from
Gly is a known route to severe backbone strain.

Does NOT re-extract DSSP's own SASA -- 06_solvent_accessibility.py already
provides a cleaner FreeSASA-based RSA; reporting two different SASA numbers
for the same residue would create ambiguity, not signal.

Install (Ubuntu/WSL):
    sudo apt-get install dssp
    # or: conda install -c bioconda dssp
    which mkdssp   # confirm it's on PATH before running this script

Run:
    python workflow/07_secondary_structure.py
    python workflow/07_secondary_structure.py --config config/2JIE_A.yaml

Output columns:
    residue_number, wt_residue          -- join keys, same convention as 03/04/05/06
    struct_ss_raw                       -- DSSP's native 8-state code (H,G,I,E,B,T,S,-)
    struct_ss_class                     -- simplified: helix / sheet / loop
    struct_phi, struct_psi              -- backbone dihedral angles (degrees)
    struct_gly_positive_phi             -- True if WT is Gly AND phi > 0
                                            (backbone conformation other residues can't access)
"""
import argparse
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

from project_utils import load_config, cif_to_pdb, resolve_path, results_dir_for, structures_dir_for

SS_SIMPLIFY = {
    "H": "helix", "G": "helix", "I": "helix",
    "E": "sheet", "B": "sheet",
    "T": "loop", "S": "loop", "-": "loop", "P": "loop",
}


def run_dssp(pdb_path: Path, chain_id: str, dssp_exe: str = "mkdssp") -> pd.DataFrame:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("query", str(pdb_path))
    model = structure[0]
    dssp = DSSP(model, str(pdb_path), dssp=dssp_exe, file_type="PDB")

    rows = []
    for key in dssp.keys():
        chain, res_id = key
        if chain != chain_id:
            continue
        _idx, aa, ss_raw, _rel_asa, phi, psi, *_ = dssp[key]

        rows.append({
            "residue_number": res_id[1],
            "wt_residue": aa,
            "struct_ss_raw": ss_raw,
            "struct_ss_class": SS_SIMPLIFY.get(ss_raw, "loop"),
            "struct_phi": phi,
            "struct_psi": psi,
        })

    df = pd.DataFrame(rows).sort_values("residue_number").reset_index(drop=True)
    # DSSP uses 360.0 as a sentinel for "undefined" (chain termini) -- convert to NaN
    for col in ("struct_phi", "struct_psi"):
        df.loc[df[col] == 360.0, col] = float("nan")

    df["struct_gly_positive_phi"] = (df["wt_residue"] == "G") & (df["struct_phi"] > 0)
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, default=None,
                         help="Override the active target set in config/base.yaml.")
    parser.add_argument("--dssp-exe", type=str, default="mkdssp")
    args = parser.parse_args()

    cfg = load_config(args.config)
    pdb_id = cfg["protein"]["pdb_id"]
    chain_id = cfg["protein"]["chain"]

    cif_path = structures_dir_for(cfg) / f"{pdb_id}.cif"
    if not cif_path.exists():
        raise FileNotFoundError(f"{cif_path} not found -- run 01_acquire_protein.py first.")

    results_dir = results_dir_for(cfg)
    out_dir = results_dir / "07_secondary_structure"
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = out_dir / f"{pdb_id}_{chain_id}.pdb"
    if not pdb_path.exists():
        cif_to_pdb(cif_path, pdb_path, chain_id)
        print(f"Converted {cif_path.name} -> {pdb_path.name}")

    df = run_dssp(pdb_path, chain_id, dssp_exe=args.dssp_exe)
    df.to_csv(out_dir / "secondary_structure.csv", index=False)
    print(f"Wrote {len(df)} residue rows to secondary_structure.csv")
    print(df["struct_ss_class"].value_counts())
    n_flagged = df["struct_gly_positive_phi"].sum()
    if n_flagged:
        print(f"{n_flagged} glycine(s) flagged in positive-phi conformation "
              f"(likely structurally required, risky to mutate away).")


if __name__ == "__main__":
    main()
