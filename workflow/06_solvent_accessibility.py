"""
06_solvent_accessibility.py

Per-residue solvent accessibility via FreeSASA. Reports both absolute SASA
(confounded with residue size -- Trp is just bigger than Ala) and relative
SASA / RSA (normalized against each residue type's theoretical max exposed
area, using FreeSASA's own built-in reference tables) -- RSA is what's
actually comparable across residue types and what most mutability/stability
literature reports, not raw SASA.

Install:
    pip install freesasa

Run:
    python workflow/06_solvent_accessibility.py
    python workflow/06_solvent_accessibility.py --config config/2JIE_A.yaml

Output columns:
    residue_number, wt_residue          -- join keys, same convention as 03/04/05
    struct_sasa_abs                     -- total SASA (Angstrom^2), size-confounded
    struct_sasa_mainchain / _sidechain  -- SASA split by backbone vs side chain
    struct_rsa                          -- relative SASA (0-1ish, >1 possible for
                                            unusual conformations/chain termini)
    struct_burial_class                 -- categorical: buried / intermediate / exposed
"""
import argparse
from pathlib import Path

import freesasa
import numpy as np
import pandas as pd

from project_utils import load_config, AA3_TO_1, cif_to_pdb, resolve_path, results_dir_for, structures_dir_for


def compute_sasa(pdb_path: Path, chain_id: str) -> pd.DataFrame:
    structure = freesasa.Structure(str(pdb_path))
    result = freesasa.calc(structure)
    residue_areas = result.residueAreas()

    if chain_id not in residue_areas:
        raise KeyError(f"Chain '{chain_id}' not found in FreeSASA results. "
                        f"Chains present: {list(residue_areas.keys())}")

    rows = []
    for resnum_str, area in residue_areas[chain_id].items():
        rows.append({
            "residue_number": int(resnum_str),
            "wt_residue": AA3_TO_1.get(area.residueType, "X"),
            "struct_sasa_abs": area.total,
            "struct_sasa_mainchain": area.mainChain,
            "struct_sasa_sidechain": area.sideChain,
            "struct_rsa": area.relativeTotal if area.hasRelativeAreas else np.nan,
        })

    return pd.DataFrame(rows).sort_values("residue_number").reset_index(drop=True)


def add_burial_class(df: pd.DataFrame, buried_max: float, exposed_min: float) -> pd.DataFrame:
    df = df.copy()
    df["struct_burial_class"] = pd.cut(
        df["struct_rsa"],
        bins=[-0.001, buried_max, exposed_min, np.inf],
        labels=["buried", "intermediate", "exposed"],
    )
    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, default=None,
                         help="Override the active target set in config/base.yaml.")
    args = parser.parse_args()

    cfg = load_config(args.config)
    pdb_id = cfg["protein"]["pdb_id"]
    chain_id = cfg["protein"]["chain"]

    cif_path = structures_dir_for(cfg) / f"{pdb_id}.cif"
    if not cif_path.exists():
        raise FileNotFoundError(f"{cif_path} not found -- run 01_acquire_protein.py first.")

    results_dir = results_dir_for(cfg)
    out_dir = results_dir / "06_solvent_accessibility"
    out_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = out_dir / f"{pdb_id}_{chain_id}.pdb"
    if not pdb_path.exists():
        cif_to_pdb(cif_path, pdb_path, chain_id)
        print(f"Converted {cif_path.name} -> {pdb_path.name}")

    df = compute_sasa(pdb_path, chain_id)

    sasa_cfg = cfg.get("structural_annotation", {})
    df = add_burial_class(df,
                           buried_max=sasa_cfg.get("sasa_buried_max_rsa", 0.05),
                           exposed_min=sasa_cfg.get("sasa_exposed_min_rsa", 0.25))

    df.to_csv(out_dir / "sasa_scores.csv", index=False)
    print(f"Wrote {len(df)} residue rows to sasa_scores.csv")
    print(df["struct_burial_class"].value_counts())


if __name__ == "__main__":
    main()
