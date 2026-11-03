"""
05_bfactor_and_geometry.py

Per-residue structural features computable directly from the crystal
structure with Biopython alone (no external binaries): B-factor, distance to
the nearest catalytic residue, distance to the catalytic-residue centroid,
distance to the nearest known substrate-binding-site residue, and local
packing density (contact number).

Deliberately uses the *original* structure from step 01, not FoldX's
repaired output -- RepairPDB only adjusts side-chain rotamers, not the
backbone (per FoldX's own docs), but it's not guaranteed to preserve the
experimental B-factor column, so B-factor specifically should come from the
real crystallographic data, not a repaired derivative.

Run:
    python workflow/05_bfactor_and_geometry.py
    python workflow/05_bfactor_and_geometry.py --config config/2JIE_A.yaml

Output columns:
    residue_number, wt_residue                    -- join keys, same convention as 03/04
    struct_bfactor_ca                              -- raw B-factor of the CA atom
    struct_bfactor_ca_zscore / _percentile          -- within-protein normalized
                                                        (same rationale as cons_conservation_percentile:
                                                        raw B-factor scale depends on crystal
                                                        resolution/refinement, not comparable
                                                        across structures)
    struct_bfactor_mean_all                         -- mean B-factor over all atoms in the residue
                                                        (captures side-chain flexibility too, not just backbone)
    struct_dist_nearest_catalytic                   -- CA-CA distance (Angstrom) to nearest catalytic residue
    struct_dist_catalytic_centroid                  -- CA distance to the centroid of all catalytic residues
    struct_dist_nearest_binding_site                -- CA-CA distance to nearest known binding-site residue
    struct_contact_density                          -- count of other CAs within a configurable radius
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, NeighborSearch

from project_utils import load_config, AA3_TO_1, parse_residue_code, resolve_path, results_dir_for, structures_dir_for


def load_chain(cif_path: Path, chain_id: str):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("query", str(cif_path))
    return structure[0][chain_id]


def get_standard_residues(chain) -> list:
    """Exclude waters/heteroatoms/ligands -- keep standard amino acid residues only."""
    return [r for r in chain if r.id[0] == " " and "CA" in r]


def resolve_functional_residues(codes: list, residues_by_number: dict, label: str) -> list:
    """
    Parse config-style residue codes ('D313') and sanity-check each against
    the actual structure -- catches numbering-convention mismatches loudly
    instead of silently computing distances to the wrong residue.
    """
    resolved = []
    for code in codes:
        expected_aa, resnum = parse_residue_code(code)
        residue = residues_by_number.get(resnum)
        if residue is None:
            print(f"Warning: {label} residue {code} not found in structure (no residue at position {resnum}).")
            continue
        actual_aa = AA3_TO_1.get(residue.get_resname(), "?")
        if actual_aa != expected_aa:
            print(f"Warning: {label} residue {code} expected {expected_aa} but structure has "
                  f"{actual_aa} at position {resnum} -- check numbering convention before trusting results.")
        resolved.append(residue)
    return resolved


def find_disulfides(residues: list, max_sg_distance: float = 2.5) -> set:
    """
    Cysteines whose SG atoms are within bonding distance of another SG.

    Flagged for the same reason as glycines: mutating one breaks a covalent
    cross-link, so FoldX will report every substitution as heavily destabilising
    and the position drops out with no viable candidate. Better to label it than
    to leave an unexplained gap.
    """
    cys_sg = [(r.id[1], r["SG"]) for r in residues
              if r.get_resname().strip() == "CYS" and "SG" in r]
    bonded = set()
    for i, (num_i, sg_i) in enumerate(cys_sg):
        for num_j, sg_j in cys_sg[i + 1:]:
            if np.linalg.norm(sg_i.coord - sg_j.coord) <= max_sg_distance:
                bonded.add(num_i)
                bonded.add(num_j)
    return bonded


def compute_bfactor(residue) -> tuple[float, float]:
    """Returns (CA B-factor, mean B-factor over all atoms in the residue)."""
    ca_bfactor = residue["CA"].get_bfactor()
    all_bfactors = [atom.get_bfactor() for atom in residue]
    return ca_bfactor, float(np.mean(all_bfactors))


def annotate_predicted_confidence(df: pd.DataFrame, plddt_min: float) -> pd.DataFrame:
    """
    For AlphaFold models the B-factor column holds pLDDT, not a crystallographic
    B-factor, and the direction is INVERTED: high pLDDT means well-predicted,
    high B-factor means mobile. Relabel it and flag unreliable positions.

    AlphaFold's own confidence bands: >90 very high, 70-90 confident,
    50-70 low, <50 very low (and in that range often disordered). A predicted
    model also covers the full UniProt sequence, so signal peptides and
    disordered termini are present but meaningless -- their geometry, SASA and
    ddG values are all unreliable, and without this flag they would enter the
    candidate ranking indistinguishable from well-resolved positions.
    """
    df = df.copy()
    df["struct_plddt"] = df["struct_bfactor_ca"]
    df["struct_confidence_ok"] = df["struct_plddt"] >= plddt_min
    df["struct_plddt_band"] = pd.cut(
        df["struct_plddt"], bins=[0, 50, 70, 90, 100],
        labels=["very_low", "low", "confident", "very_high"])

    n_fail = int((~df["struct_confidence_ok"]).sum())
    print(f"pLDDT: median {df.struct_plddt.median():.1f}, "
          f"{n_fail} of {len(df)} residues below {plddt_min}")
    print(df["struct_plddt_band"].value_counts().sort_index().to_string())
    if n_fail:
        failed = df.loc[~df["struct_confidence_ok"], "residue_number"].tolist()
        print(f"  low-confidence positions: {failed[:20]}"
              f"{'...' if len(failed) > 20 else ''}")
    # B-factor framing does not apply to a predicted model
    df = df.drop(columns=[c for c in df.columns
                           if c.startswith("struct_bfactor_ca_")], errors="ignore")
    return df


def add_relative_normalization(df: pd.DataFrame, raw_col: str, prefix: str) -> pd.DataFrame:
    """Same within-protein normalization rationale as conservation's percentile/grade:
    raw B-factor scale depends on crystal resolution/refinement and isn't comparable
    across structures, so rank each position relative to others in this same protein."""
    df = df.copy()
    df[f"{prefix}_zscore"] = (df[raw_col] - df[raw_col].mean()) / df[raw_col].std()
    df[f"{prefix}_percentile"] = 100 * df[raw_col].rank(pct=True)
    return df


def compute_geometry(residues: list, catalytic: list, binding_site: list,
                      contact_radius: float, is_predicted: bool = False,
                      plddt_min: float = 70.0) -> pd.DataFrame:
    disulfide_positions = find_disulfides(residues)
    if disulfide_positions:
        print(f"Disulfide-bonded cysteines: {sorted(disulfide_positions)} "
              f"({len(disulfide_positions)//2} bond(s))")
    ca_coords = {r.id[1]: r["CA"].coord for r in residues}
    catalytic_coords = [r["CA"].coord for r in catalytic]
    binding_coords = [r["CA"].coord for r in binding_site]
    centroid = np.mean(catalytic_coords, axis=0) if catalytic_coords else None

    all_ca_atoms = [r["CA"] for r in residues]
    neighbor_search = NeighborSearch(all_ca_atoms)

    rows = []
    for r in residues:
        resnum = r.id[1]
        ca = ca_coords[resnum]
        wt = AA3_TO_1.get(r.get_resname(), "X")

        dist_nearest_cat = (min(np.linalg.norm(ca - c) for c in catalytic_coords)
                             if catalytic_coords else np.nan)
        dist_centroid = float(np.linalg.norm(ca - centroid)) if centroid is not None else np.nan
        dist_nearest_bind = (min(np.linalg.norm(ca - c) for c in binding_coords)
                              if binding_coords else np.nan)

        neighbors = neighbor_search.search(ca, contact_radius, level="A")
        contact_density = len(neighbors) - 1  # exclude self

        ca_bf, mean_bf = compute_bfactor(r)

        rows.append({
            "residue_number": resnum,
            "wt_residue": wt,
            "struct_bfactor_ca": ca_bf,
            "struct_bfactor_mean_all": mean_bf,
            "struct_dist_nearest_catalytic": dist_nearest_cat,
            "struct_dist_catalytic_centroid": dist_centroid,
            "struct_dist_nearest_binding_site": dist_nearest_bind,
            "struct_contact_density": contact_density,
            "struct_in_disulfide": resnum in disulfide_positions,
        })

    df = pd.DataFrame(rows).sort_values("residue_number").reset_index(drop=True)
    df = add_relative_normalization(df, "struct_bfactor_ca", "struct_bfactor_ca")
    if is_predicted:
        df = annotate_predicted_confidence(df, plddt_min)
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

    chain = load_chain(cif_path, chain_id)
    residues = get_standard_residues(chain)
    residues_by_number = {r.id[1]: r for r in residues}
    print(f"Loaded {len(residues)} standard residues from {cif_path.name} chain {chain_id}.")

    func_cfg = cfg.get("functional_residues", {})
    catalytic = resolve_functional_residues(func_cfg.get("catalytic", []), residues_by_number, "catalytic")
    binding_site = resolve_functional_residues(func_cfg.get("binding_site", []), residues_by_number, "binding_site")
    print(f"Resolved {len(catalytic)} catalytic and {len(binding_site)} binding-site reference residues.")

    struct_cfg = cfg.get("structural_annotation", {})
    contact_radius = struct_cfg.get("contact_radius_angstrom", 10.0)
    is_predicted = cfg["protein"].get("structure_source", "pdb").lower() == "alphafold"
    if is_predicted:
        print("Predicted model: B-factor column treated as pLDDT (high = confident).")
    df = compute_geometry(residues, catalytic, binding_site, contact_radius,
                           is_predicted=is_predicted,
                           plddt_min=struct_cfg.get("plddt_min", 70.0))

    results_dir = results_dir_for(cfg)
    out_dir = results_dir / "05_bfactor_geometry"
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "structural_scores.csv", index=False)
    print(f"Wrote {len(df)} residue rows to structural_scores.csv")


if __name__ == "__main__":
    main()
