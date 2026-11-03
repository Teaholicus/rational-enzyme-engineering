"""
00_find_binding_site.py

Identify binding-site residues geometrically, from contacts with a ligand in a
co-crystal structure, instead of curating them from the literature.

This is the generalizable path for a novel target: literature-curated
binding-site lists (as used for 1CTN) do not exist for an uncharacterised
protein, but a co-crystal structure gives the same information directly.

Run:
    python workflow/00_find_binding_site.py --config config/2JIE_A.yaml
    python workflow/00_find_binding_site.py --config config/2JIE_A.yaml --list-ligands

Prints a ready-to-paste `functional_residues.binding_site` block for the config.

Notes:
  - Distance is measured between any ligand atom and any residue atom, not
    between centroids -- side-chain contacts are what define a binding site.
  - Waters and common crystallisation additives are skipped by default; use
    --list-ligands first to see what is actually present before trusting the
    automatic pick.
"""
import argparse
from collections import defaultdict

import numpy as np
from Bio.PDB import MMCIFParser

from project_utils import AA3_TO_1, load_config, resolve_path, results_dir_for, structures_dir_for

# Solvent and common crystallisation additives -- present in the file but not
# the biological ligand. Not exhaustive; always check --list-ligands.
SKIP_HETATMS = {
    "HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "MPD", "TRS",
    "ACT", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "IOD", "BR", "NO3",
    "FMT", "ACY", "DMS", "IMD", "BME", "CIT", "MES", "EPE",
}


def list_ligands(structure, chain_id: str):
    """Every non-water heteroatom group in the structure, with atom counts."""
    found = []
    for model in structure:
        for chain in model:
            for res in chain:
                hetflag = res.id[0]
                if hetflag == " ":
                    continue
                name = res.get_resname().strip()
                if name in {"HOH", "WAT", "DOD"}:
                    continue
                found.append({
                    "resname": name,
                    "chain": chain.id,
                    "resnum": res.id[1],
                    "n_atoms": len(list(res)),
                    "likely_solvent": name in SKIP_HETATMS,
                })
    return found


def find_contacts(structure, chain_id: str, ligand_resnames: list, cutoff: float):
    """Residues in `chain_id` with any atom within `cutoff` of any ligand atom."""
    ligand_atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] != " " and res.get_resname().strip() in ligand_resnames:
                    ligand_atoms.extend(list(res))
    if not ligand_atoms:
        raise ValueError(f"No atoms found for ligand(s) {ligand_resnames}. "
                          f"Run with --list-ligands to see what is present.")

    ligand_coords = np.array([a.coord for a in ligand_atoms])

    contacts = {}
    target_chain = structure[0][chain_id]
    for res in target_chain:
        if res.id[0] != " ":
            continue
        coords = np.array([a.coord for a in res])
        # min distance from any residue atom to any ligand atom
        dmin = float(np.min(np.linalg.norm(
            coords[:, None, :] - ligand_coords[None, :, :], axis=-1)))
        if dmin <= cutoff:
            contacts[res.id[1]] = {
                "resnum": res.id[1],
                "resname": res.get_resname().strip(),
                "one_letter": AA3_TO_1.get(res.get_resname().strip(), "X"),
                "min_dist": round(dmin, 2),
            }
    return contacts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, default=None,
                         help="Override the active target set in config/base.yaml.")
    parser.add_argument("--list-ligands", action="store_true",
                         help="List heteroatom groups present and exit.")
    parser.add_argument("--ligand", type=str, default=None,
                         help="Comma-separated ligand residue names (e.g. G2F). "
                              "Defaults to every non-solvent heteroatom group.")
    parser.add_argument("--cutoff", type=float, default=4.5,
                         help="Contact distance in Angstrom (default 4.5).")
    parser.add_argument("--structure", type=str, default=None,
                         help="Use this structure file instead of the configured one -- "
                              "e.g. the target+transferred-ligand file written by "
                              "util_transfer_ligand.py.")
    args = parser.parse_args()

    cfg = load_config(args.config)
    pdb_id = cfg["protein"]["pdb_id"]
    chain_id = cfg["protein"]["chain"]

    if args.structure:
        cif_path = resolve_path(args.structure)
        if not cif_path.exists():
            raise FileNotFoundError(f"{cif_path} not found.")
    else:
        cif_path = structures_dir_for(cfg) / f"{pdb_id}.cif"
        if not cif_path.exists():
            raise FileNotFoundError(f"{cif_path} not found -- run 01_acquire_protein.py first.")

    if cif_path.suffix.lower() == ".cif":
        structure = MMCIFParser(QUIET=True).get_structure("query", str(cif_path))
    else:
        from Bio.PDB import PDBParser
        structure = PDBParser(QUIET=True).get_structure("query", str(cif_path))

    ligands = list_ligands(structure, chain_id)
    if args.list_ligands or not ligands:
        print(f"Heteroatom groups in {pdb_id} (excluding water):\n")
        if not ligands:
            print("  none -- this structure has no ligand; binding-site residues "
                  "must come from another source (co-crystal structure, UniProt "
                  "BINDING annotations, or pocket detection).")
        for lig in ligands:
            flag = "  [likely solvent/additive]" if lig["likely_solvent"] else "  <-- candidate ligand"
            print(f"  {lig['resname']:<6} chain {lig['chain']} "
                  f"resnum {lig['resnum']:<6} {lig['n_atoms']:>3} atoms{flag}")
        return

    if args.ligand:
        ligand_names = [x.strip().upper() for x in args.ligand.split(",")]
    else:
        ligand_names = sorted({l["resname"] for l in ligands if not l["likely_solvent"]})
        if not ligand_names:
            print("Only solvent/additive heteroatoms found. Specify --ligand explicitly "
                  "if one of them is actually the ligand of interest.")
            return
        print(f"Auto-selected ligand(s): {ligand_names}  "
              f"(use --list-ligands to check, --ligand to override)")

    contacts = find_contacts(structure, chain_id, ligand_names, args.cutoff)
    ordered = sorted(contacts.values(), key=lambda c: c["resnum"])

    print(f"\n{len(ordered)} residue(s) within {args.cutoff}A of {ligand_names} "
          f"in chain {chain_id}:\n")
    for c in ordered:
        print(f"  {c['one_letter']}{c['resnum']:<5} {c['resname']:<4} "
              f"min distance {c['min_dist']:>5.2f} A")

    print("\nPaste into the config's functional_residues section:\n")
    print("  binding_site:")
    for c in ordered:
        print(f"    - {c['one_letter']}{c['resnum']}")

    out_dir = results_dir_for(cfg) / "00_binding_site"
    out_dir.mkdir(parents=True, exist_ok=True)
    import pandas as pd
    pd.DataFrame(ordered).to_csv(out_dir / "binding_site_contacts.csv", index=False)
    print(f"\nSaved -> {out_dir / 'binding_site_contacts.csv'}")
    print("\nNote: catalytic residues will appear in this list too -- they contact "
          "the ligand. Keep them in functional_residues.catalytic (they are excluded "
          "from all candidate sets) and remove them from binding_site.")


if __name__ == "__main__":
    main()
