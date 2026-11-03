"""
util_transfer_ligand.py

Transfer a ligand from a homologous co-crystal structure onto the target
structure by structural superimposition.

Why this rather than predicting the complex directly: a transferred ligand pose
is an *observed* position from an experiment, carried over by a rigid-body
transform whose quality is reported as an RMSD the reader can judge. A predicted
pose is an inference with no comparable audit trail, and co-folding methods place
roughly three-quarters of ligands correctly -- without indicating which quarter
is wrong. Use prediction only when no homologous co-crystal exists.

The transferred ligand is used *only* as a geometric probe, to identify which
target residues line the binding site. It contributes no coordinates to the
target protein: every pipeline step (distances, FoldX, SASA, DSSP) runs on the
original target structure, so no predicted ligand-induced conformational change
enters the analysis.

Method: align the two sequences, collect CA atoms of aligned identical-position
pairs, solve for the rotation/translation minimising their RMSD (Kabsch), then
apply that transform to the reference ligand.

Run:
    # see what ligands the reference contains
    python workflow/util_transfer_ligand.py --reference 1WD4 --list-ligands

    # transfer, picking the ligand copy nearest the catalytic residues
    python workflow/util_transfer_ligand.py --reference 1WD4 --ligand AHR

Then derive the binding site from the written structure:
    python workflow/util_find_binding_site.py --structure <written path>
"""
import argparse
from pathlib import Path

import numpy as np
import requests
from Bio import Align
from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Superimposer
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain

from project_utils import (AA3_TO_1, load_config, parse_residue_code,
                            resolve_path, results_dir_for, structures_dir_for)

SOLVENT = {"HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "MPD",
           "ACT", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "NAG", "BMA", "MAN"}


def fetch_structure(pdb_id: str, dest_dir: Path) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    path = dest_dir / f"{pdb_id.upper()}.cif"
    if path.exists():
        print(f"Reference structure already present: {path}")
        return path
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    path.write_bytes(r.content)
    print(f"Downloaded reference {pdb_id.upper()} -> {path}")
    return path


def load_chain(path: Path, chain_id: str = None):
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() == ".cif" else PDBParser(QUIET=True)
    structure = parser.get_structure(path.stem, str(path))
    model = structure[0]
    if chain_id:
        return structure, model[chain_id]
    chains = [c for c in model if any(r.id[0] == " " for r in c)]
    if not chains:
        raise ValueError(f"No protein chain found in {path}")
    return structure, chains[0]


def chain_sequence(chain):
    """Returns (one-letter sequence, list of residues) for standard residues with CA."""
    residues = [r for r in chain if r.id[0] == " " and "CA" in r]
    seq = "".join(AA3_TO_1.get(r.get_resname().strip(), "X") for r in residues)
    return seq, residues


def aligned_ca_pairs(target_chain, ref_chain, min_identity_warn=0.30):
    """Align sequences; return matched CA atom pairs at identical positions."""
    t_seq, t_res = chain_sequence(target_chain)
    r_seq, r_res = chain_sequence(ref_chain)

    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")
    alignment = aligner.align(t_seq, r_seq)[0]

    fixed, moving, pairs, n_aligned, n_identical = [], [], [], 0, 0
    for (t_start, t_end), (r_start, r_end) in zip(*alignment.aligned):
        for offset in range(t_end - t_start):
            ti, ri = t_start + offset, r_start + offset
            n_aligned += 1
            if t_seq[ti] == r_seq[ri]:
                n_identical += 1
                fixed.append(t_res[ti]["CA"])
                moving.append(r_res[ri]["CA"])
                pairs.append((t_res[ti], r_res[ri]))

    identity = n_identical / min(len(t_seq), len(r_seq)) if min(len(t_seq), len(r_seq)) else 0
    print(f"Sequence alignment: target {len(t_seq)} aa, reference {len(r_seq)} aa")
    print(f"  aligned columns: {n_aligned} | identical: {n_identical} "
          f"({identity:.1%} of the shorter sequence)")
    if identity < min_identity_warn:
        print(f"  *** WARNING: identity below {min_identity_warn:.0%}. A superimposition "
              f"at this level is unreliable; check the RMSD carefully before trusting "
              f"the transferred ligand.")
    return fixed, moving, pairs, identity


def refit_near_ligand(pairs, ligand_residue, radius: float = 12.0):
    """
    Re-superimpose using only residue pairs near the ligand.

    A global fit across a multi-domain protein is the wrong tool for ligand
    transfer: GH54 and CBM42 are joined by a linker, so any difference in
    inter-domain angle between model and crystal forces the fit to compromise
    between the two domains and fit neither well. What matters for transferring
    a ligand is local agreement around the binding site, so the fit is repeated
    using only the residues that actually surround it.

    Returns (Superimposer, n_pairs_used) or (None, 0) if too few pairs.
    """
    lig_coords = np.array([a.coord for a in ligand_residue])
    local_fixed, local_moving = [], []
    for t_res, r_res in pairs:
        d = float(np.min(np.linalg.norm(
            r_res["CA"].coord[None, :] - lig_coords, axis=-1)))
        if d <= radius:
            local_fixed.append(t_res["CA"])
            local_moving.append(r_res["CA"])

    if len(local_fixed) < 10:
        return None, len(local_fixed)

    sup = Superimposer()
    sup.set_atoms(local_fixed, local_moving)
    return sup, len(local_fixed)


def map_positions_to_reference(pairs, target_positions):
    """Translate target residue numbers into their reference equivalents."""
    lookup = {t_res.id[1]: r_res for t_res, r_res in pairs}
    return [lookup[p] for p in target_positions if p in lookup]


def list_ligands(structure):
    out = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.id[0] == " ":
                    continue
                name = res.get_resname().strip()
                if name in {"HOH", "WAT", "DOD"}:
                    continue
                out.append({"resname": name, "chain": chain.id, "resnum": res.id[1],
                             "n_atoms": len(list(res)), "residue": res,
                             "likely_solvent": name in SOLVENT})
    return out


def pick_ligand_copy(candidates, reference_catalytic_residues):
    """
    A reference structure may hold several copies of the same ligand in
    different sites -- 1WD4 has three arabinose molecules: one in the catalytic
    pocket and two on the CBM42 targeting module. Pick the copy closest to the
    catalytic residues, or the binding site merges unrelated pockets.

    Distances are measured in the reference structure's own frame, between the
    ligand and the reference residues equivalent to the target's catalytic
    residues, so this does not depend on the superimposition being good.
    """
    if not reference_catalytic_residues:
        print("  No catalytic residues available to choose by. Specify the copy "
              "explicitly with --ligand-resnum, or fill in functional_residues.catalytic.")
        return None

    cat_coords = np.array([r["CA"].coord for r in reference_catalytic_residues
                            if "CA" in r])
    if not len(cat_coords):
        return None
    scored = []
    for c in candidates:
        lig_coords = np.array([a.coord for a in c["residue"]])
        d = float(np.min(np.linalg.norm(
            lig_coords[:, None, :] - cat_coords[None, :, :], axis=-1)))
        scored.append((d, c))
    scored.sort(key=lambda x: x[0])

    print(f"  {len(candidates)} copies of this ligand; distance to nearest catalytic residue:")
    for d, c in scored:
        print(f"    {c['resname']} {c['chain']}{c['resnum']}: {d:.1f} A"
              f"{'   <- selected' if c is scored[0][1] else ''}")
    return scored[0][1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, default=None)
    parser.add_argument("--reference", type=str, required=True,
                         help="PDB ID of the homologous co-crystal (e.g. 1WD4).")
    parser.add_argument("--reference-chain", type=str, default=None)
    parser.add_argument("--ligand", type=str, default=None,
                         help="Ligand residue name to transfer (e.g. AHR).")
    parser.add_argument("--ligand-resnum", type=int, default=None,
                         help="Residue number of the specific ligand copy to transfer. "
                              "Use when a structure holds several copies of the same "
                              "ligand (1WD4 has three AHR: one catalytic, two on CBM42).")
    parser.add_argument("--local-radius", type=float, default=12.0,
                         help="Radius around the ligand used for the local refit "
                              "(default 12 A).")
    parser.add_argument("--list-ligands", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    pdb_id = cfg["protein"]["pdb_id"]
    chain_id = cfg["protein"]["chain"]
    catalytic_positions = [parse_residue_code(c)[1]
                            for c in cfg.get("functional_residues", {}).get("catalytic", [])]

    target_path = structures_dir_for(cfg) / f"{pdb_id}.cif"
    if not target_path.exists():
        raise FileNotFoundError(f"{target_path} not found -- run 01_acquire_protein.py first.")

    out_dir = results_dir_for(cfg) / "00_binding_site"
    out_dir.mkdir(parents=True, exist_ok=True)

    ref_path = fetch_structure(args.reference, resolve_path("data/structures") / "_references")
    ref_structure, ref_chain = load_chain(ref_path, args.reference_chain)
    target_structure, target_chain = load_chain(target_path, chain_id)

    ligands = list_ligands(ref_structure)
    if args.list_ligands or not args.ligand:
        print(f"\nHeteroatom groups in {args.reference.upper()} (excluding water):\n")
        for lig in ligands:
            flag = "  [likely solvent/glycan/additive]" if lig["likely_solvent"] else "  <-- candidate"
            print(f"  {lig['resname']:<6} chain {lig['chain']} resnum {lig['resnum']:<6} "
                  f"{lig['n_atoms']:>3} atoms{flag}")
        if not args.ligand:
            print("\nRe-run with --ligand <CODE> to transfer one.")
        return

    print()
    fixed, moving, pairs, identity = aligned_ca_pairs(target_chain, ref_chain)
    if len(fixed) < 20:
        raise ValueError(f"Only {len(fixed)} equivalent CA pairs -- too few to superimpose.")

    global_sup = Superimposer()
    global_sup.set_atoms(fixed, moving)
    print(f"\nGlobal superimposition on {len(fixed)} CA pairs: RMSD = {global_sup.rms:.2f} A")
    if global_sup.rms > 2.0:
        print("  A high global RMSD at high sequence identity usually means the domains")
        print("  sit at different relative angles, not that the folds differ. The local")
        print("  refit below is what matters for transferring the ligand.")

    # choose which copy of the ligand to transfer
    wanted = [l for l in ligands if l["resname"] == args.ligand.upper()]
    if not wanted:
        raise ValueError(f"Ligand {args.ligand} not found in {args.reference.upper()}. "
                          f"Use --list-ligands to see what is present.")

    if args.ligand_resnum is not None:
        match = [l for l in wanted if l["resnum"] == args.ligand_resnum]
        if not match:
            raise ValueError(f"No {args.ligand.upper()} copy with resnum "
                              f"{args.ligand_resnum}. Present: "
                              f"{[l['resnum'] for l in wanted]}")
        chosen = match[0]
        print(f"\nUsing {chosen['resname']} {chosen['chain']}{chosen['resnum']} (specified).")
    elif len(wanted) == 1:
        chosen = wanted[0]
        print(f"\nSingle copy: {chosen['resname']} {chosen['chain']}{chosen['resnum']}")
    else:
        print(f"\n{len(wanted)} copies of {args.ligand.upper()} present. Choosing by "
              f"distance to the catalytic residues:")
        ref_catalytic = map_positions_to_reference(pairs, catalytic_positions)
        if catalytic_positions and not ref_catalytic:
            print("  (catalytic residues did not map onto the reference)")
        chosen = pick_ligand_copy(wanted, ref_catalytic)
        if chosen is None:
            print(f"\n  Copies available: "
                  f"{[(l['chain'], l['resnum']) for l in wanted]}")
            raise SystemExit("Cannot choose a copy. Re-run with --ligand-resnum.")

    # refit locally around the chosen ligand, then transfer using that transform
    local_sup, n_local = refit_near_ligand(pairs, chosen["residue"], args.local_radius)
    if local_sup is None:
        print(f"\nOnly {n_local} residue pairs within {args.local_radius} A of the ligand "
              f"-- too few for a local refit; falling back to the global transform.")
        sup, fit_label, n_used = global_sup, "global", len(fixed)
    else:
        print(f"\nLocal refit on {n_local} CA pairs within {args.local_radius} A of the "
              f"ligand: RMSD = {local_sup.rms:.2f} A")
        if local_sup.rms > 2.0:
            print("  *** WARNING: the binding site itself does not superimpose well. "
                  "Do not trust the transferred ligand.")
        elif local_sup.rms > 1.0:
            print("  Moderate. Usable, but quote this number alongside any conclusion.")
        else:
            print("  Good local agreement.")
        sup, fit_label, n_used = local_sup, "local", n_local

    ligand_residue = chosen["residue"].copy()
    sup.apply(list(ligand_residue))

    # write target protein + transferred ligand as one file
    combined = Structure("combined")
    model = Model(0)
    combined.add(model)
    prot = Chain(chain_id)
    model.add(prot)
    for res in target_chain:
        if res.id[0] == " ":
            prot.add(res.copy())
    lig_chain = Chain("L")
    model.add(lig_chain)
    ligand_residue.id = (f"H_{chosen['resname']}", 1, " ")
    lig_chain.add(ligand_residue)

    out_path = out_dir / f"{pdb_id}_{chain_id}_with_{args.ligand.upper()}.pdb"
    io = PDBIO()
    io.set_structure(combined)
    io.save(str(out_path))
    content = out_path.read_text()
    if not content.startswith("HEADER"):
        out_path.write_text(f"HEADER    {out_path.stem[:40]:<40}01-JAN-00   {pdb_id[:4].upper():>4}\n" + content)

    print(f"\nWrote target + transferred ligand -> {out_path}")
    (out_dir / "transfer_provenance.txt").write_text(
        f"target={pdb_id}_{chain_id}\n"
        f"reference={args.reference.upper()}\n"
        f"ligand={args.ligand.upper()} copy chain {chosen['chain']} resnum {chosen['resnum']}\n"
        f"sequence_identity={identity:.3f}\n"
        f"global_rmsd={global_sup.rms:.3f} (on {len(fixed)} CA pairs)\n"
        f"transform_used={fit_label}\n"
        f"applied_rmsd={sup.rms:.3f} (on {n_used} CA pairs"
        f"{f', within {args.local_radius} A of the ligand' if fit_label == 'local' else ''})\n")
    print(f"Provenance -> {out_dir / 'transfer_provenance.txt'}")
    print(f"\nNext: python workflow/util_find_binding_site.py --structure {out_path}")


if __name__ == "__main__":
    main()
