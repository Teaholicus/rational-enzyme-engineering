#!/usr/bin/env python3
"""
01_acquire_protein.py

Acquire a protein for analysis.

Steps
-----
1. Download the mmCIF structure from the RCSB Protein Data Bank.
2. Determine the UniProt accession for the selected chain (from the mmCIF,
   or from the config if given -- cross-checked against the mmCIF).
3. Download the canonical UniProt FASTA sequence.
4. Verify that structure residue numbering agrees with UniProt numbering.

Step 4 matters: this project uses UniProt numbering throughout, and every
downstream table joins on residue_number. A construct with an N-terminal tag,
a renumbered deposition, or genuine sequence conflicts will silently produce
tables that look fine but describe different residues. Catching it here costs
seconds; catching it at the merge step (08) costs a FoldX run.

Run:
    python workflow/01_acquire_protein.py
    python workflow/01_acquire_protein.py --config config/2JIE_A.yaml

Outputs
-------
data/structures/<PDBID>_<CHAIN>/<PDBID>.cif
data/sequences/<UNIPROTID>(<PDBID>_<CHAIN>).fasta
results/<PDBID>_<CHAIN>/01_acquire/numbering_check.csv   (if mismatches found)
"""

import argparse
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from project_utils import (AA3_TO_1, load_config, protein_tag, resolve_path,
                           results_dir_for)


def download_alphafold_model(uniprot_id: str, structure_dir: Path, pdb_id: str) -> Path:
    """
    Fetch a predicted model from the AlphaFold Database.

    For a target with no experimental structure. Two things change downstream
    and are flagged in the numbering check below:
      - the B-factor column holds pLDDT, not crystallographic B-factor, and the
        direction is INVERTED (high pLDDT = confident; high B-factor = mobile),
      - the model covers the full UniProt sequence including signal peptide and
        disordered termini, which are usually very low confidence.
    """
    structure_dir.mkdir(parents=True, exist_ok=True)
    path = structure_dir / f"{pdb_id}.cif"
    if path.exists():
        print(f"Structure already exists: {path}")
        return path

    # Ask the AFDB API for the current model URL rather than guessing a version
    # number -- the database re-releases models under incrementing versions.
    api = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id.upper()}"
    meta = requests.get(api, timeout=30)
    url = None
    if meta.status_code == 200:
        entries = meta.json()
        if entries:
            url = entries[0].get("cifUrl")
            model_version = entries[0].get("latestVersion", "?")
            print(f"AFDB API reports model version {model_version}")
    if not url:
        # Fallback: probe versions newest-first in case the API is unavailable.
        for v in range(8, 0, -1):
            candidate = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id.upper()}-F1-model_v{v}.cif"
            if requests.head(candidate, timeout=15).status_code == 200:
                url = candidate
                print(f"AFDB API unavailable; found model v{v} by probing")
                break
    if not url:
        raise FileNotFoundError(
            f"No AlphaFold model found for {uniprot_id}. Checked {api} and versions v8-v1.")

    response = requests.get(url, timeout=60)
    response.raise_for_status()
    path.write_bytes(response.content)
    print(f"Downloaded AlphaFold model: {uniprot_id} -> {path}")
    print("  NOTE: B-factor column contains pLDDT (higher = more confident), "
          "which is the OPPOSITE direction to crystallographic B-factor.")
    return path


def download_structure(pdb_id: str, structure_dir: Path) -> Path:
    pdb_id = pdb_id.upper()
    structure_dir.mkdir(parents=True, exist_ok=True)
    cif_file = structure_dir / f"{pdb_id}.cif"

    if cif_file.exists():
        print(f"Structure already exists: {cif_file}")
        return cif_file

    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    cif_file.write_bytes(response.content)
    print(f"Downloaded structure: {pdb_id} -> {cif_file}")
    return cif_file


def get_uniprot_accession(cif_file: Path, chain: str) -> str:
    """Extract the UniProt accession mapped to the selected chain."""
    mmcif = MMCIF2Dict(str(cif_file))
    try:
        asym_ids = mmcif["_struct_ref_seq.pdbx_strand_id"]
        accessions = mmcif["_struct_ref_seq.pdbx_db_accession"]
    except KeyError:
        raise ValueError(
            f"{cif_file.name} has no _struct_ref_seq cross-reference block. "
            f"Set protein.uniprot_id explicitly in the config.")

    if isinstance(asym_ids, str):
        asym_ids, accessions = [asym_ids], [accessions]

    matches = {acc for chains, acc in zip(asym_ids, accessions)
               if chain in [c.strip() for c in chains.split(",")]}
    if not matches:
        raise ValueError(f"No UniProt accession found for chain '{chain}' in {cif_file.name}.")
    if len(matches) > 1:
        raise ValueError(
            f"Chain '{chain}' maps to multiple UniProt accessions {sorted(matches)} "
            f"(chimeric or multi-domain entity). Set protein.uniprot_id explicitly.")
    return matches.pop()


def resolve_uniprot_id(cfg: dict, cif_file: Path, chain: str) -> str:
    """
    Use the config's uniprot_id if provided, but cross-check it against the
    mmCIF mapping. An explicit value is authoritative -- the check exists to
    surface disagreement, not to override the user.
    """
    configured = cfg["protein"].get("uniprot_id")
    try:
        derived = get_uniprot_accession(cif_file, chain)
    except ValueError as exc:
        if configured:
            print(f"Note: could not auto-derive accession ({exc}); using configured {configured}.")
            return configured
        raise

    if configured and configured.upper() != derived.upper():
        print(f"WARNING: config says UniProt {configured}, but {cif_file.name} maps "
              f"chain {chain} to {derived}. Using the configured value -- verify which "
              f"is correct before trusting downstream results.")
        return configured
    if configured:
        print(f"UniProt {configured} (config; matches mmCIF cross-reference)")
        return configured

    print(f"Chain {chain} -> UniProt {derived} (auto-derived from mmCIF)")
    return derived


def download_sequence(uniprot_id: str, pdb_id: str, chain: str, sequence_dir: Path) -> Path:
    sequence_dir.mkdir(parents=True, exist_ok=True)
    fasta_file = sequence_dir / f"{uniprot_id}({pdb_id}_{chain}).fasta"

    if fasta_file.exists():
        print(f"Sequence already exists: {fasta_file}")
        return fasta_file

    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url, timeout=30)
    response.raise_for_status()
    fasta_file.write_bytes(response.content)
    print(f"Downloaded UniProt sequence: {uniprot_id} -> {fasta_file}")
    return fasta_file


def verify_numbering(cif_file: Path, fasta_file: Path, chain: str, out_dir: Path) -> pd.DataFrame:
    """
    Check that residue N in the structure is the same amino acid as position N
    of the UniProt sequence, for every modeled residue.

    Expected outcomes:
      - 100% agreement: numbering conventions match, nothing to do.
      - A handful of mismatches: genuine UniProt sequence conflicts (documented
        in UniProt's own annotations). Downstream, the structure wins.
      - Widespread mismatch: an offset -- an expression tag, a renumbered
        deposition, or the wrong accession. Stop and resolve before continuing.
    """
    seq = str(next(SeqIO.parse(fasta_file, "fasta")).seq)
    structure = MMCIFParser(QUIET=True).get_structure("q", str(cif_file))
    residues = [r for r in structure[0][chain] if r.id[0] == " "]

    rows = []
    for res in residues:
        num = res.id[1]
        struct_aa = AA3_TO_1.get(res.get_resname().strip(), "X")
        seq_aa = seq[num - 1] if 1 <= num <= len(seq) else None
        rows.append({
            "residue_number": num,
            "structure_residue": struct_aa,
            "uniprot_residue": seq_aa,
            "agrees": seq_aa is not None and struct_aa == seq_aa,
        })

    df = pd.DataFrame(rows)
    in_range = df[df.uniprot_residue.notna()]
    n_agree = int(in_range.agrees.sum())
    n_total = len(in_range)
    pct = 100 * n_agree / n_total if n_total else 0.0

    print(f"\nNumbering check: {n_agree}/{n_total} modeled residues agree with "
          f"UniProt numbering ({pct:.1f}%)")
    print(f"  UniProt sequence length: {len(seq)} | modeled residues: {len(df)} "
          f"| range {df.residue_number.min()}-{df.residue_number.max()}")

    out_of_range = df[df.uniprot_residue.isna()]
    if len(out_of_range):
        print(f"  {len(out_of_range)} residue(s) numbered outside the UniProt sequence "
              f"(likely tags or ligand-adjacent entities): "
              f"{out_of_range.residue_number.tolist()[:10]}")

    mismatches = in_range[~in_range.agrees]
    if pct < 90:
        print(f"  *** {len(mismatches)} mismatches. Widespread disagreement suggests a "
              f"numbering OFFSET, not sequence conflicts. Do not proceed until resolved: "
              f"check for an expression tag, a renumbered deposition, or a wrong accession.")
    elif len(mismatches):
        print(f"  {len(mismatches)} isolated mismatch(es) -- consistent with documented "
              f"UniProt sequence conflicts. Downstream steps use the structure's identity "
              f"and flag these as uniprot_structure_conflict.")
        print(mismatches.to_string(index=False))
    else:
        print("  No mismatches.")

    if len(mismatches) or len(out_of_range):
        out_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(out_dir / "numbering_check.csv", index=False)
        print(f"  Full check saved -> {out_dir / 'numbering_check.csv'}")

    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, default=None,
                         help="Override the active target set in config/base.yaml.")
    parser.add_argument("--skip-numbering-check", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    pdb_id = cfg["protein"]["pdb_id"].upper()
    chain = cfg["protein"]["chain"]

    structure_dir = resolve_path(cfg["paths"].get("structures_dir", "data/structures")) / protein_tag(cfg)
    sequence_dir = resolve_path(cfg["paths"].get("sequences_dir", "data/sequences"))

    source = cfg["protein"].get("structure_source", "pdb").lower()
    if source == "alphafold":
        uniprot_id = cfg["protein"].get("uniprot_id")
        if not uniprot_id:
            raise KeyError("structure_source 'alphafold' requires protein.uniprot_id in config.")
        cif_file = download_alphafold_model(uniprot_id, structure_dir, pdb_id)
    elif source == "pdb":
        cif_file = download_structure(pdb_id, structure_dir)
        uniprot_id = resolve_uniprot_id(cfg, cif_file, chain)
    else:
        raise ValueError(f"Unknown protein.structure_source: {source!r} (use 'pdb' or 'alphafold')")
    fasta_file = download_sequence(uniprot_id, pdb_id, chain, sequence_dir)

    if not args.skip_numbering_check:
        verify_numbering(cif_file, fasta_file, chain, results_dir_for(cfg) / "01_acquire")


if __name__ == "__main__":
    main()
