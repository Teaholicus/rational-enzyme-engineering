"""
03_conservation_analysis.py

Align query + homologs with MAFFT, then compute a per-residue conservation
score mapped back onto the reference (query) numbering, using
Henikoff & Henikoff (1994) position-based sequence weighting so that
over-represented near-identical sequences don't flatten the signal at
genuinely variable positions.

Run as a script (from WSL, needs mafft on PATH):
    python workflow/03_conservation_analysis.py
    python workflow/03_conservation_analysis.py --config config/2JIE_A.yaml

Target defaults to `active_protein` in config/base.yaml.

Can also be imported from a notebook:
    from conservation_analysis import run_from_config
    df = run_from_config()

Output columns:
    residue_number, wt_residue   -- join keys, identical across all pipeline
                                     tables (02/03/04/...), unprefixed
    cons_*                       -- everything from this step, prefixed so it
                                     can't collide with foldx_*/struct_* later
"""
import argparse
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO

from project_utils import load_config, results_dir_for

GAP_CHARS = {"-", "."}


def run_mafft(fasta_in: Path, aln_out: Path, threads: int = 8) -> Path:
    with open(aln_out, "w") as f:
        subprocess.run(["mafft", "--auto", "--thread", str(threads), str(fasta_in)],
                        stdout=f, check=True)
    return aln_out


def _sequence_weights(alignment) -> np.ndarray:
    """
    Henikoff & Henikoff (1994) position-based weighting: at each alignment
    column, residues that are rare among the sequences present get more
    weight, so a clade of 50 near-identical homologs doesn't count as 50
    independent votes for "conserved."
    """
    n_seq = len(alignment)
    n_col = alignment.get_alignment_length()
    weights = np.zeros(n_seq)

    for col in range(n_col):
        column = alignment[:, col]
        residues_here = [r for r in column if r not in GAP_CHARS]
        if not residues_here:
            continue
        counts = pd.Series(residues_here).value_counts()
        n_distinct = len(counts)
        for i, r in enumerate(column):
            if r not in GAP_CHARS:
                weights[i] += 1.0 / (n_distinct * counts[r])

    weights /= weights.sum()
    return weights


ALL_AAS = list("ACDEFGHIKLMNPQRSTVWY")


def _weighted_residue_distribution(column: str, weights: np.ndarray) -> tuple[dict, dict, float]:
    """
    Returns (weighted frequencies, unweighted frequencies, weighted Shannon entropy).

    Both frequency sets are reported because they answer different questions:
      - weighted:   corrects for redundant sampling (a clade of 50 near-identical
                    homologs shouldn't count as 50 independent votes). Correct
                    basis for a *conservation* measure, so entropy uses it.
      - unweighted: the plain observed frequency, i.e. "how often does nature
                    actually use residue X here". Correct basis for judging
                    whether a candidate substitution is novel, since weighting
                    systematically depresses the majority residue's apparent
                    share (measured on this dataset: mean WT frequency 0.42
                    unweighted vs 0.33 weighted) and would overstate novelty.
    """
    weighted = {}
    counts = {}
    total_weight = 0.0
    n_observed = 0
    for r, w in zip(column, weights):
        if r in GAP_CHARS:
            continue
        weighted[r] = weighted.get(r, 0.0) + w
        counts[r] = counts.get(r, 0) + 1
        total_weight += w
        n_observed += 1

    if total_weight == 0 or n_observed == 0:
        return {}, {}, np.nan

    weighted_freqs = {r: w / total_weight for r, w in weighted.items()}
    unweighted_freqs = {r: c / n_observed for r, c in counts.items()}
    probs = np.array(list(weighted_freqs.values()))
    entropy = float(-(probs * np.log2(probs)).sum())
    return weighted_freqs, unweighted_freqs, entropy


def align_and_score(combined_fasta: Path, reference_id: str, out_dir: Path) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)

    # combine_with_query() always writes the query record(s) first, before any
    # homologs -- record it here so we can verify MAFFT preserved that order,
    # rather than searching for reference_id inside headers (PDB-derived FASTA
    # headers don't contain the UniProt accession, so that search can silently
    # fail depending on header format).
    query_id_before_alignment = next(SeqIO.parse(combined_fasta, "fasta")).id

    aln_path = run_mafft(combined_fasta, out_dir / "alignment.fasta")
    alignment = AlignIO.read(aln_path, "fasta")

    ref_record = alignment[0]
    if ref_record.id != query_id_before_alignment:
        raise RuntimeError(
            f"Alignment record order mismatch: expected the first aligned record to be "
            f"the query '{query_id_before_alignment}', but got '{ref_record.id}'. "
            f"MAFFT is expected to preserve input order by default (no --reorder flag "
            f"used here) -- if this fires, something upstream changed the ordering."
        )

    print(f"Reference record: '{ref_record.id}' (expected UniProt {reference_id})")
    weights = _sequence_weights(alignment)

    rows = []
    ref_pos = 0  # 1-indexed residue number in the ungapped reference sequence
    for col in range(alignment.get_alignment_length()):
        ref_char = ref_record.seq[col]
        if ref_char in GAP_CHARS:
            continue
        ref_pos += 1
        column_residues = [alignment[i, col] for i in range(len(alignment))
                            if alignment[i, col] not in GAP_CHARS]
        gap_fraction = 1 - (len(column_residues) / len(alignment))
        top_residue = pd.Series(column_residues).mode().iloc[0] if column_residues else None
        freqs, unweighted_freqs, entropy = _weighted_residue_distribution(alignment[:, col], weights)

        row = {
            "residue_number": ref_pos,
            "wt_residue": ref_char,
            "cons_shannon_entropy": entropy,
            "cons_gap_fraction": gap_fraction,
            "cons_n_sequences_observed": len(column_residues),
            "cons_predominant_residue": top_residue,
            "cons_is_conserved_as_wt": top_residue == ref_char,
            "cons_wt_freq_unweighted": unweighted_freqs.get(ref_char, 0.0),
        }
        for aa in ALL_AAS:
            row[f"cons_freq_{aa}"] = freqs.get(aa, 0.0)
            row[f"cons_freq_unweighted_{aa}"] = unweighted_freqs.get(aa, 0.0)
        rows.append(row)

    df = pd.DataFrame(rows)
    df = add_relative_conservation(df)
    df.to_csv(out_dir / "conservation_scores.csv", index=False)
    return df


def add_relative_conservation(df: pd.DataFrame, min_n_seq: int = 10) -> pd.DataFrame:
    """
    Raw Shannon entropy isn't comparable across proteins or alignments -- it
    depends on depth and diversity of the homolog set. Score each position
    *relative to other positions in this same protein* instead, the way
    ConSurf's grade is a per-protein normalized rank, not a fixed threshold.
    """
    df = df.copy()
    valid = df.cons_shannon_entropy.notna()

    df.loc[valid, "cons_conservation_percentile"] = (
        100 * (1 - df.loc[valid, "cons_shannon_entropy"].rank(pct=True))
    )
    df.loc[valid, "cons_consurf_like_grade"] = pd.qcut(
        df.loc[valid, "cons_shannon_entropy"], q=9, labels=range(9, 0, -1), duplicates="drop"
    ).astype("Int64")
    df["cons_low_confidence"] = df.cons_n_sequences_observed < min_n_seq
    return df


def load_query_metadata(homolog_search_out_dir: Path) -> dict:
    """Read the uniprot_id/pdb_id/chain that 02_search_homologs.py wrote out."""
    meta_file = homolog_search_out_dir / "query_metadata.txt"
    if not meta_file.exists():
        raise FileNotFoundError(f"{meta_file} not found -- run 02_search_homologs.py first.")
    meta = {}
    for line in meta_file.read_text().splitlines():
        key, _, value = line.partition("=")
        meta[key] = value
    return meta


def run_from_config(config_path: str = None) -> pd.DataFrame:
    cfg = load_config(config_path)
    results_dir = results_dir_for(cfg)
    homolog_out_dir = results_dir / "02_homolog_search"

    meta = load_query_metadata(homolog_out_dir)
    combined_fasta = homolog_out_dir / "query_plus_homologs.fasta"
    out_dir = results_dir / "03_conservation"

    return align_and_score(combined_fasta, reference_id=meta["uniprot_id"], out_dir=out_dir)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, default=None,
                         help="Override the active target set in config/base.yaml.")
    args = parser.parse_args()

    df = run_from_config(args.config)
    print(f"Wrote {len(df)} residue rows to conservation_scores.csv")
    print(df[["residue_number", "wt_residue", "cons_shannon_entropy",
              "cons_conservation_percentile", "cons_low_confidence"]].describe(include="all"))


if __name__ == "__main__":
    main()
