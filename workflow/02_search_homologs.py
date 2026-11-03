"""
02_search_homologs.py

Search a query protein against a local MMseqs2 UniRef90 database, filter hits
by coverage/identity, deduplicate near-identical sequences, and write a
combined FASTA (query + homologs) ready for MAFFT alignment in step 03.

Usage:
    python workflow/02_search_homologs.py
    python workflow/02_search_homologs.py --config config/2JIE_A.yaml

Target defaults to `active_protein` in config/base.yaml. Can be run from any
working directory -- all paths resolve relative to the repo root.
"""
import argparse
import subprocess
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from project_utils import load_config, find_query_fasta, resolve_path, results_dir_for, uniprot_id_from_fasta_path

M8_COLUMNS = ["query", "target", "pident", "alnlen", "mismatch", "gapopen",
              "qstart", "qend", "tstart", "tend", "evalue", "bits", "qcov", "tcov"]


def run_mmseqs_search(query_fasta: Path, db_path: Path, out_dir: Path,
                       sensitivity: float = 7.5, evalue: float = 1e-5,
                       max_seqs: int = 1000) -> Path:
    """Run mmseqs easy-search, return path to the .m8 results table."""
    out_dir.mkdir(parents=True, exist_ok=True)
    result_m8 = out_dir / "search_results.m8"
    tmp_dir = out_dir / "tmp"  # scratch space lives under results/, never /tmp
    tmp_dir.mkdir(exist_ok=True)

    cmd = [
        "mmseqs", "easy-search", str(query_fasta), str(db_path), str(result_m8), str(tmp_dir),
        "-s", str(sensitivity),
        "-e", str(evalue),
        "--max-seqs", str(max_seqs),
        "--format-output", ",".join(M8_COLUMNS),
        "--threads", "8",
    ]
    subprocess.run(cmd, check=True)
    return result_m8


def filter_hits(m8_path: Path, min_pident: float = 0.40, min_qcov: float = 0.7) -> pd.DataFrame:
    """
    Keep hits that are plausible orthologs rather than distant family members.

    On min_pident: the 20-35% identity band is the "twilight zone", where
    pairwise identity no longer reliably indicates homology and alignment
    quality degrades badly. Including it pulls in the whole enzyme family
    (for a GH18 chitinase: insect/fungal chitinases and chitinase-like
    proteins with different substrates), which answers "what does this fold
    require?" rather than "what has evolution tuned in this specific
    enzyme?". The latter is the question relevant to engineering targets,
    so the default sits above the twilight zone.

    Redundancy is NOT handled here -- see cluster_hits(), which does it by
    actual sequence clustering. (An earlier version compared hits' pident
    *values* to each other, which cannot detect redundancy at all: two
    unrelated sequences that happen to share the same identity to the query
    looked like duplicates, while a clade of near-identical homologs at
    slightly different identities all passed through.)
    """
    df = pd.read_csv(m8_path, sep="\t", names=M8_COLUMNS)
    n_raw = len(df)
    df = df[(df.pident >= min_pident * 100) & (df.qcov >= min_qcov)]
    print(f"  {n_raw} raw hits -> {len(df)} after pident>={min_pident:.0%}, qcov>={min_qcov:.0%}")
    return df.sort_values("bits", ascending=False)


def cluster_hits(homologs_fasta: Path, out_dir: Path, min_seq_id: float = 0.90) -> Path:
    """
    Collapse near-identical homologs by actual sequence clustering, keeping
    one representative per cluster. Redundant clades otherwise dominate the
    alignment and distort per-column frequencies.
    """
    tmp_dir = out_dir / "cluster_tmp"
    tmp_dir.mkdir(exist_ok=True)
    rep_fasta = out_dir / "homologs_clustered.fasta"

    subprocess.run([
        "mmseqs", "easy-cluster", str(homologs_fasta), str(out_dir / "clust"), str(tmp_dir),
        "--min-seq-id", str(min_seq_id), "-c", "0.8", "--cov-mode", "0", "--threads", "8",
    ], check=True)

    produced = out_dir / "clust_rep_seq.fasta"
    if not produced.exists():
        raise FileNotFoundError(f"Expected clustering output not found: {produced}")
    produced.replace(rep_fasta)

    n_before = sum(1 for _ in SeqIO.parse(homologs_fasta, "fasta"))
    n_after = sum(1 for _ in SeqIO.parse(rep_fasta, "fasta"))
    print(f"  Clustering at {min_seq_id:.0%} identity: {n_before} -> {n_after} representatives")
    return rep_fasta


def extract_sequences(filtered_df: pd.DataFrame, db_path: Path, out_dir: Path) -> Path:
    """Pull full-length sequences for the filtered target IDs out of the local DB."""
    id_file = out_dir / "keep_ids.txt"
    id_file.write_text("\n".join(filtered_df.target.unique()))

    subset_db = out_dir / "subset_db"
    subprocess.run(["mmseqs", "createsubdb", str(id_file), str(db_path), str(subset_db),
                     "--id-mode", "1"], check=True)

    homologs_fasta = out_dir / "homologs.fasta"
    subprocess.run(["mmseqs", "convert2fasta", str(subset_db), str(homologs_fasta)], check=True)
    return homologs_fasta


def combine_with_query(query_fasta: Path, homologs_fasta: Path, out_path: Path) -> Path:
    records = list(SeqIO.parse(query_fasta, "fasta")) + list(SeqIO.parse(homologs_fasta, "fasta"))
    SeqIO.write(records, out_path, "fasta")
    return out_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, default=None,
                         help="Override the active target set in config/base.yaml.")
    args = parser.parse_args()

    cfg = load_config(args.config)

    query_fasta = find_query_fasta(
        cfg["protein"]["pdb_id"], cfg["protein"]["chain"],
        sequences_dir=cfg["paths"].get("sequences_dir", "data/sequences"))
    uniprot_id = uniprot_id_from_fasta_path(query_fasta)
    print(f"Query: {query_fasta.name}  (UniProt {uniprot_id})")

    db_path = resolve_path(cfg["homolog_search"]["mmseqs_db_path"])
    out_dir = results_dir_for(cfg) / "02_homolog_search"

    m8 = run_mmseqs_search(query_fasta, db_path, out_dir,
                            max_seqs=cfg["homolog_search"].get("max_seqs", 2000))
    filtered = filter_hits(m8,
                            min_pident=cfg["homolog_search"].get("min_pident", 0.40),
                            min_qcov=cfg["homolog_search"].get("min_qcov", 0.7))
    filtered.to_csv(out_dir / "filtered_hits.csv", index=False)

    homologs_fasta = extract_sequences(filtered, db_path, out_dir)
    clustered_fasta = cluster_hits(
        homologs_fasta, out_dir,
        min_seq_id=cfg["homolog_search"].get("cluster_min_seq_id", 0.90))

    combined = combine_with_query(query_fasta, clustered_fasta, out_dir / "query_plus_homologs.fasta")
    n_final = sum(1 for _ in SeqIO.parse(combined, "fasta"))
    print(f"Wrote combined FASTA for alignment ({n_final} sequences incl. query): {combined}")

    # Small metadata file so step 03 (and the notebook) can find the same
    # reference_id without re-deriving it from config.
    (out_dir / "query_metadata.txt").write_text(
        f"uniprot_id={uniprot_id}\npdb_id={cfg['protein']['pdb_id']}\nchain={cfg['protein']['chain']}\n"
    )


if __name__ == "__main__":
    main()
