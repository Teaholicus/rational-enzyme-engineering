"""
util_clean_foldx_artifacts.py

Remove the per-mutation PDB and energies files PositionScan leaves in the FoldX
directory. A full-protein scan writes one of each per mutation (~9000 files for
a 450-residue protein), all regenerable from the PS_*_scanning_output.txt that
the pipeline actually parses.

04_foldx_ddg.py now cleans up after itself; this is for clearing leftovers from
earlier runs.

Dry run first (default -- lists what would go, deletes nothing):
    python workflow/util_clean_foldx_artifacts.py

Then actually delete:
    python workflow/util_clean_foldx_artifacts.py --execute

Or keep them, moved out of the way:
    python workflow/util_clean_foldx_artifacts.py --execute --archive-to results/foldx_artifacts

Only files matching FoldX's own output patterns are considered. rotabase.txt,
the binary, molecules/, and anything else in the directory are never touched.
"""
import argparse
import re
from pathlib import Path

from project_utils import load_config, resolve_path

# Files FoldX writes per mutation / per run. Anchored patterns so that unrelated
# files sharing a substring are not matched.
ARTIFACT_PATTERNS = [
    re.compile(r"^[A-Z]{3}\d+[A-Za-z]?_.*_Repair\.pdb$"),   # ALA4_2JIE_A_Repair.pdb
    re.compile(r"^WT_[A-Z]{3}\d+[A-Za-z]?_.*_Repair\.pdb$"),
    re.compile(r"^energies_\d+_.*_Repair\.txt$"),            # energies_4_2JIE_A_Repair.txt
    re.compile(r"^(Raw|Dif|Average|PdbList)_.*\.(fxout|txt)$"),
]

# Never remove these even if a pattern somehow matched.
PROTECTED = {"rotabase.txt", "foldx.exe", "foldx5.exe", "yasaraPlugin.zip"}


def find_artifacts(foldx_dir: Path) -> list:
    hits = []
    for p in foldx_dir.iterdir():
        if not p.is_file() or p.name in PROTECTED:
            continue
        if any(pat.match(p.name) for pat in ARTIFACT_PATTERNS):
            hits.append(p)
    return sorted(hits)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, default=None)
    parser.add_argument("--foldx-dir", type=str, default=None,
                         help="Override the foldx directory from config.")
    parser.add_argument("--execute", action="store_true",
                         help="Actually remove files. Without this, only reports.")
    parser.add_argument("--archive-to", type=str, default=None,
                         help="Move files here instead of deleting them.")
    args = parser.parse_args()

    if args.foldx_dir:
        foldx_dir = resolve_path(args.foldx_dir)
    else:
        cfg = load_config(args.config)
        foldx_dir = resolve_path(cfg["foldx"]["foldx_dir"])

    if not foldx_dir.exists():
        raise FileNotFoundError(f"FoldX directory not found: {foldx_dir}")

    artifacts = find_artifacts(foldx_dir)
    total_mb = sum(p.stat().st_size for p in artifacts) / (1024 * 1024)

    print(f"FoldX directory: {foldx_dir}")
    print(f"Artifact files found: {len(artifacts)} ({total_mb:.1f} MB)")
    if not artifacts:
        return

    print("\nSample (first 10):")
    for p in artifacts[:10]:
        print(f"  {p.name}")
    if len(artifacts) > 10:
        print(f"  ... and {len(artifacts) - 10} more")

    kept = [p.name for p in foldx_dir.iterdir() if p.is_file() and p not in artifacts]
    print(f"\nWill be left untouched ({len(kept)} file(s)): {sorted(kept)[:12]}")

    if not args.execute:
        print("\nDry run -- nothing removed. Re-run with --execute to proceed.")
        return

    if args.archive_to:
        dest = resolve_path(args.archive_to)
        dest.mkdir(parents=True, exist_ok=True)
        for p in artifacts:
            p.replace(dest / p.name)
        print(f"\nMoved {len(artifacts)} file(s) -> {dest}")
    else:
        for p in artifacts:
            p.unlink()
        print(f"\nDeleted {len(artifacts)} file(s), freed {total_mb:.1f} MB.")


if __name__ == "__main__":
    main()
