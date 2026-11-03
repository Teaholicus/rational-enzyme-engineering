"""
04_foldx_ddg.py

Compute FoldX ddG (stability change on mutation) for every amino-acid
substitution at a set of target positions, using RepairPDB + PositionScan.

FoldX ships as a native Windows executable, so this script auto-detects
whether it's running under WSL and, if so, converts paths passed to the
FoldX subprocess into Windows-style paths via `wslpath` -- FoldX has no
concept of /mnt/n/..., only WSL's process-launch interop is transparent,
not path arguments. Works unmodified whether launched from WSL or from
native Windows Python.

Run (from WSL, same pattern as 02/03):
    python workflow/04_foldx_ddg.py
    python workflow/04_foldx_ddg.py --config config/2JIE_A.yaml

BEFORE TRUSTING A FULL RUN:
  1. CLI flags (--pdb-dir, --output-dir, --positions, --output-file) are
     taken from FoldX's own documentation but flag support drifts across
     versions/builds. If RepairPDB fails immediately, run
     `foldx5.exe --help` (from tools/foldx5) to confirm flag names on your
     build and adjust the cmd lists below.
  2. PositionScan's output column layout is parsed defensively -- anything
     that doesn't match the expected mutation-code pattern is logged to
     results/04_foldx_ddg/unparsed_lines.txt rather than dropped silently.
  3. Do a smoke test with 2-3 positions (see target_positions in
     config.yaml) before committing to a full-protein scan -- this also
     confirms the WSL/Windows cwd + path handling works on your machine
     before you spend hours on a full scan.

Output columns:
    residue_number, wt_residue           -- join keys, same convention as 03
    foldx_ddg_A ... foldx_ddg_Y           -- ddG (kcal/mol) per target amino
                                             acid, positive = destabilizing
"""
import argparse
import platform
import re
import subprocess
from pathlib import Path

import pandas as pd

from project_utils import load_config, resolve_path, results_dir_for, structures_dir_for, AA3_TO_1, cif_to_pdb

ALL_AAS = list("ACDEFGHIKLMNPQRSTVWY")

# FoldX distinguishes histidine protonation states as *mutation targets* too:
# 'o' = delta tautomer, 'e' = epsilon tautomer. Same H1S/H2S concept, just in
# the single-letter mutant-code field instead of the 3-letter WT-code field.
MUTANT_CODE_TO_AA = {"o": "H", "e": "H"}


def build_mutation_regex(chain: str) -> re.Pattern:
    """
    FoldX echoes mutation codes as: <3-char WT code><chain><position><mutant>,
    e.g. 'ALAA24N' or 'H1SA384o'. The WT code is usually 3 uppercase letters
    but can contain a digit for histidine tautomers (H1S/H2S) -- hence
    [A-Z0-9]{3}, not [A-Z]{3}. Chain must be matched against the known
    literal value (from config), not treated as a wildcard -- both chain and
    position can contain digits/letters, so an ambiguous pattern silently
    misparses the position (verified: it swallows the leading digit of a
    2-digit position into the "chain" group otherwise).
    """
    return re.compile(rf"^([A-Z0-9]{{3}})({re.escape(chain)})(\d+)([A-Za-z])$")


def _running_under_wsl() -> bool:
    if platform.system() != "Linux":
        return False
    try:
        return "microsoft" in Path("/proc/version").read_text().lower()
    except FileNotFoundError:
        return False


def to_foldx_path(path: Path) -> str:
    """Convert to a Windows-style path if we're calling FoldX from WSL."""
    path = Path(path).resolve()
    if _running_under_wsl():
        result = subprocess.run(["wslpath", "-w", str(path)],
                                 capture_output=True, text=True, check=True)
        return result.stdout.strip()
    return str(path)


def run_repair(pdb_path: Path, foldx_dir: Path, foldx_exe: str, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    repaired = out_dir / f"{pdb_path.stem}_Repair.pdb"
    if repaired.exists():
        print(f"Repaired structure already exists, skipping RepairPDB: {repaired}")
        return repaired

    cmd = [
        str(foldx_dir / foldx_exe),
        "--command=RepairPDB",
        f"--pdb={pdb_path.name}",
        f"--pdb-dir={to_foldx_path(pdb_path.parent)}",
        f"--output-dir={to_foldx_path(out_dir)}",
    ]
    subprocess.run(cmd, cwd=str(foldx_dir), check=True)

    if not repaired.exists():
        raise FileNotFoundError(
            f"Expected repaired PDB not found: {repaired}. "
            f"Check the FoldX output above -- likely a CLI flag mismatch, see module docstring."
        )
    return repaired


def get_structure_residues(pdb_path: Path, chain_id: str) -> dict:
    """
    Map residue_number -> 1-letter wt_residue directly from the structure.
    This is the ONLY correct source for building FoldX position codes:
    using an external sequence source (e.g. UniProt) here can silently
    generate a WT letter that disagrees with what's actually in the
    structure -- FoldX then rejects that position outright with zero
    output and no visible error, which is worse than a wrong answer
    (confirmed: this is exactly what caused 6 positions to vanish
    entirely from an earlier full-protein run without any error).
    """
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("repaired", str(pdb_path))
    chain = structure[0][chain_id]
    return {res.id[1]: AA3_TO_1.get(res.get_resname(), "X")
            for res in chain if res.id[0] == " "}


def build_position_codes(target_positions, repaired_pdb: Path, chain: str) -> list[str]:
    """
    Build FoldX PositionScan codes (e.g. 'MA1a') using wt_residue read
    directly from the repaired structure.

    target_positions: "all" (every modeled residue) or an explicit list of
    residue numbers (e.g. for a targeted rescan/patch).
    """
    structure_residues = get_structure_residues(repaired_pdb, chain)

    if target_positions == "all":
        position_numbers = sorted(structure_residues.keys())
    else:
        position_numbers = list(target_positions)
        missing = set(position_numbers) - set(structure_residues.keys())
        if missing:
            print(f"Warning: requested positions not found in the structure "
                  f"(unmodeled/missing density): {sorted(missing)}")
            position_numbers = [p for p in position_numbers if p in structure_residues]

    return [f"{structure_residues[pos]}{chain}{pos}a" for pos in position_numbers]


def cleanup_foldx_artifacts(foldx_dir: Path, pdb_stem: str, out_dir: Path,
                             mode: str = "delete") -> int:
    """
    Remove (or archive) the per-mutation files PositionScan writes into its
    working directory.

    FoldX writes one mutant PDB and one energies file per mutation regardless
    of --output-dir, so a full scan leaves thousands of files in tools/foldx5.
    They are all regenerable from the scan output that we actually parse, so
    the default is to delete them.

    Matching is deliberately narrow -- only files carrying this run's PDB stem
    are touched, so rotabase.txt, the binary, and any other target's leftovers
    are never at risk.

    mode: "delete" | "archive" (move under out_dir) | "keep" (no-op)
    """
    if mode == "keep":
        return 0

    patterns = [
        f"*_{pdb_stem}.pdb",            # e.g. ALA4_2JIE_A_Repair.pdb
        f"energies_*_{pdb_stem}.txt",   # e.g. energies_4_2JIE_A_Repair.txt
        f"Raw_*_{pdb_stem}.*",
        f"Dif_*_{pdb_stem}.*",
        f"Average_*_{pdb_stem}.*",
        f"PdbList_*_{pdb_stem}.*",
    ]
    targets = []
    for pattern in patterns:
        targets.extend(foldx_dir.glob(pattern))
    # never touch the repaired structure we depend on
    targets = [p for p in targets if p.name != f"{pdb_stem}.pdb"]

    if not targets:
        return 0

    if mode == "archive":
        archive_dir = out_dir / "foldx_artifacts"
        archive_dir.mkdir(parents=True, exist_ok=True)
        for p in targets:
            p.replace(archive_dir / p.name)
        print(f"Archived {len(targets)} FoldX artifact file(s) -> {archive_dir}")
    else:
        for p in targets:
            p.unlink()
        print(f"Deleted {len(targets)} FoldX artifact file(s) from {foldx_dir} "
              f"(regenerable from the scan output).")
    return len(targets)


def run_position_scan(repaired_pdb: Path, foldx_dir: Path, foldx_exe: str,
                       positions: list[str], out_dir: Path, tag: str,
                       force: bool = False) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    result_file = out_dir / f"PS_{tag}_scanning_output.txt"

    if result_file.exists() and not force:
        print(f"PositionScan output already exists, skipping the scan (pass force=True to redo): {result_file}")
        return result_file

    cmd = [
        str(foldx_dir / foldx_exe),
        "--command=PositionScan",
        f"--pdb={repaired_pdb.name}",
        f"--pdb-dir={to_foldx_path(repaired_pdb.parent)}",
        f"--output-dir={to_foldx_path(out_dir)}",
        f"--positions={','.join(positions)}",
        f"--output-file={tag}",
    ]
    subprocess.run(cmd, cwd=str(foldx_dir), check=True)

    if not result_file.exists():
        candidates = list(out_dir.glob("PS_*scanning_output*"))
        if not candidates:
            raise FileNotFoundError(f"No PositionScan output found in {out_dir}. See module docstring point 1.")
        result_file = candidates[0]
    return result_file


def parse_position_scan_output(result_file: Path, log_dir: Path, chain: str,
                                unparsed_filename: str = "unparsed_lines.txt") -> pd.DataFrame:
    mutation_re = build_mutation_regex(chain)
    rows, unparsed = [], []

    for line in result_file.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        parts = re.split(r"\s+", line)
        if len(parts) < 2:
            unparsed.append(line)
            continue
        code, ddg_str = parts[0], parts[1]
        match = mutation_re.match(code)
        if not match:
            unparsed.append(line)
            continue

        wt3, _chain, pos, mut_raw = match.groups()
        wt = AA3_TO_1.get(wt3)
        mut = MUTANT_CODE_TO_AA.get(mut_raw, mut_raw.upper())
        if wt is None or mut not in ALL_AAS:
            unparsed.append(f"{line}  (unrecognized residue code: wt={wt3!r} mut={mut_raw!r})")
            continue
        try:
            ddg = float(ddg_str)
        except ValueError:
            unparsed.append(line)
            continue
        rows.append({"residue_number": int(pos), "wt_residue": wt, "mut_residue": mut, "foldx_ddg": ddg})

    if unparsed:
        (log_dir / unparsed_filename).write_text("\n".join(unparsed))
        print(f"Warning: {len(unparsed)} lines could not be parsed -- see {unparsed_filename}")

    return pd.DataFrame(rows)


def pivot_to_wide(long_df: pd.DataFrame) -> pd.DataFrame:
    # aggfunc='min': collapses histidine's two protonation-state rows (see
    # MUTANT_CODE_TO_AA) into the more stable (lower ddG) value for the H column.
    wide = long_df.pivot_table(index=["residue_number", "wt_residue"],
                                columns="mut_residue", values="foldx_ddg", aggfunc="min")
    wide = wide.reindex(columns=ALL_AAS)
    wide.columns = [f"foldx_ddg_{aa}" for aa in wide.columns]
    return wide.reset_index()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=str, default=None,
                         help="Override the active target set in config/base.yaml.")
    parser.add_argument("--target-positions-override", type=str, default=None,
                         help="Comma-separated residue numbers to scan instead of config's "
                              "target_positions -- for patching specific positions without "
                              "redoing a full-protein scan.")
    parser.add_argument("--output-suffix", type=str, default="",
                         help="Appended to output filenames/tag so a targeted rerun doesn't "
                              "overwrite the main results, e.g. --output-suffix _patch")
    parser.add_argument("--artifacts", choices=["delete", "archive", "keep"], default="delete",
                         help="What to do with the per-mutation PDB/energies files FoldX "
                              "writes into its own directory. Default deletes them (they are "
                              "regenerable from the scan output); 'archive' moves them under "
                              "the results directory instead.")
    args = parser.parse_args()

    cfg = load_config(args.config)
    pdb_id = cfg["protein"]["pdb_id"]
    chain = cfg["protein"]["chain"]

    foldx_dir = resolve_path(cfg["foldx"]["foldx_dir"])
    foldx_exe = cfg["foldx"].get("foldx_exe", "foldx5.exe")

    results_dir = results_dir_for(cfg)
    out_dir = results_dir / "04_foldx_ddg"
    out_dir.mkdir(parents=True, exist_ok=True)

    cif_path = structures_dir_for(cfg) / f"{pdb_id}.cif"
    if not cif_path.exists():
        raise FileNotFoundError(f"{cif_path} not found -- run 01_acquire_protein.py first, "
                                 f"or check paths.structures_dir in config.yaml.")

    pdb_path = out_dir / f"{pdb_id}_{chain}.pdb"
    if not pdb_path.exists():
        cif_to_pdb(cif_path, pdb_path, chain)
        print(f"Converted {cif_path.name} -> {pdb_path.name}")

    repaired = run_repair(pdb_path, foldx_dir, foldx_exe, out_dir / "repair")
    print(f"Repaired structure: {repaired}")

    if args.target_positions_override:
        target_positions = [int(x) for x in args.target_positions_override.split(",")]
    else:
        target_positions = cfg["foldx"].get("target_positions", "all")

    positions = build_position_codes(target_positions, repaired, chain)
    print(f"Scanning {len(positions)} position(s) x 20 amino acids "
          f"({len(positions) * 20} mutant evaluations)...")

    tag = f"{pdb_id}{args.output_suffix}"
    result_file = run_position_scan(repaired, foldx_dir, foldx_exe, positions, out_dir / "scan", tag=tag)
    long_df = parse_position_scan_output(result_file, out_dir, chain,
                                          unparsed_filename=f"unparsed_lines{args.output_suffix}.txt")
    wide_df = pivot_to_wide(long_df)

    out_csv = out_dir / f"foldx_ddg_scores{args.output_suffix}.csv"
    wide_df.to_csv(out_csv, index=False)
    print(f"Wrote {len(wide_df)} residue rows to {out_csv.name}")

    # Only after the results are safely parsed and written.
    cleanup_foldx_artifacts(foldx_dir, repaired.stem, out_dir, mode=args.artifacts)


if __name__ == "__main__":
    main()
