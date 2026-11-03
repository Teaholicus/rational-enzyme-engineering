"""
project_utils.py

Shared helpers so every pipeline script resolves paths the same way,
anchored to the repo root (not the caller's working directory), and so
the query FASTA is found from pdb_id/chain alone -- config.yaml never
needs to know the UniProt accession or a hardcoded file path.
"""
import re
from pathlib import Path

# workflow/project_utils.py -> repo root is one level up
PROJECT_ROOT = Path(__file__).resolve().parents[1]


def resolve_path(cfg_value: str) -> Path:
    """
    Resolve a config value relative to the project root, regardless of the
    working directory the script was launched from.
    """
    p = Path(cfg_value)
    return p if p.is_absolute() else (PROJECT_ROOT / p).resolve()


def protein_tag(cfg: dict) -> str:
    """Short identifier for the current target, used to namespace outputs."""
    return f"{cfg['protein']['pdb_id'].upper()}_{cfg['protein']['chain'].upper()}"


def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override into base; override wins on conflicts."""
    out = dict(base)
    for key, value in override.items():
        if key in out and isinstance(out[key], dict) and isinstance(value, dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def load_config(config_path: str = None) -> dict:
    """
    Load configuration for the currently active target.

    Layout:
        config/base.yaml        shared settings + `active_protein: "1CTN_A"`
        config/1CTN_A.yaml      per-protein: pdb_id, chain, functional residues,
        config/2JIE_A.yaml      and any per-protein overrides

    Switching targets means editing ONE line (`active_protein`) in base.yaml --
    no script edits, no config file swapping. Every step then runs against the
    same target automatically:

        python workflow/03_conservation_analysis.py

    An explicit path still overrides, for one-off runs against a target that
    isn't the active one:

        python workflow/03_conservation_analysis.py --config config/2JIE_A.yaml

    Per-protein values win over base values on conflict.
    """
    import yaml

    base_path = resolve_path("config/base.yaml")
    base_cfg = {}
    if base_path.exists():
        base_cfg = yaml.safe_load(base_path.read_text()) or {}

    if config_path:
        protein_cfg_path = resolve_path(config_path)
    else:
        # Which target is active lives in its own file, separate from the shared
        # settings in base.yaml. They change for different reasons and by
        # different hands: base.yaml carries pipeline parameters, active.yaml
        # carries local working state. Keeping them apart means editing one
        # cannot clobber the other.
        active_path = base_path.parent / "active.yaml"
        active = None
        if active_path.exists():
            active_cfg = yaml.safe_load(active_path.read_text()) or {}
            active = active_cfg.get("active_protein")
        if not active:
            active = base_cfg.get("active_protein")   # older layout
        if not active:
            available = sorted(p.stem for p in base_path.parent.glob("*.yaml")
                                if p.stem not in {"base", "active"})
            raise KeyError(
                f"No target selected. Set `active_protein` in config/active.yaml "
                f"(e.g. active_protein: \"2JIE_A\") or pass --config explicitly.\n"
                f"Available targets: {available}")
        protein_cfg_path = base_path.parent / f"{active}.yaml"

    if not protein_cfg_path.exists():
        available = sorted(p.stem for p in base_path.parent.glob("*.yaml") if p.stem != "base")
        raise FileNotFoundError(
            f"Config not found: {protein_cfg_path}\nAvailable targets: {available}")

    protein_cfg = yaml.safe_load(protein_cfg_path.read_text()) or {}
    merged = _deep_merge(base_cfg, protein_cfg) if base_cfg else protein_cfg
    merged.pop("active_protein", None)   # bookkeeping key, not pipeline config

    for required in ("protein", "paths"):
        if required not in merged:
            raise KeyError(f"Config missing required section '{required}': {protein_cfg_path}")

    print(f"[config] target: {protein_tag(merged)}  (from {protein_cfg_path.name})")
    return merged


def results_dir_for(cfg: dict) -> Path:
    """
    Protein-scoped results directory: results/<PDBID>_<CHAIN>/.
    Keeps multiple targets from overwriting each other's tables, so the
    pipeline can be re-run on a new protein by changing only pdb_id/chain.
    """
    d = resolve_path(cfg["paths"]["results_dir"]) / protein_tag(cfg)
    d.mkdir(parents=True, exist_ok=True)
    return d


def structures_dir_for(cfg: dict) -> Path:
    """Protein-scoped structure directory: data/structures/<PDBID>_<CHAIN>/.
    Falls back to the shared directory if the protein-scoped one doesn't
    exist yet, so existing single-protein layouts keep working."""
    base = resolve_path(cfg["paths"].get("structures_dir", "data/structures"))
    scoped = base / protein_tag(cfg)
    return scoped if scoped.exists() else base


def find_query_fasta(pdb_id: str, chain: str,
                      sequences_dir: str = "data/sequences") -> Path:
    """
    Locate the query FASTA written by 01_acquire_protein.py, from pdb_id and
    chain alone -- so switching targets needs no path edits anywhere.

    Expects filenames of the form: UNIPROTID(PDBID_CHAIN).fasta
    e.g. P07254(1CTN_A).fasta, P22505(2JIE_A).fasta
    """
    seq_dir = resolve_path(sequences_dir)
    pattern = f"*({pdb_id.upper()}_{chain.upper()}).fasta"
    matches = list(seq_dir.glob(pattern))

    if not matches:
        available = [p.name for p in seq_dir.glob("*.fasta")]
        raise FileNotFoundError(
            f"No FASTA matching pattern '{pattern}' in {seq_dir}.\n"
            f"Available files: {available}"
        )
    if len(matches) > 1:
        raise ValueError(f"Ambiguous match for '{pattern}': {[m.name for m in matches]}")
    return matches[0]


def uniprot_id_from_fasta_path(fasta_path: Path) -> str:
    """Extract 'P07254' from 'P07254(1CTN_A).fasta'."""
    match = re.match(r"^([A-Za-z0-9]+)\(", fasta_path.stem)
    if not match:
        raise ValueError(f"Could not parse UniProt ID from filename: {fasta_path.name}")
    return match.group(1)


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Histidine tautomer codes some tools (FoldX included) assign based on
    # the local H-bond environment -- X-ray structures don't resolve
    # hydrogens, so the protonation state is inferred, not observed.
    "H1S": "H", "H2S": "H", "HID": "H", "HIE": "H", "HIP": "H",
    "HSD": "H", "HSE": "H", "HSP": "H",
}


def parse_residue_code(code: str) -> tuple[str, int]:
    """Parse a config-style residue code like 'D313' into ('D', 313)."""
    match = re.match(r"^([A-Za-z])(\d+)$", code.strip())
    if not match:
        raise ValueError(f"Could not parse residue code: {code!r} (expected e.g. 'D313')")
    return match.group(1).upper(), int(match.group(2))


class ChainOnlySelector:
    """
    Bio.PDB.PDBIO Select implementation: keep one chain, standard amino acid
    residues only (excludes waters/ligands/other heteroatoms -- residue.id[0]
    == " " is Biopython's convention for "not a HETATM record").
    """
    def __init__(self, chain_id: str):
        self.chain_id = chain_id

    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def accept_residue(self, residue):
        return residue.id[0] == " "

    def accept_atom(self, atom):
        return True

    def accept_model(self, model):
        return True


def cif_to_pdb(cif_path: Path, pdb_out: Path, chain_id: str) -> Path:
    """Convert mmCIF to legacy PDB format, one chain, standard residues only.
    Several tools (FoldX, FreeSASA-by-file) expect legacy PDB, not mmCIF."""
    from Bio.PDB import MMCIFParser, PDBIO
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("query", str(cif_path))
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdb_out), ChainOnlySelector(chain_id))

    # Biopython's PDBIO doesn't write a HEADER line. Several tools (confirmed:
    # mkdssp) auto-detect PDB vs mmCIF by checking whether the first line
    # starts with "HEADER" -- without it, they wrongly assume mmCIF and fail
    # to parse what actually is a valid PDB file. Prepend a minimal one.
    content = pdb_out.read_text()
    if not content.startswith("HEADER"):
        header_line = f"HEADER    {pdb_out.stem[:40]:<40}01-JAN-00   {pdb_out.stem[:4].upper():>4}\n"
        pdb_out.write_text(header_line + content)

    return pdb_out
