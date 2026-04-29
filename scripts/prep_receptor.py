#!/usr/bin/env python3
"""
Prepare GHS-R receptor structures for RAPiDock.

This script:
- Downloads PDB 7F9Z from RCSB if missing.
- Extracts the human GHSR chain (auth chain R) and saves data/ghsr_receptor.pdb.
- Optionally runs RAPiDock's pocket_trunction.py to produce data/ghsr_pocket.pdb.

It is designed to be idempotent and can be re-run safely.
"""

import argparse
import subprocess
from pathlib import Path
from typing import Optional

import requests
from Bio.PDB import PDBIO, PDBParser, Select


RCSB_PDB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


class ChainSelect(Select):
    def __init__(self, target_chain_id: str) -> None:
        super().__init__()
        self.target_chain_id = target_chain_id

    def accept_chain(self, chain) -> bool:  # type: ignore[override]
        # auth IDs map directly to PDB chain IDs in downloaded file for 7F9Z
        return chain.id == self.target_chain_id


def download_pdb(pdb_id: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        return
    url = RCSB_PDB_URL.format(pdb_id=pdb_id.upper())
    resp = requests.get(url, timeout=30)
    resp.raise_for_status()
    dest.write_text(resp.text, encoding="utf-8")


def extract_chain(src_pdb: Path, chain_id: str, dest_pdb: Path) -> None:
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("ghsr", src_pdb)
    io = PDBIO()
    io.set_structure(structure)
    dest_pdb.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(dest_pdb), select=ChainSelect(chain_id))


def run_pocket_trunction(
    rapidock_dir: Path, receptor_pdb: Path, peptide_seq_or_pdb: Optional[str], out_pdb: Path
) -> None:
    """
    Call RAPiDock's pocket_trunction.py to generate a pocket-truncated receptor.

    peptide_seq_or_pdb is passed through as-is to the script; it can be a sequence
    string or a path to a peptide PDB, depending on how the user wants to define
    the pocket. For GHSR–GHRP-6, the full receptor PDB already encodes the bound
    peptide so we can omit this argument and let the script infer the pocket if
    supported; otherwise we fall back to using the full receptor.
    """
    pocket_script = rapidock_dir / "pocket_trunction.py"
    if not pocket_script.exists():
        raise SystemExit(f"Expected RAPiDock pocket_trunction.py at {pocket_script}, but it was not found.")

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "python",
        str(pocket_script),
        "--protein",
        str(receptor_pdb),
        "--output",
        str(out_pdb),
    ]
    if peptide_seq_or_pdb:
        cmd.extend(["--peptide", peptide_seq_or_pdb])
    subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare GHS-R receptor PDBs for RAPiDock.")
    parser.add_argument("--pdb-id", default="7F9Z", help="RCSB PDB ID to download for GHSR.")
    parser.add_argument(
        "--chain-id",
        default="R",
        help="Chain ID for the human GHSR receptor in the downloaded PDB.",
    )
    parser.add_argument(
        "--rapidock-dir",
        type=Path,
        default=None,
        help="Path to RAPiDock repo root (required if generating pocket-truncated receptor).",
    )
    parser.add_argument(
        "--skip-pocket",
        action="store_true",
        help="Skip running pocket_trunction.py and only write the full receptor PDB.",
    )
    args = parser.parse_args()

    project_root = Path(__file__).resolve().parents[1]
    data_dir = project_root / "data"
    structures_dir = data_dir / "structures"

    raw_pdb_path = structures_dir / f"{args.pdb_id.lower()}.pdb"
    receptor_pdb = data_dir / "ghsr_receptor.pdb"
    pocket_pdb = data_dir / "ghsr_pocket.pdb"

    print(f"Downloading PDB {args.pdb_id} to {raw_pdb_path} (if needed)...")
    download_pdb(args.pdb_id, raw_pdb_path)

    print(f"Extracting chain {args.chain_id} to {receptor_pdb}...")
    extract_chain(raw_pdb_path, args.chain_id, receptor_pdb)
    print(f"Wrote full receptor to {receptor_pdb}")

    if args.skip-pocket:
        print("Skipping pocket_trunction step (--skip-pocket set).")
        return

    if args.rapidock_dir is None:
        print(
            "No --rapidock-dir provided; skipping pocket_trunction. "
            f"You can still use {receptor_pdb} as the protein_description for RAPiDock."
        )
        return

    print(f"Running pocket_trunction.py from {args.rapidock_dir}...")
    try:
        run_pocket_trunction(args.rapidock_dir, receptor_pdb, peptide_seq_or_pdb=None, out_pdb=pocket_pdb)
        print(f"Wrote pocket-truncated receptor to {pocket_pdb}")
    except Exception as exc:  # pragma: no cover - defensive logging
        print(f"Warning: pocket_trunction.py failed with {exc}. You may use {receptor_pdb} instead.")


if __name__ == "__main__":
    main()

