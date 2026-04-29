#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path

import yaml


def load_fasta_sequence(path: Path) -> str:
    lines = path.read_text().splitlines()
    seq = "".join(line.strip() for line in lines if line and not line.startswith(">"))
    if not seq:
        raise ValueError(f"No sequence found in FASTA: {path}")
    return seq


def slugify(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def build_payload(receptor_seq: str, ligand_type: str, value: str) -> dict:
    payload = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": receptor_seq}},
        ],
    }

    if ligand_type == "small_molecule":
        payload["sequences"].append({"ligand": {"id": "B", "smiles": value}})
        payload["properties"] = [{"affinity": {"binder": "B"}}]
    elif ligand_type == "peptide":
        payload["sequences"].append({"protein": {"id": "B", "sequence": value}})
    else:
        raise ValueError(f"Unsupported ligand type: {ligand_type}")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate Boltz-2 YAML inputs from CSV.")
    parser.add_argument("--library", default="data/ligand_library.csv")
    parser.add_argument("--target-fasta", default="data/target_ghsr.fasta")
    parser.add_argument("--out-dir", default="inputs")
    args = parser.parse_args()

    library_path = Path(args.library)
    target_fasta_path = Path(args.target_fasta)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    receptor_seq = load_fasta_sequence(target_fasta_path)

    rows = list(csv.DictReader(library_path.open(newline="", encoding="utf-8")))
    written = 0
    skipped = 0

    for row in rows:
        name = (row.get("name") or "").strip()
        ligand_type = (row.get("type") or "").strip().lower()
        smiles = (row.get("smiles") or "").strip()
        sequence = (row.get("sequence") or "").strip()
        value = smiles if ligand_type == "small_molecule" else sequence

        if not name or not ligand_type:
            skipped += 1
            continue
        if ligand_type not in {"small_molecule", "peptide"}:
            skipped += 1
            continue
        if not value:
            skipped += 1
            continue

        ligand_id = slugify(name)
        payload = build_payload(receptor_seq, ligand_type, value)
        out_path = out_dir / f"{ligand_id}.yaml"
        out_path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
        written += 1

    print(f"Wrote {written} YAML files to {out_dir}")
    print(f"Skipped {skipped} entries with missing/invalid fields")


if __name__ == "__main__":
    main()
