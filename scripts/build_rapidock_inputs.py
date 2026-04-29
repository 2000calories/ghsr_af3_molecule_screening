#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate RAPiDock protein_peptide.csv from ligand_library.csv peptide rows."
    )
    parser.add_argument("--library", default="data/ligand_library.csv")
    parser.add_argument("--receptor-pdb", default="data/ghsr_pocket.pdb")
    parser.add_argument("--out-dir", default="inputs_rapidock")
    parser.add_argument("--out-csv-name", default="protein_peptide.csv")
    args = parser.parse_args()

    library_path = Path(args.library)
    receptor_pdb = Path(args.receptor_pdb)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / args.out_csv_name

    rows = list(csv.DictReader(library_path.open(newline="", encoding="utf-8")))

    out_rows = []
    skipped = 0
    for row in rows:
        ligand_type = (row.get("type") or "").strip().lower()
        if ligand_type != "peptide":
            continue

        name = (row.get("name") or "").strip()
        if not name:
            skipped += 1
            continue

        seq = (row.get("peptide_seq_rapidock") or "").strip() or (row.get("sequence") or "").strip()
        if not seq:
            skipped += 1
            continue

        # Flag clearly unsupported ghrelin octanoyl proxy or other PTMs if they slip in as peptides.
        if "octanoyl" in (row.get("notes") or "").lower():
            skipped += 1
            continue

        complex_name = name.lower().replace("-", "_").replace(" ", "_")
        out_rows.append(
            {
                "complex_name": complex_name,
                "protein_description": str(receptor_pdb.resolve()),
                "peptide_description": seq,
            }
        )

    fieldnames = ["complex_name", "protein_description", "peptide_description"]
    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in out_rows:
            writer.writerow(r)

    print(f"Wrote {len(out_rows)} RAPiDock peptide rows to {out_csv}")
    if skipped:
        print(f"Skipped {skipped} peptide entries lacking sequences or unsupported PTMs.")


if __name__ == "__main__":
    main()

