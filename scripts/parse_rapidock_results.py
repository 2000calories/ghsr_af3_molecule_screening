#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple


R_KCAL = 1.987204258e-3
T_K = 298.15


def reu_to_ic50_uM(ref2015_reu: float) -> Tuple[float, float, float]:
    """
    Convert a Rosetta REF2015 score (assumed kcal/mol) into:
    - dG_kcal
    - log10(IC50_uM)
    - IC50_uM

    Using:
        dG = REF2015
        Kd = exp(dG / (R * T))
        IC50 ≈ Kd  (μM)
    """
    dG_kcal = float(ref2015_reu)
    Kd_M = math.exp(dG_kcal / (R_KCAL * T_K))
    IC50_uM = Kd_M * 1e6
    log10_IC50_uM = math.log10(IC50_uM)
    return dG_kcal, log10_IC50_uM, IC50_uM


def read_ref2015_scores(csv_path: Path) -> List[Dict[str, Any]]:
    with csv_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        return list(reader)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse RAPiDock ref2015_score.csv into peptide affinity IC50 estimates."
    )
    parser.add_argument("--library", default="data/ligand_library.csv")
    parser.add_argument("--outputs-dir", default="outputs_rapidock")
    parser.add_argument("--results-dir", default="results")
    args = parser.parse_args()

    library = list(csv.DictReader(Path(args.library).open(newline="", encoding="utf-8")))
    outputs_dir = Path(args.outputs_dir)
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    out_csv = results_dir / "ranking_peptides.csv"

    rows: List[Dict[str, Any]] = []
    missing = 0

    for rec in library:
        ligand_type = (rec.get("type") or "").strip().lower()
        if ligand_type != "peptide":
            continue

        name = (rec.get("name") or "").strip()
        if not name:
            continue
        complex_name = name.lower().replace("-", "_").replace(" ", "_")
        out_root = outputs_dir / complex_name
        score_csv = out_root / "ref2015_score.csv"
        run_seconds_path = out_root / "run_seconds.txt"

        if not score_csv.exists():
            missing += 1
            continue

        scores = read_ref2015_scores(score_csv)
        if not scores:
            missing += 1
            continue

        # Expect RAPiDock's ref2015_score.csv to contain a numeric 'ref2015' column
        best_row = None
        best_val = None
        for s in scores:
            try:
                val = float(s.get("ref2015", ""))
            except ValueError:
                continue
            if best_val is None or val < best_val:
                best_val = val
                best_row = s

        if best_row is None or best_val is None:
            missing += 1
            continue

        dG_kcal, log10_ic50_uM, ic50_uM = reu_to_ic50_uM(best_val)

        run_seconds = None
        if run_seconds_path.exists():
            try:
                run_seconds = float(run_seconds_path.read_text(encoding="utf-8").strip())
            except ValueError:
                run_seconds = None

        rows.append(
            {
                "name": name,
                "type": "peptide",
                "source": rec.get("source"),
                "sequence": rec.get("sequence"),
                "sequence_rapidock": rec.get("peptide_seq_rapidock"),
                "rapidock_ref2015_score_REU": best_val,
                "rapidock_dG_kcalmol": dG_kcal,
                "affinity_pred_value_log10IC50uM": log10_ic50_uM,
                "affinity_pred_IC50_uM": ic50_uM,
                "affinity_probability_binary": None,
                "iptm": None,
                "ptm": None,
                "ligand_plddt_mean": None,
                "pocket_plddt_mean": None,
                "interface_pae_mean": None,
                "run_seconds": run_seconds,
                "source_method": "rapidock_ref2015",
                "top_pose_path": str(out_root),
            }
        )

    fieldnames = [
        "name",
        "type",
        "source",
        "sequence",
        "sequence_rapidock",
        "rapidock_ref2015_score_REU",
        "rapidock_dG_kcalmol",
        "affinity_pred_value_log10IC50uM",
        "affinity_pred_IC50_uM",
        "affinity_probability_binary",
        "iptm",
        "ptm",
        "ligand_plddt_mean",
        "pocket_plddt_mean",
        "interface_pae_mean",
        "run_seconds",
        "source_method",
        "top_pose_path",
    ]

    with out_csv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

    print(f"Wrote {len(rows)} peptide affinity rows to {out_csv}")
    if missing:
        print(f"Missing or empty ref2015_score.csv for {missing} peptide entries.")


if __name__ == "__main__":
    main()

