#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path
from typing import Any, Dict, List


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
    parser.add_argument(
        "--vina-csv",
        default=None,
        help="Optional path to ranking_peptides_vina.csv. Defaults to <results-dir>/ranking_peptides_vina.csv.",
    )
    args = parser.parse_args()

    library = list(csv.DictReader(Path(args.library).open(newline="", encoding="utf-8")))
    outputs_dir = Path(args.outputs_dir)
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    out_csv = results_dir / "ranking_peptides.csv"
    vina_csv = Path(args.vina_csv) if args.vina_csv else (results_dir / "ranking_peptides_vina.csv")

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

        # RAPiDock outputs can use either 'ref2015' or 'ref2015score' headers.
        best_row = None
        best_val = None
        for s in scores:
            try:
                raw_val = s.get("ref2015", "")
                if raw_val in (None, ""):
                    raw_val = s.get("ref2015score", "")
                val = float(raw_val)
            except ValueError:
                continue
            if best_val is None or val < best_val:
                best_val = val
                best_row = s

        if best_row is None or best_val is None:
            missing += 1
            continue

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
                "rapidock_ref2015_rank": None,
                "rapidock_ref2015_delta_from_best_REU": None,
                # Keep these columns for downstream schema compatibility; disabled by design.
                "rapidock_dG_kcalmol": None,
                "affinity_pred_value_log10IC50uM": None,
                "affinity_pred_IC50_uM": None,
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

    # REF2015 is used as a relative ranking signal only.
    rows.sort(key=lambda r: float(r["rapidock_ref2015_score_REU"]))
    if rows:
        best_score = float(rows[0]["rapidock_ref2015_score_REU"])
        for idx, row in enumerate(rows, start=1):
            score = float(row["rapidock_ref2015_score_REU"])
            row["rapidock_ref2015_rank"] = idx
            row["rapidock_ref2015_delta_from_best_REU"] = score - best_score

    # Optional merge: Vina peptide rescoring rows keyed by peptide name.
    vina_by_name: Dict[str, Dict[str, Any]] = {}
    if vina_csv.exists():
        with vina_csv.open(newline="", encoding="utf-8") as fh:
            for vr in csv.DictReader(fh):
                key = (vr.get("name") or "").strip()
                if key:
                    vina_by_name[key] = vr

        for row in rows:
            vr = vina_by_name.get((row.get("name") or "").strip())
            if not vr:
                continue
            row["vina_best_kcalmol"] = vr.get("vina_best_kcalmol")
            row["vina_delta_from_best_kcalmol"] = vr.get("vina_delta_from_best_kcalmol")
            row["vina_rank"] = vr.get("vina_rank")
            row["vina_source_pose"] = vr.get("vina_source_pose")

    fieldnames = [
        "name",
        "type",
        "source",
        "sequence",
        "sequence_rapidock",
        "rapidock_ref2015_score_REU",
        "rapidock_ref2015_rank",
        "rapidock_ref2015_delta_from_best_REU",
        "vina_best_kcalmol",
        "vina_delta_from_best_kcalmol",
        "vina_rank",
        "vina_source_pose",
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
    if vina_csv.exists():
        print(f"Merged Vina rows from {vina_csv}")
    else:
        print(f"Vina CSV not found (skip merge): {vina_csv}")
    if missing:
        print(f"Missing or empty ref2015_score.csv for {missing} peptide entries.")


if __name__ == "__main__":
    main()

