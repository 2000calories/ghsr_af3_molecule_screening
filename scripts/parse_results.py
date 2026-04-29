#!/usr/bin/env python3
import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd


def maybe_read_json(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def find_first(pattern: str, root: Path) -> Optional[Path]:
    matches = sorted(root.rglob(pattern))
    return matches[0] if matches else None


def mean_or_none(values: List[float]) -> Optional[float]:
    if not values:
        return None
    return float(sum(values) / len(values))


def flatten_numbers(node: Any) -> List[float]:
    out: List[float] = []
    if isinstance(node, (int, float)):
        out.append(float(node))
    elif isinstance(node, list):
        for x in node:
            out.extend(flatten_numbers(x))
    elif isinstance(node, dict):
        for x in node.values():
            out.extend(flatten_numbers(x))
    return out


def parse_confidence(conf: Dict[str, Any]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    iptm = conf.get("iptm")
    ptm = conf.get("ptm")

    pae_fields = [k for k in conf.keys() if "pae" in k.lower()]
    pae_vals: List[float] = []
    for k in pae_fields:
        pae_vals.extend(flatten_numbers(conf.get(k)))
    interface_pae = mean_or_none(pae_vals)
    return iptm, ptm, interface_pae


def parse_affinity(aff: Dict[str, Any]) -> Tuple[Optional[float], Optional[float]]:
    value = aff.get("affinity_pred_value")
    prob = aff.get("affinity_probability_binary")
    return value, prob


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse Boltz-2 outputs into ranking CSV.")
    parser.add_argument("--library", default="data/ligand_library.csv")
    parser.add_argument("--outputs-dir", default="outputs")
    parser.add_argument("--results-dir", default="results")
    args = parser.parse_args()

    library = list(csv.DictReader(Path(args.library).open(newline="", encoding="utf-8")))
    outputs_dir = Path(args.outputs_dir)
    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    rows: List[Dict[str, Any]] = []
    for rec in library:
        name = (rec.get("name") or "").strip()
        if not name:
            continue
        ligand_id = name.lower().replace("-", "_").replace(" ", "_")
        out_root = outputs_dir / ligand_id

        conf_path = find_first("confidence*.json", out_root)
        aff_path = find_first("affinity*.json", out_root)
        run_seconds_path = out_root / "run_seconds.txt"

        conf = maybe_read_json(conf_path) if conf_path else {}
        aff = maybe_read_json(aff_path) if aff_path else {}
        iptm, ptm, interface_pae = parse_confidence(conf)
        affinity_log, affinity_prob = parse_affinity(aff)
        affinity_uM = (10 ** affinity_log) if affinity_log is not None else None

        run_seconds = None
        if run_seconds_path.exists():
            try:
                run_seconds = float(run_seconds_path.read_text(encoding="utf-8").strip())
            except ValueError:
                run_seconds = None

        rows.append(
            {
                "name": name,
                "type": rec.get("type"),
                "source": rec.get("source"),
                "smiles_or_seq": rec.get("smiles") or rec.get("sequence"),
                "affinity_pred_value_log10IC50uM": affinity_log,
                "affinity_pred_IC50_uM": affinity_uM,
                "affinity_probability_binary": affinity_prob,
                "iptm": iptm,
                "ptm": ptm,
                "ligand_plddt_mean": None,
                "pocket_plddt_mean": None,
                "interface_pae_mean": interface_pae,
                "run_seconds": run_seconds,
            }
        )

    df = pd.DataFrame(rows)
    # Affinity lower is better when available; fallback to ipTM (higher better)
    df["sort_affinity"] = df["affinity_pred_value_log10IC50uM"].fillna(math.inf)
    df["sort_iptm"] = df["iptm"].fillna(-1.0)
    df = df.sort_values(by=["sort_affinity", "sort_iptm"], ascending=[True, False]).drop(
        columns=["sort_affinity", "sort_iptm"]
    )
    out_csv = results_dir / "ranking.csv"
    df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv}")


if __name__ == "__main__":
    main()
