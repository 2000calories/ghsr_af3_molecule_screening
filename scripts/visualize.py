#!/usr/bin/env python3
import argparse
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def save_top_pose_files(df: pd.DataFrame, outputs_dir: Path, top_pose_dir: Path) -> None:
    top_pose_dir.mkdir(parents=True, exist_ok=True)
    small = df[df["type"] == "small_molecule"].head(5)
    # Prefer peptides that actually have affinity predictions; otherwise fall back to iptm.
    pep = df[df["type"] == "peptide"].copy()
    if "affinity_pred_value_log10IC50uM" in pep.columns:
        pep = pep.sort_values(
            ["affinity_pred_value_log10IC50uM", "iptm"],
            ascending=[True, False],
        )
    else:
        pep = pep.sort_values("iptm", ascending=False)
    pep = pep.head(5)
    merged = pd.concat([small, pep], ignore_index=True)

    for idx, row in merged.iterrows():
        ligand_id = str(row["name"]).lower().replace("-", "_").replace(" ", "_")
        out_root = outputs_dir / ligand_id
        pdb_candidates = sorted(out_root.rglob("*.pdb"))
        if not pdb_candidates:
            continue
        src = pdb_candidates[0]
        dst = top_pose_dir / f"{idx+1:02d}_{ligand_id}.pdb"
        shutil.copy2(src, dst)


def main() -> None:
    parser = argparse.ArgumentParser(description="Create ranking plots and export top poses.")
    parser.add_argument("--ranking-csv", default="results/ranking.csv")
    parser.add_argument("--outputs-dir", default="outputs")
    parser.add_argument("--plots-dir", default="results/plots")
    parser.add_argument("--top-poses-dir", default="results/top_poses")
    args = parser.parse_args()

    ranking_csv = Path(args.ranking_csv)
    outputs_dir = Path(args.outputs_dir)
    plots_dir = Path(args.plots_dir)
    top_poses_dir = Path(args.top_poses_dir)
    plots_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(ranking_csv)

    aff_df = df[df["affinity_pred_value_log10IC50uM"].notna()].copy()
    aff_df = aff_df.sort_values("affinity_pred_value_log10IC50uM", ascending=True).head(20)
    if not aff_df.empty:
        plt.figure(figsize=(12, 6))
        colors = None
        if "source_method" in aff_df.columns:
            cmap = {
                "boltz2_affinity": "#1f77b4",
                "rapidock_ref2015": "#ff7f0e",
            }
            colors = [cmap.get(m, "#2ca02c") for m in aff_df["source_method"]]
        plt.bar(aff_df["name"], aff_df["affinity_pred_value_log10IC50uM"], color=colors)
        plt.xticks(rotation=65, ha="right")
        plt.ylabel("Predicted log10(IC50 uM)")
        plt.title("Top 20 affinity predictions")
        if "source_method" in aff_df.columns:
            from matplotlib.lines import Line2D

            handles = [
                Line2D([0], [0], color="#1f77b4", lw=8, label="Boltz-2 small molecule"),
                Line2D([0], [0], color="#ff7f0e", lw=8, label="RAPiDock peptide"),
            ]
            plt.legend(handles=handles)
        plt.tight_layout()
        plt.savefig(plots_dir / "top20_affinity.png", dpi=180)
        plt.close()

    scatter_df = df[df["affinity_pred_value_log10IC50uM"].notna() & df["iptm"].notna()]
    if not scatter_df.empty:
        plt.figure(figsize=(8, 6))
        # Group by source method if available, otherwise by type.
        group_key = "source_method" if "source_method" in scatter_df.columns else "type"
        for label, group in scatter_df.groupby(group_key):
            plt.scatter(
                group["iptm"],
                group["affinity_pred_value_log10IC50uM"],
                label=label,
                alpha=0.8,
            )
        plt.xlabel("ipTM")
        plt.ylabel("Predicted log10(IC50 uM)")
        plt.title("ipTM vs predicted affinity")
        plt.legend()
        plt.tight_layout()
        plt.savefig(plots_dir / "iptm_vs_affinity.png", dpi=180)
        plt.close()

    if "type" in df.columns:
        plt.figure(figsize=(6, 4))
        df["type"].value_counts().plot(kind="bar")
        plt.ylabel("Count")
        plt.title("Ligand class distribution")
        plt.tight_layout()
        plt.savefig(plots_dir / "class_distribution.png", dpi=180)
        plt.close()

    save_top_pose_files(df, outputs_dir, top_poses_dir)
    print(f"Wrote plots to {plots_dir} and top poses to {top_poses_dir}")


if __name__ == "__main__":
    main()
