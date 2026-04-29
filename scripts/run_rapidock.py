#!/usr/bin/env python3
import argparse
import subprocess
import time
from pathlib import Path


def is_completed(output_dir: Path) -> bool:
    # RAPiDock stores scores in ref2015_score.csv when scoring_function=ref2015
    return (output_dir / "ref2015_score.csv").exists()


def main() -> None:
    parser = argparse.ArgumentParser(description="Run RAPiDock over peptide inputs.")
    parser.add_argument("--inputs-csv", default="inputs_rapidock/protein_peptide.csv")
    parser.add_argument("--outputs-dir", default="outputs_rapidock")
    parser.add_argument(
        "--rapidock-dir",
        required=True,
        help="Path to RAPiDock repo root (must contain inference.py and train_models/).",
    )
    parser.add_argument("--ckpt", default="rapidock_local.pt", help="Checkpoint filename in train_models directory.")
    parser.add_argument("--N", type=int, default=10, help="Number of samples per complex.")
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--cpu", type=int, default=8)
    parser.add_argument("--force", action="store_true", help="Re-run even if ref2015_score.csv already exists.")
    args = parser.parse_args()

    inputs_csv = Path(args.inputs_csv)
    outputs_dir = Path(args.outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    rapidock_dir = Path(args.rapidock_dir)
    inference_py = rapidock_dir / "inference.py"
    if not inference_py.exists():
        raise SystemExit(f"Could not find RAPiDock inference.py at {inference_py}")

    # We invoke RAPiDock in multi-task CSV mode; it manages per-complex outputs.
    cmd = [
        "python",
        str(inference_py),
        "--protein_peptide_csv",
        str(inputs_csv.resolve()),
        "--output_dir",
        str(outputs_dir.resolve()),
        "--model_dir",
        str((rapidock_dir / "train_models" / "CGTensorProductEquivariantModel").resolve()),
        "--ckpt",
        args.ckpt,
        "--scoring_function",
        "ref2015",
        "--N",
        str(args.N),
        "--batch_size",
        str(args.batch_size),
        "--no_final_step_noise",
        "--inference_steps",
        "16",
        "--actual_steps",
        "16",
        "--conformation_partial",
        "1:1:1",
        "--cpu",
        str(args.cpu),
    ]

    # RAPiDock writes all complexes in one go, so we only have a coarse-grained timing here.
    out_run_seconds = outputs_dir / "run_seconds.txt"
    start = time.time()
    print("[RUN ] RAPiDock peptide docking")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        print(f"[FAIL] RAPiDock run exited with code {exc.returncode}")
        raise
    elapsed = time.time() - start
    out_run_seconds.write_text(f"{elapsed:.3f}\n", encoding="utf-8")
    print(f"[DONE] RAPiDock in {elapsed:.1f}s, wrote outputs to {outputs_dir}")


if __name__ == "__main__":
    main()

