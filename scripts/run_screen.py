#!/usr/bin/env python3
import argparse
import subprocess
import time
from pathlib import Path


def is_completed(output_dir: Path) -> bool:
    pred_dir = output_dir / "predictions"
    return pred_dir.exists() and any(pred_dir.rglob("*.json"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Boltz-2 screen over YAML inputs.")
    parser.add_argument("--inputs-dir", default="inputs")
    parser.add_argument("--outputs-dir", default="outputs")
    parser.add_argument("--recycling-steps", type=int, default=3)
    parser.add_argument("--diffusion-samples", type=int, default=1)
    parser.add_argument("--msa-mode", choices=["server", "none"], default="server")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    inputs_dir = Path(args.inputs_dir)
    outputs_dir = Path(args.outputs_dir)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    input_files = sorted(inputs_dir.glob("*.yaml"))
    if not input_files:
        raise SystemExit(f"No YAML inputs found in {inputs_dir}")

    for input_file in input_files:
        ligand_id = input_file.stem
        out_dir = outputs_dir / ligand_id
        out_dir.mkdir(parents=True, exist_ok=True)

        if not args.force and is_completed(out_dir):
            print(f"[SKIP] {ligand_id}: already completed")
            continue

        cmd = [
            "boltz",
            "predict",
            str(input_file),
            "--out_dir",
            str(out_dir),
            "--recycling_steps",
            str(args.recycling_steps),
            "--diffusion_samples",
            str(args.diffusion_samples),
        ]
        if args.msa_mode == "server":
            cmd.append("--use_msa_server")

        print(f"[RUN ] {ligand_id}")
        start = time.time()
        try:
            subprocess.run(cmd, check=True)
            elapsed = time.time() - start
            print(f"[DONE] {ligand_id} in {elapsed:.1f}s")
            (out_dir / "run_seconds.txt").write_text(f"{elapsed:.3f}\n", encoding="utf-8")
        except subprocess.CalledProcessError as exc:
            print(f"[FAIL] {ligand_id}: return code {exc.returncode}")


if __name__ == "__main__":
    main()
