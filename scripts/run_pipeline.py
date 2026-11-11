#!/usr/bin/env python3
"""Run local WRKY characterization steps using a prefix derived from FASTA name.

External InterPro and WoLF PSORT submissions remain manual. Save their outputs
under results/<prefix>/ using the filenames printed by this script, then rerun.
"""

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def derive_prefix(fasta_path: Path) -> str:
    stem = fasta_path.stem
    suffixes = (
        "_candidate_proteins", "_predicted_proteins", "_protein_sequences",
        "_candidate_sequences", "_proteins", "_protein", "_sequences",
        "_sequence", "_proteome", "_proteins_aa", "_protein_aa",
    )
    for suffix in suffixes:
        if stem.lower().endswith(suffix):
            stem = stem[:-len(suffix)]
            break
    # Ensure the prefix can be used safely in file/directory names.
    prefix = re.sub(r"[^A-Za-z0-9._-]+", "_", stem).strip("._-")
    if not prefix:
        raise ValueError("Could not derive a prefix from the FASTA filename. Use --prefix.")
    return prefix


def resolve_from_root(value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def run(command: list[str]) -> None:
    print("\n$ " + " ".join(str(part) for part in command), flush=True)
    subprocess.run(command, check=True, cwd=ROOT)


def is_stale(output: Path, inputs: list[Path], force: bool) -> bool:
    if force or not output.exists():
        return True
    output_mtime = output.stat().st_mtime
    return any(path.exists() and path.stat().st_mtime > output_mtime for path in inputs)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run WRKY characterization using the FASTA filename as the sample prefix."
    )
    parser.add_argument(
        "--fasta", required=True,
        help="Input protein FASTA, e.g. data/Csinensis_proteins.fasta",
    )
    parser.add_argument(
        "--prefix", default=None,
        help="Optional prefix override. By default it is derived from the FASTA filename.",
    )
    parser.add_argument(
        "--results-root", default="results",
        help="Results directory relative to the project root (default: results).",
    )
    parser.add_argument(
        "--model-dir", default="software/signalp6/signalp6_fast/signalp-6-package/models",
        help="SignalP model directory, relative to project root unless absolute.",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-run local steps even when their outputs appear up to date.",
    )
    args = parser.parse_args()

    fasta = resolve_from_root(args.fasta)
    if not fasta.is_file():
        parser.error(f"FASTA file not found: {fasta}")

    prefix = args.prefix or derive_prefix(fasta)
    prefix = re.sub(r"[^A-Za-z0-9._-]+", "_", prefix).strip("._-")
    if not prefix:
        parser.error("Prefix is empty after sanitization.")

    results_root = resolve_from_root(args.results_root)
    sample_dir = results_root / prefix
    signalp_dir = sample_dir / "signalp"
    interpro_dir = sample_dir / "interpro_parsed"
    signalp_prediction = signalp_dir / "prediction_results.txt"
    base_table = sample_dir / f"{prefix}_results.tsv"
    interpro_json = sample_dir / f"{prefix}_interpro.json"
    interpro_tsv = sample_dir / f"{prefix}_interpro.tsv"
    wolf_file = sample_dir / f"{prefix}_wolfpsort.txt"
    domains_file = interpro_dir / f"{prefix}_domains.tsv"
    final_table = sample_dir / f"{prefix}_final_characterization.tsv"
    model_dir = resolve_from_root(args.model_dir)

    sample_dir.mkdir(parents=True, exist_ok=True)
    signalp_dir.mkdir(parents=True, exist_ok=True)

    # Locate the characterization script; support the spelling already used
    # in this project and the English alternative.
    characterizer_candidates = [
        ROOT / "scripts" / "caracterize_proteins.py",
        ROOT / "scripts" / "characterize_proteins.py",
    ]
    characterizer = next((p for p in characterizer_candidates if p.is_file()), None)
    if characterizer is None:
        parser.error(
            "Could not find scripts/caracterize_proteins.py or "
            "scripts/characterize_proteins.py."
        )

    signalp_executable = shutil.which("signalp6")
    if signalp_executable is None:
        parser.error("signalp6 not found. Activate the environment where SignalP 6 is installed.")
    if not model_dir.is_dir():
        parser.error(f"SignalP model directory not found: {model_dir}. Use --model-dir to set it.")

    print(f"FASTA:       {fasta}")
    print(f"Prefix:      {prefix}")
    print(f"Results dir: {sample_dir}")

    # Step 1: SignalP. Reuse results if they are newer than the FASTA.
    if is_stale(signalp_prediction, [fasta], args.force):
        run([
            signalp_executable,
            "--fastafile", str(fasta),
            "--output_dir", str(signalp_dir),
            "--organism", "eukarya",
            "--mode", "fast",
            "--format", "none",
            "--bsize", "10",
            "--write_procs", "2",
            "--torch_num_threads", "4",
            "--model_dir", str(model_dir),
        ])
    else:
        print(f"\nSignalP result is up to date; reusing {signalp_prediction}")

    if not signalp_prediction.is_file():
        parser.error(f"SignalP finished but prediction file is missing: {signalp_prediction}")

    # Step 2: physicochemical properties + SignalP results.
    if is_stale(base_table, [fasta, signalp_prediction], args.force):
        run([
            sys.executable, str(characterizer),
            "-i", str(fasta),
            "--signalp", str(signalp_prediction),
            "-o", str(base_table),
        ])
    else:
        print(f"\nProperties table is up to date; reusing {base_table}")

    # Step 3: InterPro and WoLF PSORT are web submissions and must be saved by the user.
    missing = []
    if not interpro_json.is_file():
        missing.append(f"InterPro JSON: {interpro_json}")
    if not interpro_tsv.is_file():
        missing.append(f"InterPro TSV (archive/reference): {interpro_tsv}")
    if not wolf_file.is_file():
        missing.append(f"WoLF PSORT text output: {wolf_file}")

    if not interpro_json.is_file() or not wolf_file.is_file():
        print("\nLocal steps completed. The following web results are still needed:")
        for item in missing:
            print(f"  - {item}")
        print("\nInterPro: submit the same FASTA and save its JSON and TSV downloads using the filenames above.")
        print("WoLF PSORT: choose Plant, submit the FASTA, and save the complete text output using the filename above.")
        print("\nAfter saving those files, run this same command again to build the final table.")
        return 0

    # Step 4: parse InterPro JSON and rename outputs using the sample prefix.
    parse_script = ROOT / "scripts" / "parse_interpro.py"
    if not parse_script.is_file():
        parser.error(f"Missing script: {parse_script}")

    if is_stale(domains_file, [interpro_json], args.force):
        interpro_dir.mkdir(parents=True, exist_ok=True)
        run([
            sys.executable, str(parse_script),
            "-i", str(interpro_json),
            "-o", str(interpro_dir),
            "--prefix", prefix,
        ])
    else:
        print(f"\nParsed InterPro domains are up to date; reusing {domains_file}")

    if not domains_file.is_file():
        parser.error(f"InterPro parser did not produce expected file: {domains_file}")

    # Step 5: combine all results.
    build_script = ROOT / "scripts" / "build_final_table.py"
    if not build_script.is_file():
        parser.error(f"Missing script: {build_script}")

    if is_stale(final_table, [base_table, wolf_file, domains_file], args.force):
        run([
            sys.executable, str(build_script),
            "--base", str(base_table),
            "--wolfpsort", str(wolf_file),
            "--domains", str(domains_file),
            "-o", str(final_table),
        ])
    else:
        print(f"\nFinal table is up to date; reusing {final_table}")

    print(f"\nCompleted. Final table: {final_table}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except subprocess.CalledProcessError as exc:
        print(f"\nERROR: command failed with exit code {exc.returncode}.", file=sys.stderr)
        raise SystemExit(exc.returncode)
