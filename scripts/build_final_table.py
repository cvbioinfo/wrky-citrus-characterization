#!/usr/bin/env python3
"""Combine physicochemical, SignalP, WoLF PSORT and InterPro domain data."""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def read_base_table(path):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    # A technical sequence-validation status is not a biological feature.
    if "status" in fields:
        fields.remove("status")
        for row in rows:
            row.pop("status", None)
    return rows, fields


def read_wolfpsort(path):
    results = {}
    pattern = re.compile(r"^(\S+)\s+details\s+(.+)$")
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            match = pattern.match(line)
            if not match:
                continue
            pid, prediction_text = match.group(1), match.group(2)
            predictions = []
            for item in prediction_text.split(","):
                item = item.strip()
                if ":" not in item:
                    continue
                location, score = item.split(":", 1)
                try:
                    score = float(score.strip())
                except ValueError:
                    continue
                predictions.append((location.strip(), score))
            if not predictions:
                continue
            predictions.sort(key=lambda item: item[1], reverse=True)
            primary, score = predictions[0]
            results[pid] = {
                "wolfpsort_primary": primary,
                "wolfpsort_primary_score": score,
                "wolfpsort_all": ";".join(f"{loc}:{value:g}" for loc, value in predictions),
            }
    return results


def read_domains(path):
    domains = defaultdict(list)
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            domains[row["protein_id"]].append(row)
    return domains


def summarize_domains(domain_rows):
    wrky = [row for row in domain_rows if row["interpro_accession"] == "IPR003657"]
    others = []
    all_domains = []
    for row in domain_rows:
        label = (
            f"{row['domain_description']} ({row['interpro_accession']}) "
            f"[{row['domain_span_start']}-{row['domain_span_end']}]"
        )
        all_domains.append(label)
        if row["interpro_accession"] != "IPR003657":
            others.append(label)

    wrky_coordinates = ";".join(
        f"{row['domain_span_start']}-{row['domain_span_end']}" for row in wrky
    )
    wrky_support = ";".join(
        f"{row['domain_span_start']}-{row['domain_span_end']}:{row['supporting_databases']}"
        for row in wrky
    )
    if not wrky:
        status = "not_detected"
    elif any(int(row["n_supporting_signatures"]) == 1 for row in wrky):
        status = "review"
    else:
        status = "supported"

    return {
        "wrky_domain_count": len(wrky),
        "wrky_coordinates": wrky_coordinates or "NA",
        "wrky_support": wrky_support or "NA",
        "wrky_status": status,
        "other_domains": " | ".join(others) if others else "None",
        "interpro_domains": " | ".join(all_domains) if all_domains else "None",
    }


def main():
    parser = argparse.ArgumentParser(description="Build the final WRKY characterization table.")
    parser.add_argument("--base", required=True, help="Physicochemical/SignalP TSV")
    parser.add_argument("--wolfpsort", required=True, help="WoLF PSORT text output")
    parser.add_argument("--domains", required=True, help="Consolidated InterPro domains TSV")
    parser.add_argument("-o", "--output", required=True, help="Final TSV path")
    args = parser.parse_args()

    base_rows, base_fields = read_base_table(args.base)
    wolf = read_wolfpsort(args.wolfpsort)
    domains = read_domains(args.domains)
    extra_fields = [
        "wolfpsort_primary", "wolfpsort_primary_score", "wolfpsort_all",
        "wrky_domain_count", "wrky_coordinates", "wrky_support", "wrky_status",
        "other_domains", "interpro_domains",
    ]

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=base_fields + extra_fields, delimiter="\t")
        writer.writeheader()
        for row in base_rows:
            pid = row["protein_id"]
            row.update(wolf.get(pid, {
                "wolfpsort_primary": "NA",
                "wolfpsort_primary_score": "NA",
                "wolfpsort_all": "NA",
            }))
            row.update(summarize_domains(domains.get(pid, [])))
            writer.writerow(row)

    print(f"Final table written to: {output}")
    print(f"Proteins written: {len(base_rows)}")


if __name__ == "__main__":
    main()
