#!/usr/bin/env python3

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


def read_base_table(path):
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)

    return rows, reader.fieldnames


def read_wolfpsort(path):
    results = {}

    pattern = re.compile(
        r"^(\S+)\s+details\s+(.+)$"
    )

    with open(path) as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            match = pattern.match(line)

            if not match:
                continue

            protein_id = match.group(1)
            predictions_text = match.group(2)

            predictions = []

            for item in predictions_text.split(","):
                item = item.strip()

                if ":" not in item:
                    continue

                location, score = item.split(":", 1)

                try:
                    score = float(score.strip())
                except ValueError:
                    continue

                predictions.append(
                    (location.strip(), score)
                )

            if not predictions:
                continue

            predictions.sort(
                key=lambda x: x[1],
                reverse=True
            )

            primary_location, primary_score = predictions[0]

            results[protein_id] = {
                "wolfpsort_primary": primary_location,
                "wolfpsort_primary_score": primary_score,
                "wolfpsort_all": ";".join(
                    f"{loc}:{score:g}"
                    for loc, score in predictions
                ),
            }

    return results


def read_domains(path):
    domains = defaultdict(list)

    with open(path) as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t"
        )

        for row in reader:
            domains[row["protein_id"]].append(row)

    return domains


def summarize_domains(domain_rows):
    wrky = []
    others = []

    for row in domain_rows:

        description = row["domain_description"]
        accession = row["interpro_accession"]

        start = row["domain_span_start"]
        end = row["domain_span_end"]

        label = (
            f"{description} "
            f"({accession}) "
            f"[{start}-{end}]"
        )

        if accession == "IPR003657":
            wrky.append(row)
        else:
            others.append(label)

    wrky_coordinates = ";".join(
        f"{row['domain_span_start']}-"
        f"{row['domain_span_end']}"
        for row in wrky
    )

    wrky_support = ";".join(
        (
            f"{row['domain_span_start']}-"
            f"{row['domain_span_end']}:"
            f"{row['supporting_databases']}"
        )
        for row in wrky
    )

    if not wrky:
        wrky_status = "not_detected"

    elif any(
        int(row["n_supporting_signatures"]) == 1
        for row in wrky
    ):
        wrky_status = "review"

    else:
        wrky_status = "supported"

    all_domains = []

    for row in domain_rows:
        all_domains.append(
            (
                f"{row['domain_description']} "
                f"({row['interpro_accession']}) "
                f"[{row['domain_span_start']}-"
                f"{row['domain_span_end']}]"
            )
        )

    return {
        "wrky_domain_count": len(wrky),
        "wrky_coordinates":
            wrky_coordinates if wrky else "NA",
        "wrky_support":
            wrky_support if wrky else "NA",
        "wrky_status": wrky_status,
        "other_domains":
            " | ".join(others) if others else "None",
        "interpro_domains":
            " | ".join(all_domains)
            if all_domains else "None",
    }


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--base",
        required=True,
        help="RE_results.tsv"
    )

    parser.add_argument(
        "--wolfpsort",
        required=True,
        help="WoLF PSORT text output"
    )

    parser.add_argument(
        "--domains",
        required=True,
        help="RE_domains.tsv"
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True
    )

    args = parser.parse_args()

    base_rows, base_fields = read_base_table(
        args.base
    )

    wolf = read_wolfpsort(
        args.wolfpsort
    )

    domains = read_domains(
        args.domains
    )

    extra_fields = [
        "wolfpsort_primary",
        "wolfpsort_primary_score",
        "wolfpsort_all",
        "wrky_domain_count",
        "wrky_coordinates",
        "wrky_support",
        "wrky_status",
        "other_domains",
        "interpro_domains",
    ]

    output_file = Path(args.output)

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    with open(
        output_file,
        "w",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=base_fields + extra_fields,
            delimiter="\t",
        )

        writer.writeheader()

        for row in base_rows:

            pid = row["protein_id"]

            wolf_result = wolf.get(
                pid,
                {
                    "wolfpsort_primary": "NA",
                    "wolfpsort_primary_score": "NA",
                    "wolfpsort_all": "NA",
                }
            )

            domain_summary = summarize_domains(
                domains.get(pid, [])
            )

            row.update(wolf_result)
            row.update(domain_summary)

            writer.writerow(row)

    print(
        f"Final table written to: {output_file}"
    )

    print(
        f"Proteins written: {len(base_rows)}"
    )


if __name__ == "__main__":
    main()
