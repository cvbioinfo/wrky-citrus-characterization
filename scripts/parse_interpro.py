#!/usr/bin/env python3

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def protein_id(result):
    xrefs = result.get("xref", [])

    if xrefs:
        return xrefs[0].get("id", xrefs[0].get("name", "UNKNOWN"))

    return "UNKNOWN"


def value_or_na(value):
    if value is None:
        return "NA"
    return value


def overlap_fraction(a, b):
    """
    Fraction of the shorter interval that overlaps.
    Coordinates are 1-based and inclusive.
    """

    overlap = max(
        0,
        min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1
    )

    if overlap == 0:
        return 0.0

    len_a = a["end"] - a["start"] + 1
    len_b = b["end"] - b["start"] + 1

    return overlap / min(len_a, len_b)


def cluster_locations(locations, threshold=0.30):
    """
    Group predictions for the same domain that substantially overlap.

    This allows different databases to support the same biological domain
    while preserving separate occurrences of repeated domains.
    """

    n = len(locations)

    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        root_a = find(a)
        root_b = find(b)

        if root_a != root_b:
            parent[root_b] = root_a

    for i in range(n):
        for j in range(i + 1, n):
            if overlap_fraction(locations[i], locations[j]) >= threshold:
                union(i, j)

    groups = defaultdict(list)

    for i, location in enumerate(locations):
        groups[find(i)].append(location)

    return list(groups.values())


def main():
    parser = argparse.ArgumentParser(
        description="Parse InterPro JSON results."
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="InterPro JSON output",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        required=True,
        help="Directory for parsed results",
    )

    args = parser.parse_args()

    input_file = Path(args.input)
    output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    with open(input_file) as handle:
        data = json.load(handle)

    results = data["results"]

    #
    # ------------------------------------------------------------
    # Metadata
    # ------------------------------------------------------------
    #

    metadata_file = output_dir / "RE_interpro_metadata.txt"

    with open(metadata_file, "w") as handle:
        handle.write(
            f"InterProScan version\t"
            f"{data.get('interproscan-version', 'NA')}\n"
        )

        handle.write(
            f"InterPro version\t"
            f"{data.get('interpro-version', 'NA')}\n"
        )

    #
    # ------------------------------------------------------------
    # Raw matches
    # ------------------------------------------------------------
    #

    raw_rows = []

    #
    # Domain candidates used for consolidation:
    #
    # domains[protein][domain_key] = [...]
    #

    domains = defaultdict(lambda: defaultdict(list))

    protein_lengths = {}

    for result in results:

        pid = protein_id(result)

        sequence = result.get("sequence", "")

        protein_lengths[pid] = len(sequence)

        for match in result.get("matches", []):

            source = match.get("source", "NA")

            signature = match.get("signature") or {}

            signature_accession = signature.get(
                "accession", "NA"
            )

            signature_name = signature.get(
                "name"
            )

            signature_description = signature.get(
                "description"
            )

            signature_type = signature.get(
                "type"
            )

            release = signature.get(
                "signatureLibraryRelease"
            ) or {}

            database_version = release.get(
                "version", "NA"
            )

            entry = signature.get("entry")

            if entry:

                interpro_accession = entry.get(
                    "accession", "NA"
                )

                interpro_name = entry.get(
                    "name"
                )

                interpro_description = entry.get(
                    "description"
                )

                interpro_type = entry.get(
                    "type"
                )

            else:

                interpro_accession = "NA"
                interpro_name = None
                interpro_description = None
                interpro_type = None

            for location in match.get(
                "locations", []
            ):

                start = location.get("start")
                end = location.get("end")

                location_evalue = location.get(
                    "evalue",
                    match.get("evalue")
                )

                location_score = location.get(
                    "score",
                    match.get("score")
                )

                representative = location.get(
                    "representative",
                    False
                )

                raw_rows.append(
                    {
                        "protein_id": pid,
                        "protein_length": protein_lengths[pid],
                        "source": source,
                        "database_version": database_version,
                        "signature_accession": signature_accession,
                        "signature_name": value_or_na(
                            signature_name
                        ),
                        "signature_description": value_or_na(
                            signature_description
                        ),
                        "signature_type": value_or_na(
                            signature_type
                        ),
                        "interpro_accession": interpro_accession,
                        "interpro_name": value_or_na(
                            interpro_name
                        ),
                        "interpro_description": value_or_na(
                            interpro_description
                        ),
                        "interpro_type": value_or_na(
                            interpro_type
                        ),
                        "start": start,
                        "end": end,
                        "evalue": value_or_na(
                            location_evalue
                        ),
                        "score": value_or_na(
                            location_score
                        ),
                        "representative": representative,
                    }
                )

                #
                # InterPro-N is retained above but not used
                # for defining consolidated domains.
                #

                if source == "InterPro-N":
                    continue

                #
                # Integrated InterPro Domain
                #

                if (
                    entry
                    and interpro_type == "Domain"
                ):

                    domain_key = interpro_accession

                    domain_description = (
                        interpro_description
                        or interpro_name
                        or signature_description
                        or signature_name
                        or domain_key
                    )

                    integrated = True

                #
                # Domain signature that has not been integrated
                # into an InterPro entry.
                #

                elif signature_type == "Domain":

                    domain_key = (
                        f"{source}:{signature_accession}"
                    )

                    domain_description = (
                        signature_description
                        or signature_name
                        or signature_accession
                    )

                    integrated = False

                else:
                    continue

                domains[pid][domain_key].append(
                    {
                        "start": start,
                        "end": end,
                        "source": source,
                        "signature_accession":
                            signature_accession,
                        "domain_description":
                            domain_description,
                        "interpro_accession":
                            interpro_accession,
                        "integrated": integrated,
                    }
                )

    raw_file = output_dir / "RE_interpro_matches.tsv"

    raw_fields = [
        "protein_id",
        "protein_length",
        "source",
        "database_version",
        "signature_accession",
        "signature_name",
        "signature_description",
        "signature_type",
        "interpro_accession",
        "interpro_name",
        "interpro_description",
        "interpro_type",
        "start",
        "end",
        "evalue",
        "score",
        "representative",
    ]

    with open(
        raw_file,
        "w",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=raw_fields,
            delimiter="\t",
        )

        writer.writeheader()

        writer.writerows(raw_rows)

    #
    # ------------------------------------------------------------
    # Consolidated domains
    # ------------------------------------------------------------
    #

    domain_rows = []

    for pid in sorted(domains):

        for domain_key, predictions in domains[pid].items():

            clusters = cluster_locations(
                predictions,
                threshold=0.30
            )

            clusters.sort(
                key=lambda cluster: min(
                    item["start"]
                    for item in cluster
                )
            )

            for occurrence, cluster in enumerate(
                clusters,
                start=1
            ):

                start = min(
                    item["start"]
                    for item in cluster
                )

                end = max(
                    item["end"]
                    for item in cluster
                )

                sources = sorted(
                    {
                        item["source"]
                        for item in cluster
                    }
                )

                signatures = sorted(
                    {
                        (
                            f"{item['source']}:"
                            f"{item['signature_accession']}"
                        )
                        for item in cluster
                    }
                )

                descriptions = [
                    item["domain_description"]
                    for item in cluster
                    if item["domain_description"]
                ]

                domain_description = (
                    descriptions[0]
                    if descriptions
                    else domain_key
                )

                interpro_ids = {
                    item["interpro_accession"]
                    for item in cluster
                    if item["interpro_accession"]
                    != "NA"
                }

                interpro_accession = (
                    sorted(interpro_ids)[0]
                    if interpro_ids
                    else "NA"
                )

                integrated = any(
                    item["integrated"]
                    for item in cluster
                )

                domain_rows.append(
                    {
                        "protein_id": pid,
                        "protein_length":
                            protein_lengths[pid],
                        "domain_occurrence": occurrence,
                        "domain_key": domain_key,
                        "interpro_accession":
                            interpro_accession,
                        "domain_description":
                            domain_description,
                        "domain_span_start": start,
                        "domain_span_end": end,
                        "supporting_databases":
                            ";".join(sources),
                        "supporting_signatures":
                            ";".join(signatures),
                        "n_supporting_signatures":
                            len(signatures),
                        "integrated_interpro":
                            integrated,
                    }
                )

    domain_rows.sort(
        key=lambda row: (
            row["protein_id"],
            row["domain_span_start"],
            row["domain_span_end"],
        )
    )

    domain_file = output_dir / "RE_domains.tsv"

    domain_fields = [
        "protein_id",
        "protein_length",
        "domain_occurrence",
        "domain_key",
        "interpro_accession",
        "domain_description",
        "domain_span_start",
        "domain_span_end",
        "supporting_databases",
        "supporting_signatures",
        "n_supporting_signatures",
        "integrated_interpro",
    ]

    with open(
        domain_file,
        "w",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=domain_fields,
            delimiter="\t",
        )

        writer.writeheader()

        writer.writerows(domain_rows)

    #
    # ------------------------------------------------------------
    # One-row-per-protein domain summary
    # ------------------------------------------------------------
    #

    domains_by_protein = defaultdict(list)

    for row in domain_rows:

        label = (
            f"{row['domain_description']} "
            f"({row['interpro_accession']}) "
            f"[{row['domain_span_start']}-"
            f"{row['domain_span_end']}]"
        )

        domains_by_protein[
            row["protein_id"]
        ].append(label)

    summary_file = (
        output_dir
        / "RE_domains_summary.tsv"
    )

    with open(
        summary_file,
        "w",
        newline=""
    ) as handle:

        writer = csv.writer(
            handle,
            delimiter="\t"
        )

        writer.writerow(
            [
                "protein_id",
                "n_domains",
                "domain_architecture",
            ]
        )

        for pid in sorted(protein_lengths):

            labels = domains_by_protein.get(
                pid,
                []
            )

            writer.writerow(
                [
                    pid,
                    len(labels),
                    " | ".join(labels)
                    if labels
                    else "None",
                ]
            )

    print(f"Proteins: {len(results)}")
    print(f"Raw matches: {len(raw_rows)}")
    print(
        f"Consolidated domain occurrences: "
        f"{len(domain_rows)}"
    )

    print(f"Written: {raw_file}")
    print(f"Written: {domain_file}")
    print(f"Written: {summary_file}")
    print(f"Written: {metadata_file}")


if __name__ == "__main__":
    main()
