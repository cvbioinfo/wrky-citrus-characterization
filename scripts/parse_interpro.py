#!/usr/bin/env python3
"""Parse InterPro JSON into raw matches and consolidated domain tables."""

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


def protein_id(result):
    xrefs = result.get("xref", [])
    if xrefs:
        return xrefs[0].get("id", xrefs[0].get("name", "UNKNOWN"))
    return "UNKNOWN"


def value_or_na(value):
    return "NA" if value is None else value


def overlap_fraction(a, b):
    overlap = max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1)
    if overlap == 0:
        return 0.0
    len_a = a["end"] - a["start"] + 1
    len_b = b["end"] - b["start"] + 1
    return overlap / min(len_a, len_b)


def cluster_locations(locations, threshold=0.30):
    """Cluster overlapping predictions for one integrated domain accession."""
    parent = list(range(len(locations)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(len(locations)):
        for j in range(i + 1, len(locations)):
            if overlap_fraction(locations[i], locations[j]) >= threshold:
                union(i, j)

    groups = defaultdict(list)
    for i, location in enumerate(locations):
        groups[find(i)].append(location)
    return list(groups.values())


def derive_prefix(input_path):
    name = Path(input_path).name
    for suffix in ("_interpro.json", ".json"):
        if name.lower().endswith(suffix):
            return name[:-len(suffix)]
    return Path(input_path).stem


def main():
    parser = argparse.ArgumentParser(description="Parse InterPro JSON results.")
    parser.add_argument("-i", "--input", required=True, help="InterPro JSON output")
    parser.add_argument("-o", "--output-dir", required=True, help="Directory for parsed results")
    parser.add_argument("--prefix", default=None, help="Output filename prefix; defaults to the JSON filename")
    args = parser.parse_args()

    input_file = Path(args.input)
    output_dir = Path(args.output_dir)
    prefix = args.prefix or derive_prefix(input_file)
    prefix = re.sub(r"[^A-Za-z0-9._-]+", "_", prefix).strip("._-")
    if not prefix:
        parser.error("Output prefix is empty.")
    output_dir.mkdir(parents=True, exist_ok=True)

    with input_file.open(encoding="utf-8") as handle:
        data = json.load(handle)
    results = data["results"]

    metadata_file = output_dir / f"{prefix}_interpro_metadata.txt"
    with metadata_file.open("w", encoding="utf-8") as handle:
        handle.write(f"InterProScan version\t{data.get('interproscan-version', 'NA')}\n")
        handle.write(f"InterPro version\t{data.get('interpro-version', 'NA')}\n")

    raw_rows = []
    domains = defaultdict(lambda: defaultdict(list))
    protein_lengths = {}

    for result in results:
        pid = protein_id(result)
        protein_lengths[pid] = len(result.get("sequence", ""))

        for match in result.get("matches", []):
            source = match.get("source", "NA")
            signature = match.get("signature") or {}
            signature_accession = signature.get("accession", "NA")
            signature_name = signature.get("name")
            signature_description = signature.get("description")
            signature_type = signature.get("type")
            release = signature.get("signatureLibraryRelease") or {}
            database_version = release.get("version", "NA")
            entry = signature.get("entry")

            if entry:
                interpro_accession = entry.get("accession", "NA")
                interpro_name = entry.get("name")
                interpro_description = entry.get("description")
                interpro_type = entry.get("type")
            else:
                interpro_accession = "NA"
                interpro_name = interpro_description = interpro_type = None

            for location in match.get("locations", []):
                start, end = location.get("start"), location.get("end")
                location_evalue = location.get("evalue", match.get("evalue"))
                location_score = location.get("score", match.get("score"))
                representative = location.get("representative", False)

                raw_rows.append({
                    "protein_id": pid,
                    "protein_length": protein_lengths[pid],
                    "source": source,
                    "database_version": database_version,
                    "signature_accession": signature_accession,
                    "signature_name": value_or_na(signature_name),
                    "signature_description": value_or_na(signature_description),
                    "signature_type": value_or_na(signature_type),
                    "interpro_accession": interpro_accession,
                    "interpro_name": value_or_na(interpro_name),
                    "interpro_description": value_or_na(interpro_description),
                    "interpro_type": value_or_na(interpro_type),
                    "start": start,
                    "end": end,
                    "evalue": value_or_na(location_evalue),
                    "score": value_or_na(location_score),
                    "representative": representative,
                })

                # Preserve InterPro-N in the raw table but do not use it as an
                # independent domain prediction for the consolidated table.
                if source == "InterPro-N":
                    continue

                if entry and interpro_type == "Domain":
                    domain_key = interpro_accession
                    domain_description = interpro_description or interpro_name or signature_description or signature_name or domain_key
                    integrated = True
                elif signature_type == "Domain":
                    domain_key = f"{source}:{signature_accession}"
                    domain_description = signature_description or signature_name or signature_accession
                    integrated = False
                else:
                    continue

                domains[pid][domain_key].append({
                    "start": start,
                    "end": end,
                    "source": source,
                    "signature_accession": signature_accession,
                    "domain_description": domain_description,
                    "interpro_accession": interpro_accession,
                    "integrated": integrated,
                })

    raw_file = output_dir / f"{prefix}_interpro_matches.tsv"
    raw_fields = [
        "protein_id", "protein_length", "source", "database_version",
        "signature_accession", "signature_name", "signature_description",
        "signature_type", "interpro_accession", "interpro_name",
        "interpro_description", "interpro_type", "start", "end", "evalue",
        "score", "representative",
    ]
    with raw_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=raw_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(raw_rows)

    domain_rows = []
    for pid in sorted(domains):
        for domain_key, predictions in domains[pid].items():
            clusters = cluster_locations(predictions, threshold=0.30)
            clusters.sort(key=lambda cluster: min(item["start"] for item in cluster))
            for occurrence, cluster in enumerate(clusters, start=1):
                start = min(item["start"] for item in cluster)
                end = max(item["end"] for item in cluster)
                sources = sorted({item["source"] for item in cluster})
                signatures = sorted({f"{item['source']}:{item['signature_accession']}" for item in cluster})
                descriptions = [item["domain_description"] for item in cluster if item["domain_description"]]
                domain_description = descriptions[0] if descriptions else domain_key
                interpro_ids = {item["interpro_accession"] for item in cluster if item["interpro_accession"] != "NA"}
                interpro_accession = sorted(interpro_ids)[0] if interpro_ids else "NA"
                integrated = any(item["integrated"] for item in cluster)
                domain_rows.append({
                    "protein_id": pid,
                    "protein_length": protein_lengths[pid],
                    "domain_occurrence": occurrence,
                    "domain_key": domain_key,
                    "interpro_accession": interpro_accession,
                    "domain_description": domain_description,
                    "domain_span_start": start,
                    "domain_span_end": end,
                    "supporting_databases": ";".join(sources),
                    "supporting_signatures": ";".join(signatures),
                    "n_supporting_signatures": len(signatures),
                    "integrated_interpro": integrated,
                })

    domain_rows.sort(key=lambda row: (row["protein_id"], row["domain_span_start"], row["domain_span_end"]))
    domain_file = output_dir / f"{prefix}_domains.tsv"
    domain_fields = [
        "protein_id", "protein_length", "domain_occurrence", "domain_key",
        "interpro_accession", "domain_description", "domain_span_start",
        "domain_span_end", "supporting_databases", "supporting_signatures",
        "n_supporting_signatures", "integrated_interpro",
    ]
    with domain_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=domain_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(domain_rows)

    summary_file = output_dir / f"{prefix}_domains_summary.tsv"
    domains_by_protein = defaultdict(list)
    for row in domain_rows:
        label = (
            f"{row['domain_description']} ({row['interpro_accession']}) "
            f"[{row['domain_span_start']}-{row['domain_span_end']}]"
        )
        domains_by_protein[row["protein_id"]].append(label)

    # Include proteins with no recognized domains, using the JSON results as
    # the full list of sequence IDs.
    all_proteins = sorted(protein_lengths)
    with summary_file.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["protein_id", "n_domains", "domain_architecture"])
        for pid in all_proteins:
            labels = domains_by_protein.get(pid, [])
            writer.writerow([pid, len(labels), " | ".join(labels) if labels else "None"])

    print(f"Proteins: {len(results)}")
    print(f"Raw matches: {len(raw_rows)}")
    print(f"Consolidated domain occurrences: {len(domain_rows)}")
    print(f"Written: {metadata_file}")
    print(f"Written: {raw_file}")
    print(f"Written: {domain_file}")
    print(f"Written: {summary_file}")


if __name__ == "__main__":
    main()
