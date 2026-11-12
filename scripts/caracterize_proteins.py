#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis


VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_signalp(signalp_file):
    """
    Read SignalP 6 prediction_results.txt.

    Returns a dictionary indexed by protein ID.
    """

    results = {}

    if signalp_file is None:
        return results

    signalp_file = Path(signalp_file)

    if not signalp_file.exists():
        raise FileNotFoundError(
            f"SignalP result file not found: {signalp_file}"
        )

    with open(signalp_file) as handle:
        for line in handle:

            if line.startswith("#") or not line.strip():
                continue

            fields = [
    		field.strip()
    		for field in line.rstrip("\r\n").split("\t")
	    ]

            protein_id = fields[0]
            prediction = fields[1]
            other_probability = float(fields[2])
            sp_probability = float(fields[3])

            cleavage_position = "NA"
            cleavage_probability = "NA"

            if len(fields) > 4 and fields[4].strip():

                cs = fields[4]

                position_match = re.search(
                    r"CS pos:\s*(\d+-\d+)",
                    cs
                )

                probability_match = re.search(
                    r"Pr:\s*([0-9.]+)",
                    cs
                )

                if position_match:
                    cleavage_position = position_match.group(1)

                if probability_match:
                    cleavage_probability = float(
                        probability_match.group(1)
                    )

            results[protein_id] = {
                "signalp_prediction": prediction,
                "signal_peptide": "Yes" if prediction == "SP" else "No",
                "signalp_other_probability": other_probability,
                "signalp_sp_probability": sp_probability,
                "signalp_cleavage_site": cleavage_position,
                "signalp_cleavage_probability": cleavage_probability,
            }

    return results


def analyze_protein(record, signalp_results):

    sequence = str(record.seq).upper().rstrip("*")

    length = len(sequence)

    invalid = sorted(
        set(sequence) - VALID_AA
    )

    result = {
        "protein_id": record.id,
        "description": record.description,
        "length_aa": length,
        "molecular_weight_da": "NA",
        "molecular_weight_kda": "NA",
        "pI": "NA",
        "signal_peptide": "NA",
        "signalp_prediction": "NA",
        "signalp_other_probability": "NA",
        "signalp_sp_probability": "NA",
        "signalp_cleavage_site": "NA",
        "signalp_cleavage_probability": "NA",
        "status": "OK",
    }

    if invalid:

        result["status"] = (
            f"Ambiguous residues: {','.join(invalid)}"
        )

    else:

        analysis = ProteinAnalysis(sequence)

        molecular_weight = analysis.molecular_weight()
        pi = analysis.isoelectric_point()

        result["molecular_weight_da"] = round(
            molecular_weight, 2
        )

        result["molecular_weight_kda"] = round(
            molecular_weight / 1000, 3
        )

        result["pI"] = round(pi, 2)

    if record.id in signalp_results:

        result.update(
            signalp_results[record.id]
        )

    return result


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Calculate physicochemical properties and "
            "integrate SignalP predictions."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input protein FASTA file",
    )

    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output TSV file",
    )

    parser.add_argument(
        "--signalp",
        help="SignalP 6 prediction_results.txt",
    )

    args = parser.parse_args()

    input_file = Path(args.input)
    output_file = Path(args.output)

    if not input_file.exists():
        raise FileNotFoundError(
            f"Input file not found: {input_file}"
        )

    output_file.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    signalp_results = read_signalp(
        args.signalp
    )

    records = list(
        SeqIO.parse(input_file, "fasta")
    )

    if not records:
        raise ValueError(
            "No protein sequences were found in the FASTA file."
        )

    fieldnames = [
        "protein_id",
        "description",
        "length_aa",
        "molecular_weight_da",
        "molecular_weight_kda",
        "pI",
        "signal_peptide",
        "signalp_prediction",
        "signalp_other_probability",
        "signalp_sp_probability",
        "signalp_cleavage_site",
        "signalp_cleavage_probability",
        "status",
    ]

    with open(
        output_file,
        "w",
        newline=""
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
        )

        writer.writeheader()

        for record in records:

            result = analyze_protein(
                record,
                signalp_results
            )

            writer.writerow(result)

    print(
        f"Proteins analyzed: {len(records)}"
    )

    print(
        f"SignalP predictions loaded: "
        f"{len(signalp_results)}"
    )

    print(
        f"Results written to: {output_file}"
    )


if __name__ == "__main__":
    main()
