import argparse
import csv
import json
import subprocess
import requests
import sys
from pathlib import Path


OUTDIR = Path("feature_tables")
OUTDIR.mkdir(exist_ok=True)



def parse_arguments():
    desc = "Donwload genome feature tables"
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_output_dir = '''(Required) Output dir'''
    parser.add_argument("--output", "-o", type=str,
                        help=help_output_dir,
                        required=True)
    
    help_output_dir = '''(Required) Metadata file'''
    parser.add_argument("--metadata", "-m", type=str,
                        help=help_output_dir,
                        required=True)
    
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    output = Path(parser.output)
    if not output.exists():
        output.mkdir(parents=True)
    return {"out": output,
            "metadata": parser.metadata}





def main():
    args = get_arguments()
    out_dir = args["out"]
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)

    with open(args["metadata"]) as f:
        accessions = set(row["Genome"] for row in csv.DictReader(f, delimiter=","))

    for acc in accessions:

        print(f"Processing {acc}...")

        try:
            # Obtener resumen del ensamblado desde datasets
            cmd = [
                "datasets",
                "summary",
                "genome",
                "accession",
                acc,
                "--as-json"
            ]

            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            data = json.loads(result.stdout)

            report = data["reports"][0]

            ftp_path = report["assembly_info"]["assembly_accession"]

            # Buscar el nombre exacto del ensamblado
            assembly_name = report["assembly_info"]["assembly_name"]

            assembly_dir = f"{acc}_{assembly_name}".replace(" ", "_")

            # Construir URL FTP
            prefix = (
                f"https://ftp.ncbi.nlm.nih.gov/genomes/all/"
                f"{acc[0:3]}/"
                f"{acc[4:7]}/"
                f"{acc[7:10]}/"
                f"{acc[10:13]}/"
                f"{acc.split('.')[0][13:]}/"
                f"{assembly_dir}"
            )

            feature_url = (
                f"{prefix}/{assembly_dir}_feature_table.txt.gz"
            )

            outfile = out_dir / f"{acc}_feature_table.txt.gz"

            print(f"Downloading {feature_url}")

            r = requests.get(feature_url, stream=True, timeout=60)

            if r.status_code == 200:
                with open(outfile, "wb") as fh:
                    for chunk in r.iter_content(chunk_size=8192):
                        fh.write(chunk)

                print(f"Saved: {outfile}")

            else:
                print(f"Not found ({r.status_code})")

        except Exception as e:
            print(f"ERROR: {acc}: {e}")


if __name__ == "__main__":
    main()