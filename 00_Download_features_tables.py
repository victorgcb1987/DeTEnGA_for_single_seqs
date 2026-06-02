import argparse
import csv
import gzip
import re
import requests
import sys

from subprocess import run
from pathlib import Path



def parse_arguments():
    desc = "Generate NCBI feature tables for multispecies"
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_output_dir = '''(Required) Output dir'''
    parser.add_argument("--output", "-o", type=str,
                        help=help_output_dir,
                        required=True)
    
    help_metadata_file = '''(Required) Metadata file'''
    parser.add_argument("--metadata", "-m", type=str,
                        help=help_metadata_file,
                        required=True)
    
    help_genbank_file = '''(Required) Genbank feature file'''
    parser.add_argument("--genbank", "-b", type=str,
                        help=help_genbank_file,
                        required=True)
    
    help_refseq_file = '''(Required) Refseq feature file'''
    parser.add_argument("--refseq", "-r", type=str,
                        help=help_refseq_file,
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
            "metadata": parser.metadata,
            "genbank": parser.genbank,
            "refseq": parser.refseq}

def get_sequences_ids_by_accession(metadata_file):
    seqsID = {}
    with open(metadata_file) as metadata_fhand:
        for row in csv.DictReader(metadata_fhand, delimiter=","):
            accession = row["Genome"]
            protein = row["Protein"]
            if accession not in seqsID:
                seqsID[accession] = [protein]
            else:
                seqsID[accession].append(protein)
    return seqsID


def get_ftp_link_by_accession(file, filter=[]):
    ftp_urls = {}
    with open(file) as fhand:
        for row in csv.DictReader(fhand, delimiter="\t"):
            accession = row["assembly_accession"]
            if accession not in filter:
                continue
            ftp = row["ftp_path"]
            suffix = ftp.split("/")[-2]
            feature_link = f'{ftp}{suffix}_feature_table.txt.gz'
            ftp_urls[accession] = Path(feature_link)
    return ftp_urls


def download_feature_tables(ftp_links, out_dir):
        downloaded_files = {}
        feature_tables = out_dir / "feature_tables"
        if not feature_tables.exists():
            feature_tables.mkdir(parents=True)
        with open(out_dir / "download_log.txt", "w") as log_fhand:
            for accession, ftp_link in ftp_links.items():
                out_file = feature_tables / ftp_link.name
                if not out_file.is_file():
                    cmd = f'curl {ftp_link} -o {out_file}'
                    cmd = run(cmd, shell=True, capture_output=True)
                    returncode = cmd.returncode
                    if returncode == 0:
                        msg = f'{accession} {ftp_link} downloaded successfully\n'
                    else:
                        msg = f'{accession} {ftp_link} download failed {cmd.stderr}\n'
                        returncode = 1
                else:
                    msg = f'{accession} {ftp_link} download done already \n'
                    returncode = 0
                log_fhand.write(msg)
                downloaded_files[accession] = {"file": out_file, 
                                               "returncode": returncode}  

        return downloaded_files


def retrive_equivalence_info(feature_table):
    equivalences = {}
    with gzip.open(feature_table, "rt") as feat_fhand:
        for line in feat_fhand:
            if line.startswith("#") or not line:
                continue
            else:
                line = line.rstrip().split("\t")
                feat = line[0]
                if feat == "mRNA":
                    mrnaID = line[10]
                    proteinID = line[11]
                    equivalences[proteinID] = mrnaID
    return equivalences


def get_seqs_equivalences(feature_tables, seqsIDs_by_accession):
    equivalences = {}
    for accession, seqIDs in seqsIDs_by_accession.items():
        try: 
            feature_table_file = feature_tables[accession]
        except:
            print(f'accession {accession} not found')
            continue
        if feature_table_file["returncode"] == 0:
            equivalence = retrive_equivalence_info(feature_table_file["file"])
            equivalences[accession] = equivalence
    return equivalences


def get_feature_table_url(accession, parent_url):
    html = requests.get(parent_url).text

    m = re.search(
        rf'({re.escape(accession)}[^"/]*)/',
        html
    )

    if not m:
        return None

    dirname = m.group(1)

    return (
        f"{parent_url}/{dirname}/"
        f"{dirname}_feature_table.txt.gz"
    )

def find_suppressed_accessions(seqIDs_by_accession, downloaded_files, out_dir):

    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/"
    for accession in seqIDs_by_accession:
        if accession not in downloaded_files:
            letter_code = accession[0:3]
            first_part = accession[3:6]
            second_part = accession[6:9]
            third_part = accession[9:]
            url_part = f'{letter_code}/{first_part}/{second_part}/{third_part}/'
            parent_url = url+url_part
            print(parent_url)
            table_url = get_feature_table_url(accession, parent_url)
            print(table_url)



def main():
    args = get_arguments()
    print("#1: Getting accessions from metadata")
    seqIDs_by_accession = get_sequences_ids_by_accession(args["metadata"])
    filter = [accession for accession in seqIDs_by_accession]
    print("#2: getting ftp urls from genbanks")
    #accession_genbank_ftp = get_ftp_link_by_accession(args["genbank"], filter=filter)
    print("#3: gettingftp urls from refseqs")
    accession_refseq_ftp = get_ftp_link_by_accession(args["refseq"], filter=filter)
    print("#4: merging data")
    #merged_ftp_links = accession_genbank_ftp | accession_refseq_ftp
    print("#5: downloading feature tables")
    downloaded_files = download_feature_tables(accession_refseq_ftp, args["out"])
    print(f'Found ftp URLs for {len(downloaded_files)} of {len(seqIDs_by_accession)}')
    print("#6 Trying to reconstruct ftp URLs from supressed accessions")
    suppressed_accessions_links = find_suppressed_accessions(seqIDs_by_accession, downloaded_files, args["out"])
    print("#7: get protein-mrna equivalence")
    equivalences = get_seqs_equivalences(downloaded_files, seqIDs_by_accession)
    print(equivalences)

if __name__ == "__main__":
    main()
