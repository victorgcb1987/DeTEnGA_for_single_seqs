import os
import re

from collections import defaultdict
from csv import DictReader
from pathlib import Path


def read_metadata(metadata_fhand):
    metadata = {}
    for row in DictReader(metadata_fhand, delimiter=","):
        hog = row["HOG"]
        prot_metadata = {"proteinID": row["Protein"],
                         "species": row["SpName"],
                         "kingdom": row["Kingdom"].lower(),
                         "category": row["Category"],
                         "mRNAID": row["mRNA"]}
        if not hog in metadata:
            metadata[row["HOG"]] = [prot_metadata]
        else:
            metadata[hog].append(prot_metadata)
    return metadata


def get_pfams_from_db(fpath):
    pfams = {}
    with open(fpath) as fhand:
        for line in fhand:
            if line.startswith("#"):
                continue
            line = line.split()
            pfams[line[0]] = " ".join(line[1:])
    return pfams


def get_pfams_from_interpro_query(fhand):
    genes = defaultdict(list)
    for line in fhand:
        line = line.split("\t")
        if line[3] == "Pfam":
            gen, code, description, start, end = line[0], line[4], line[5], line[6], line[7]
            genes[gen].append([code, description, start, end])
    sorted_genes = {key: sorted(value, key=lambda x: int(x[2]))  # Ordena por el tercer valor (convertido a entero)
                    for key, value in genes.items()}
    return sorted_genes


def classify_pfams(interpro, te_pfams):
    for gene, pfams in interpro.items():
        for pfam in pfams:
            if pfam[0] in te_pfams:
                pfam.append("TE")
            else:
                pfam.append("NT")
    return interpro


def parse_TEsort_output(fhand):
    output = defaultdict(list)
    for line in DictReader(fhand,delimiter="\t"):
        output[line["#TE"]] = {"domains": line["Domains"], 
                               "complete": line["Complete"],
                               "classification": "{}|{}|{}".format(line["Order"],
                                                                   line["Superfamily"],
                                                                   line["Clade"]),
                               "strand": line["Strand"]}
    return output


def classify_protein(interpro_classified, tesort_output, equivalences):
    summary = []
    for protein, mrna in equivalences.items():
        row = {}
        status = ""
        pfams_ids = ""
        pfams_descriptions = ""
        transposable = False
        no_transposable = False
        
        row = {"ProtID": protein,
               "mRNAID": mrna}
        if protein not in interpro_classified:
            status = "NA"
        else:
            pfams_ids = []
            pfams_descriptions = []
            for value in interpro_classified[protein]:
                pfams_ids.append(value[0])
                pfams_descriptions.append(value[1])
                if "TE" in value:
                    transposable = True
                if "NT" in value:
                    no_transposable = True
            if transposable and not no_transposable:
                status = "transposable_element"
            if not transposable and no_transposable:
                status = "coding_sequence"
            if transposable and no_transposable:
                status = "mixed"

        transcript_tesort = tesort_output.get(mrna, None)
        if transcript_tesort is not None:
            row["tesort_domains"] = transcript_tesort["domains"]
            row["tesort_complete"] = transcript_tesort["complete"]
            row["tesort_class"] = transcript_tesort["classification"]
            row["tesort_strand"] = transcript_tesort["strand"]
        else:
            row["tesort_domains"] = "NA"
            row["tesort_complete"] = "NA"
            row["tesort_class"] = "NA"
            row["tesort_strand"] = "NA"

        row["interpro_status"] = status
        row["pfams_ids"] = "|".join(pfams_ids)
        if not pfams_descriptions:
            row["pfams_descriptions"] = "NA" 
        else:
            row["pfams_descriptions"] = "|".join(pfams_descriptions)     
        row["detenga_status"] = detenga_status(row)
        summary.append(row)
    return summary


def detenga_status(row):
    status = "NA"
    if row["interpro_status"] == "coding_sequence" and row["tesort_domains"] == "NA":
        status = "PcpM0" 
    if row["interpro_status"] == "transposable_element" and row["tesort_domains"] == "NA":
        status = "PteM0"
    if row["interpro_status"] == "mixed" and row["tesort_domains"] == "NA":
        status = "PchM0" 
    if row["interpro_status"] == "coding_sequence" and row["tesort_domains"] != "NA":
        status = "PcpMte"
    if row["interpro_status"] == "transposable_element" and row["tesort_domains"] != "NA":
        status = "PteMte"
    if row["interpro_status"] == "mixed" and row["tesort_domains"] != "NA":
        status = "PchMte"
    if row["interpro_status"] == "NA" and row["tesort_domains"] != "NA":
        status = "P0Mte"
    if row["interpro_status"] == "NA" and row["tesort_domains"] == "NA":
        status = "P0M0"
    return status


def write_summary(summary, out_fhand, hogs):
    out_fhand.write("ProteinID,mRNAID,HOG,Interpro_status,TEsort_class;PFAM_domains,")
    out_fhand.write("PFAM_descriptions,TEsort_domains,TEsort_completness,")
    out_fhand.write("TEsort_strand,DeTEnGA_status\n")
    out_fhand.flush()
    for row in summary:
        line_total = "" 
        line_total += f'{row["ProtID"]},{row["mRNAID"]},'
        line_total += f'{hogs[row["ProtID"]]},'
        line_total += f'{row["interpro_status"]},{row["tesort_class"]},'
        line_total += f'{row["pfams_ids"]},{row["pfams_descriptions"].replace(";", " ").replace(",", " ")};'
        line_total += f'{row["tesort_domains"]},{row["tesort_complete"]},'
        line_total += f'{row["tesort_strand"]},{row["detenga_status"]}\n'
        out_fhand.write(line_total)
        out_fhand.flush()


def get_stats(agat_stats, summary):
    with open(agat_stats) as agat_fhand:
        text = agat_fhand.read()
        try:
            match = re.search(r"Number of mrna\s+(\d+)", text, re.IGNORECASE)
            num_transcripts = int(match.group(1))
        except:
            match = re.search(r"Number of transcript\s+(\d+)", text, re.IGNORECASE)
            num_transcripts = int(match.group(1))
    stats = {"PcpM0": 0, "PteM0": 0, "PchM0": 0, 
             "PcpMte": 0, "PteMte": 0, "PchMte": 0, 
             "P0Mte":0, "P0M0":0, "num_transcripts": num_transcripts}
    for row in DictReader(open(summary), delimiter=";"):
        stats[row["DeTEnGA_status"]] += 1
    return stats
