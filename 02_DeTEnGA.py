#!/usr/bin/env python

import argparse
import os
import sys

from pathlib import Path

from src.parsers import (get_pfams_from_db, get_pfams_from_interpro_query, 
                         parse_TEsort_output, classify_pfams, classify_protein, write_summary,
                         get_stats, read_metadata)
from src.run import run_TEsorter, remove_stop_codons, run_interpro, run_agat

from src.utils import search_sequences, generate_input_files

REXDB_PFAMS = {"rexdb-plant": Path(os.path.dirname(os.path.realpath(__file__))) / "data" / "Viridiplantae_2.0_pfams.txt",
               "rexdb-metazoa": Path(os.path.dirname(os.path.realpath(__file__))) / "data" / "Metazoa_3.1_pfams.txt",
               "rexdb": Path(os.path.dirname(os.path.realpath(__file__))) / "data" / "Combined_pfams.txt"}


CATEGORIES = {"No_TE(PcpM0)": "PcpM0", "Protein_TE_only(PteM0)": "PteM0",
              "Chimeric_Protein_Only(PchM0)": "PchM0", "mRNA_TE_Only(PcpMte)": "PcpMte",
              "Protein_and_mRNA_TE(PteMte)": "PteMte", "Chimeric_Protein_and_mRNA_TE(PchMte)": "PchMte",
              "No_Protein_Domains_mRNA_TE(P0Mte)": "P0Mte"}

#Generating program options
def parse_arguments():
    desc = "DeTEnGA fork for analyze TEs in HOGs"
    parser = argparse.ArgumentParser(description=desc)
    

    help_input = '''(Required) Input dir'''
    parser.add_argument("--input", "-i", type=str,
                        help=help_input,
                        required=True)
    
    help_output_dir = '''(Required) Output dir'''
    parser.add_argument("--output", "-o", type=str,
                        help=help_output_dir,
                        required=True)
    
    help_output_dir = '''(Required) Metadata file'''
    parser.add_argument("--metadata", "-m", type=str,
                        help=help_output_dir,
                        required=True)
    
    help_threads = "(Optional) number of threads. 1 by default"
    parser.add_argument("--threads", "-t", type=int,
                        help=help_threads, default=1)
    
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


def get_arguments():
    parser = parse_arguments()
    output = Path(parser.output)
    if not output.exists():
        output.mkdir(parents=True)
    return {"input": Path(parser.input),
            "out": output,
            "threads": parser.threads,
            "metadata": parser.metadata}


def create_header():
    header = ["Run", "Genome", "Annotation", "Annotated_transcripts"]
    for key in CATEGORIES:
        header.append(f"{key}_N")
    for key in CATEGORIES:
        header.append(f"{key}_%")
    header += ["Summary_N", "Summary_%"]
    return "\t".join(header)+"\n"


def get_row(label, genome, annotation, stats):
    inverse_categories = {value: key for key, value in CATEGORIES.items()}
    categories = ["T"] + [key for key in inverse_categories]
    values = [str(stats["num_transcripts"])] + [str(stats[key]) for key in inverse_categories]
    per_values = [str(stats["num_transcripts"])] + [str(round(float(stats[key]/stats["num_transcripts"])*100, 2)) for key in inverse_categories]
    summary = "{0}: {1};{2}: {3};{4}: {5};{6}: {7};{8}: {9};{10}: {11};{12}: {13}"
    row = [label, genome, annotation]
    row += values
    row += per_values[1:]
    row += [summary.format(*[item for pair in zip(categories, values) for item in pair])]
    row += [summary.format(*[item for pair in zip(categories, per_values) for item in pair])]    
    return "\t".join(row)+"\n"


def emit_message(msg, fhand):
    print(msg)
    fhand.write(msg+"\n")
    fhand.flush()

def main():
    args = get_arguments()
    out_dir = args["out"]
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
    #Create log file
    log = out_dir / "DeTEnGA_log.txt"
    log_fhand = open(log, "a")
    log_fhand.write("#Command used: {}\n".format(" ".join(sys.argv)))

    with open(args["metadata"]) as metadata_fhand:
        metadata = read_metadata(metadata_fhand)

    msg = "Checking if all sequences are avaiable"
    emit_message(msg, log_fhand)
    found_sequences, not_found_sequences = search_sequences(metadata, args["input"])

    input_fpaths = generate_input_files(found_sequences, out_dir)


    msg = "Analising RNA transposable elements with TEsorter"
    emit_message(msg, log_fhand)
    for kind, fpath in input_fpaths.items():
        if kind.startswith("mrna_"):
            kingdom = kind.split("_")[-1]
            if kingdom == "viridiplantae":
                rex_db = "rexdb-plant"
            elif kingdom == "metazoa":
                rex_db = "rexdb-metazoa"
            else:
                rex_db = "rex-db"
            TEsorter_results = run_TEsorter(fpath, rex_db, args["threads"])
            if TEsorter_results["returncode"] == 99:
                msg = f'TEsorter already done for {kingdom}. Skipping Tesorter analysis'
            elif TEsorter_results["returncode"] == 0:
                msg = f'TEsorter succesfully run for {kingdom}.'
            else:
                msg = f'TEsorter failed for protein {kingdom}: {TEsorter_results["msg"]}'
            
            emit_message(msg, log_fhand)
            

    msg = "Truncating protein sequence at first stop codon"
    emit_message(msg, log_fhand)
    stop_codons_out = remove_stop_codons(input_fpaths["protein"])

       

    #Run interproscan
    msg = "Analyze protein domains with interproscan"
    emit_message(msg, log_fhand)
    interpro_results = run_interpro(stop_codons_out, args["threads"])
    if interpro_results["returncode"] == 99:
        msg = 'InterproScan already done. Skipping Interpro analysis'
    elif interpro_results["returncode"] == 0:
        msg = 'InteproScan succesfully run'
    else:
        msg = f'InteproScan failed: {interpro_results["msg"]}'
    emit_message(msg, log_fhand)
    database = REXDB_PFAMS[database]
    TE_pfams = get_pfams_from_db(database)

    #         with open(TEsorter_results["out_fpath"]) as TEsorter_fhand:
    #             te_sorter_output = parse_TEsort_output(TEsorter_fhand)
        
    #         with open(interpro_results["out_fpath"]) as interpro_fhand:
    #             interpro = get_pfams_from_interpro_query(interpro_fhand)
    #             classified_pfams = classify_pfams(interpro, TE_pfams)
    
    #         protein_class = classify_protein(classified_pfams, te_sorter_output)
            
    #         out_fpath = Path(out_dir / label / "{}_TE_summary.csv".format(label))
    #         with open(out_fpath, "w") as out_fhand:
    #             write_summary(te_summary, out_fhand)
    #             summaries[label] = out_fpath
    #             msg = "TE Summary for {} written in {}\n".format(label, out_fpath)
    #             log_fhand.write(msg)
    #             log_fhand.flush()
    
    # msg = "##STEP 6: Running stats on annotation files\n"
    # print(msg)
    # log_fhand.write(msg)
    # log_fhand.flush()
    # agat_results = run_agat(summaries, files)
    # with open(args["out"]/ "combined_summaries.tsv", "w") as combined_summaries_fhand:
    #     header = create_header()
    #     combined_summaries_fhand.write(header)
    #     for label, results in agat_results.items():
    #         stats = get_stats(results["out_fpath"], summaries[label])
    #         genome = files[label]["assembly"].stem
    #         annotation = files[label]["annotation"].stem
    #         row = get_row(label, genome, annotation, stats)
    #         combined_summaries_fhand.write(row)


if __name__ == "__main__":
    main()