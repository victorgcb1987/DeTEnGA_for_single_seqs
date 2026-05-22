import os


from pathlib import Path
from subprocess import run

def run_gffread(fof, output):
    results_catalog = {}
    for label, values in fof.items():
        out_dir = output / label
        if not out_dir.exists():
            out_dir.mkdir(parents=True, exist_ok=True)
        prefix  = values["assembly"].stem
        mrna_out = out_dir / "{}.mRNA.fasta".format(prefix)
        pep_out = out_dir / "{}.pep.fasta".format(prefix)
        
        
        cmd = "gffread -w {} -g {} {}".format(str(mrna_out), 
                                              str(values["assembly"]),
                                              str(values["annotation"]))
        if mrna_out.exists():
            returncode = 99
            msg = "File {} already exists\n".format(str(mrna_out))
        else:
            run_ = run(cmd, capture_output=True, shell=True)
            msg = run_.stderr.decode()
            returncode = run_.returncode
        
        results = {"command": {"mrna": cmd}, "returncode": {"mrna": returncode},
                   "msg": {"mrna": msg}, "out_fpath": {"mrna": mrna_out.absolute()}}
        
        cmd = "gffread -y {} -g {} {}".format(str(pep_out), 
                                              str(values["assembly"]),
                                              str(values["annotation"]))
        if pep_out.exists():
            returncode = 99
            msg = "File {} already exists\n".format(str(pep_out))
        else:
            run_ = run(cmd, capture_output=True, shell=True)
            msg = run_.stderr.decode()
            returncode = run_.returncode

        results["command"]["protein"] = cmd
        results["returncode"]["protein"] = returncode
        results["msg"]["protein"] = msg
        results["out_fpath"]["protein"] = pep_out.absolute()
        results_catalog[label] = results

    return results_catalog


def run_TEsorter(input_mrna, database, threads):
    base_dir = Path(os.getcwd())
    results = {}
    out_mrna = Path("{}.{}.cls.tsv".format(input_mrna, database))
    os.chdir(out_mrna.parents[0].absolute())
    if database == "rex-db":
        #it will run an combination of viridiplantae and Metazoa
        cmd = f'TEsorter {input_mrna.name} -p {threads}'
    else:
        cmd = f'TEsorter {input_mrna.name} -db {database} -p {threads}'
    print(out_mrna)
    if Path(out_mrna.name).exists():
        returncode = 99
        msg = "File {} already exists\n".format(str(out_mrna))
    else:
        run_ = run(cmd, capture_output=True, shell=True)
        returncode = run_.returncode
        if returncode == 0:
            msg = "Done, check TEsorter.log.txt for details \n"
        else:
            msg = run_.stderr.decode()
        with open("TEsorter.log.txt", "w") as log_fhand:
            log_fhand.write(run_.stderr.decode())
    results = {"command": cmd, "returncode": returncode,
               "msg": msg, "out_fpath": out_mrna}
    os.chdir(base_dir)
    return results



def remove_stop_codons(protein_sequence):
    out_fpath = Path("{}/{}.nostop.fasta".format(protein_sequence.parents[0], protein_sequence.stem))
    if out_fpath.exists():
        return {"returncode": 99,
                "msg": "protein sequence already truncated",
                "out_fpath": out_fpath}

    else:
        stop_codons = [".", "*"]
        original_sequence = ""
        truncated_sequence = ""
        with open(protein_sequence) as protein_fhand:
            for line in protein_fhand:
                if line.startswith(">"):
                    header = line
                    protein_name = header.replace(">", "").split()
                else:
                    original_sequence += line.rstrip()
                    truncated_sequence += line.rstrip()

                for stop_codon in stop_codons:
                    if stop_codon in truncated_sequence:
                        truncated_sequence = truncated_sequence.split(stop_codon)[0]
                        with open(out_fpath, "w") as out_fhand:
                            out_fhand.write(header)
                            for line in (truncated_sequence[0+i:60+i] for i in range(0, len(truncated_sequence), 60)):
                                out_fhand.write(line+"\n")
                        msg = f'truncated, original length {len(original_sequence)}, truncated length {len(truncated_sequence)}'
                        return {"msg": msg, "returncode": 0, 
                                "out_fpath": out_fpath}
                return {"msg": f'no stop codons found in protein', 
                        "returncode": 0, 
                        "out_fpath": protein_sequence}


def run_interpro(sequence, threads):
    exclude = ["AntiFam", "CDD", "Coils", "FunFam",
               "Gene3D", "Hamap", "MobiDBLite",
               "NCBIfam", "PANTHER", "PIRSF", 
               "PIRSR", "PRINTS", "ProSitePatterns",
               "ProSiteProfiles", "SFLD", "SMART", 
               "SUPERFAMILY"]
    base_dir = Path(os.getcwd())
    
    os.chdir(sequence.parents[0].absolute())
    out_fpath = Path("{}.tsv".format(sequence))
    log_fpath = Path("interpro.log.txt")
    cmd = f'interproscan.sh -i {sequence.name} -cpu {threads} '
    cmd += f'-exclappl {",".join(exclude)} --disable-precalc > {log_fpath}'
    if Path(out_fpath.name).exists():
            returncode = 99
            msg = "InterProScan analysis already done, Skipping"
    else:
        run_ = run(cmd, shell=True, capture_output=True)
        returncode = run_.returncode
        if returncode == 0:
            msg = "InterProScan ran succesfully"
        else:
            msg = run_.stderr.decode() + " " + run_.stdout.decode()
    interpro_results = {"command": cmd, "msg": msg,
                        "out_fpath": out_fpath, "returncode": returncode}
        
    os.chdir(base_dir)
    return interpro_results
    
        
def run_agat(summaries, annotations):
    agat_results = {}
    for label, summary in summaries.items():
        base_dir = summary.parents[0].absolute()
        agat_out = base_dir / "{}.agat.stats.txt".format(label)
        annot_file = annotations[label]["annotation"]
        cmd = "agat_sp_statistics.pl --gff {} -o {}".format(str(annot_file), 
                                                            agat_out)
        if agat_out.exists():
            returncode = 99
            msg = "File {} already exists".format(str(agat_out))
        else:
            run_ = run(cmd, shell=True, capture_output=True)
            returncode = run_.returncode
            if returncode == 0:
                msg = "Done"
            else:
                msg = run_.stdout.decode()
        agat_results[label] = {"command": cmd, "msg": msg,
                               "out_fpath": agat_out, "returncode": returncode}
    return agat_results 