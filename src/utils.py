from Bio import SeqIO


def _group_sequences(sequences, out_fpaths): 
    fhands = {taxa: open(fpath, "w") for taxa, fpath in out_fpaths.items()}
    for hog, members in sequences.items():
        for member in members:

            fhand_mrna = fhands.get(f'mrna_{member["kingdom"]}', fhands["mrna_other"])
            fhand_protein = fhands.get(f'protein_{member["kingdom"]}', fhands["protein_other"])
            
            mrna_record = SeqIO.read(member["mrna"], "fasta")
            SeqIO.write(mrna_record, fhand_mrna, "fasta")

            protein_record = SeqIO.read(member["protein"], "fasta")
            SeqIO.write(protein_record, fhand_protein, "fasta")

    for kingdom, fhand in fhands.items():
        fhand.flush()
        fhand.close()


def generate_input_files(sequences, out_fpath):
    outs = {"protein_viridiplantae": out_fpath / "protein_viridiplantae_sequences.faa",
            "protein_metazoa": out_fpath / "protein_metazoa_sequences.faa",
            "protein_fungi": out_fpath / "protein_fungi_sequences.faa",
            "protein_other": out_fpath / "protein_other_sequences.faa",
            "mrna_viridiplantae": out_fpath / "mrna_viridiplantae_sequences.fna",
            "mrna_metazoa": out_fpath / "mrna_metazoa_sequences.fna",
            "mrna_fungi": out_fpath / "mrna_fungi_sequences.fna",
            "mrna_other": out_fpath / "mrna_other_sequences.fna"}
    _group_sequences(sequences, outs)

    return outs


def select_isoform(sequence_dir, protein_sequence, 
                           mrna_sequence, protein_id, mrna_id=""):

    protein_records = SeqIO.parse(protein_sequence, "fasta")
    for record in protein_records:
        if record.id == protein_id:
            selected_protein_record = record
            selected_prot_outpath = sequence_dir / "protein_selected_isoform.faa"
            with open(selected_prot_outpath, "w") as prot_out_fhand:
                SeqIO.write(selected_protein_record, prot_out_fhand, "fasta")
                break
    
    mrna_records = SeqIO.parse(mrna_sequence, "fasta")
    selected_mrna_outpath = sequence_dir / "mrna_selected_isoform.fna"
    mrna_found = False
    if mrna_id:
        mrna_id_check = mrna_id.split(".")[0]
        for record in mrna_records:
            if record.id == mrna_id:
                selected_mrna_record = record
                mrna_found = True
                break
    if not mrna_found:
        mrna_lengths = []
        mrna_seqs = []
        #Sometimes, for only one protein can appear non coding transcripts
        #we are going to get the coding ones (starts with XM)
        mrna_records = SeqIO.parse(mrna_sequence, "fasta")
        for record in mrna_records:
            print(record.id)
            if record.id.startswith("XM"):
                mrna_seqs.append(record)
                mrna_lengths.append(len(record.seq))
        longest_idx = mrna_lengths.index(max(mrna_lengths))
        selected_mrna_record = mrna_seqs[longest_idx]

    with open(selected_mrna_outpath, "w") as prot_out_fhand:
        SeqIO.write(selected_mrna_record, prot_out_fhand, "fasta")

    return selected_prot_outpath, selected_mrna_outpath, selected_mrna_record.id


def search_sequences(metadata, input_dir):
    messages = []
    found_sequences = {hog: [] for hog in metadata}
    not_found_sequences = found_sequences.copy()
    for hog, members in metadata.items():
        hog_dir = input_dir / hog
        for member in members:
            print(hog, member)
            sequence_dir = hog_dir / member["proteinID"] / "seqs"
            protein_sequence = sequence_dir / "protein.faa"
            mrna_sequence = sequence_dir / "rna.fna"

            if mrna_sequence.is_file() and protein_sequence.is_file():
                selected_protein_sequence, selected_mrna_sequence,  mrnaID = select_isoform(sequence_dir, 
                                                                                            protein_sequence, 
                                                                                            mrna_sequence,
                                                                                            member["proteinID"],
                                                                                            mrna_id=member["mRNAID"])
               
                member.update({"protein": selected_protein_sequence,
                               "mrna": selected_mrna_sequence,
                               "main_dir": sequence_dir,
                               "mrnaID": mrnaID})
                
                found_sequences[hog].append(member)

            else:
                not_found_sequences[hog].append(member)
    return found_sequences, not_found_sequences, messages