from Bio import SeqIO


def _group_sequences(sequences, out_fpaths): 
    fhands = {taxa: open(fpath, "w") for taxa, fpath in out_fpaths.items()}
    for hog, members in sequences.items():
        for member in members:
            print(member)
            fhand_mrna = fhands.get(f'mrna_{member["kingdom"]}', fhands["mrna_other"])
            fhand_protein = fhands.get(f'protein{member["kingdom"]}', fhands["protein_other"])
            
            records_mrna = SeqIO.index(member["mrna"], "fasta")
            SeqIO.write(records_mrna[member["mrnaID"]], fhand_mrna, "fasta")

            records_protein = SeqIO.index(member["protein"], "fasta")
            print(member)
            SeqIO.write(records_protein[member["proteinID"]], fhand_protein, "fasta")
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


def select_longest_isoform(sequence_dir, protein_sequence, mrna_sequence):
    protein_sequences_lengths = []
    records = SeqIO.parse(protein_sequence, "fasta")
    for record in records:
        protein_sequences_lengths.append(len(record.seq))
    if len(protein_sequences_lengths) == 1:
        print(protein_sequence, mrna_sequence)
        #Sometimes, for only one protein can appear 2 or more transcripts
        #we are goin to get the first one that is coding (starts with XM)
        records = SeqIO.parse(mrna_sequence, "fasta")
        for record in records:
            if record.id.startswith("XM"):
                mrna_id = record.id
                break
        return protein_sequence, mrna_sequence, mrna_id, ""
    
    else:
        longest_isoform = protein_sequences_lengths.index(max(protein_sequences_lengths))
        longest_prot_path = sequence_dir / "protein_longest_isoform.faa"
        with open(longest_prot_path, "w") as prot_out_fhand:
            records = [record for record in SeqIO.parse(protein_sequence, "fasta")]
            longest_isoform = records[longest_isoform]
            prot_id = longest_isoform.id
            SeqIO.write(longest_isoform, prot_out_fhand, "fasta")
            print(record.description)
            transcript = f'transcript={record.description.split()[-1].split("=")[-1][:-1]}'
        longest_mrna_path = sequence_dir / "mrna_longest_isoform.fna"
        with open(longest_mrna_path, "w") as mrna_out_fhand:
            records = SeqIO.parse(mrna_sequence, "fasta")
            for record in records:
                print(transcript, record.description)
                if transcript in record.description:
                    SeqIO.write(record, mrna_out_fhand, "fasta")
                    mrna_id = record.id
        return longest_prot_path, longest_mrna_path, mrna_id, prot_id


def search_sequences(metadata, input_dir):
    messages = []
    found_sequences = {hog: [] for hog in metadata}
    not_found_sequences = found_sequences.copy()
    for hog, members in metadata.items():
        hog_dir = input_dir / hog
        for member in members:
            sequence_dir = hog_dir / member["proteinID"] / "seqs"
            protein_sequence = sequence_dir / "protein.faa"
            mrna_sequence = sequence_dir / "rna.fna"

            if mrna_sequence.is_file() and protein_sequence.is_file():
                protein_sequence, mrna_sequence, mrnaID, protID = select_longest_isoform(sequence_dir, 
                                                                                 protein_sequence, 
                                                                                 mrna_sequence)
                member.update({"protein": protein_sequence,
                               "mrna": mrna_sequence,
                               "main_dir": sequence_dir,
                               "mrnaID": mrnaID})
                if protID:
                    msg = f'{member["proteinID"]} was changed for {protID} because it was a longer insoform'
                    messages.append(msg)
                    member["proteinID"] = protID
                found_sequences[hog].append(member)

            else:
                not_found_sequences[hog].append(member)
    return found_sequences, not_found_sequences, messages