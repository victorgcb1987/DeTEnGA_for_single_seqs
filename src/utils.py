from Bio import SeqIO


def _group_sequences(sequences, out_fpaths):
    with open(out_fpaths["protein"], "w") as protein_fhand:
        for hog, members in sequences.items():
            for member in members:
                record = SeqIO.read(member["protein"], "fasta")
                SeqIO.write(record, protein_fhand, "fasta")
    for mode, out_fpath in out_fpaths.items():
        if "mrna_" in mode:
            kingdom = mode.split("_")
            with open(out_fpath) as out_fhand:
                for hog, members in sequences.items():
                    for member in members:
                        if member["Kingdom"] == mode:
                            record = SeqIO.read(member["mrna"], "fasta")
                            SeqIO.write(record, out_fhand, "fasta")


def generate_input_files(sequences, out_fpath):
    mrna_out_viridiplantae = out_fpath / "mrna_viridiplantae_sequences.fna"
    mrna_out_metazoa = out_fpath / "mrna_metazoa_sequences.fna"
    mrna_out_other = out_fpath / "mrna_other_sequences.fna"
    prot_out = out_fpath / "protein_sequences.faa"
    outs = {"proteins": prot_out,
            "mrna_viridiplantae": mrna_out_viridiplantae,
            "mrna_metazoa": mrna_out_metazoa,
            "mrna_other": mrna_out_other}
    _group_sequences(sequences, outs)

    return outs


def select_longest_isoform(sequence_dir, protein_sequence, mrna_sequence):
    protein_sequences_lengths = []
    records = SeqIO.parse(protein_sequence, "fasta")
    for record in records:
        protein_sequences_lengths.append(len(record.seq))
    if len(protein_sequences_lengths) == 1:
        return protein_sequence, mrna_sequence
    else:
        longest_isoform = protein_sequences_lengths.index(max(protein_sequences_lengths))
        longest_prot_path = sequence_dir / "protein_longest_isoform.faa"
        with open(longest_prot_path, "w") as prot_out_fhand:
            records = [record for record in SeqIO.parse(protein_sequence, "fasta")]
            SeqIO.write(records[longest_isoform], prot_out_fhand, "fasta")
        longest_mrna_path = sequence_dir / "mrna_longest_isoform.fna"
        with open(longest_mrna_path, "w") as mrna_out_fhand:
            records = [record for record in SeqIO.parse(mrna_sequence, "fasta")]
            SeqIO.write(records[longest_isoform], mrna_out_fhand, "fasta")
            mrnaid = records[longest_isoform].id
        return longest_prot_path, longest_mrna_path, mrnaid


def search_sequences(metadata, input_dir):
    found_sequences = {hog: [] for hog in metadata}
    not_found_sequences = found_sequences.copy()
    for hog, members in metadata.items():
        hog_dir = input_dir / hog
        for member in members:
            sequence_dir = hog_dir / member["proteinID"] / "seqs"
            protein_sequence = sequence_dir / "protein.faa"
            mrna_sequence = sequence_dir / "rna.fna"

            if mrna_sequence.is_file() and protein_sequence.is_file():
                protein_sequence, mrna_sequence, mrnaID = select_longest_isoform(sequence_dir, 
                                                                                 protein_sequence, 
                                                                                 mrna_sequence)
                member.update({"protein": protein_sequence,
                               "mrna": mrna_sequence,
                               "main_dir": sequence_dir,
                               "mrnaID": mrnaID})
                found_sequences[hog].append(member)

            else:
                not_found_sequences[hog].append(member)
    return found_sequences, not_found_sequences