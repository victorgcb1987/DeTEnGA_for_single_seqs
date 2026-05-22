from Bio import SeqIO



def select_longest_isoform(sequence_dir,protein_sequence, mrna_sequence):
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
            records = SeqIO.parse(protein_sequence, "fasta")
            SeqIO.write(records[longest_isoform], prot_out_fhand, "fasta")
        longest_mrna_path = sequence_dir / "mrna_longest_isoform.fna"
        with open(longest_mrna_path, "w") as mrna_out_fhand:
            records = SeqIO.parse(mrna_sequence, "fasta")
            SeqIO.write(records[longest_isoform], mrna_out_fhand, "fasta")
        return longest_prot_path, longest_mrna_path


def search_sequences(metadata, input_dir):
    found_sequences = {hog: [] for hog in metadata}
    not_found_sequences = found_sequences.copy()
    for hog, members in metadata.items():
        hog_dir = input_dir / hog
        for member in members:
            sequence_dir = hog_dir / member["proteinID"] / "seqs"
            protein_sequence = sequence_dir / "protein.faa"
            mrna_sequence = sequence_dir / "mrna.fna"
            if mrna_sequence.is_file() and protein_sequence.is_file():
                protein_sequence, mrna_sequence = select_longest_isoform(sequence_dir, 
                                                                         protein_sequence, 
                                                                         mrna_sequence)
                member.update({"protein": protein_sequence,
                               "mrna": mrna_sequence,
                               "main_dir": sequence_dir})
                found_sequences[hog].append(member)
            else:
                not_found_sequences[hog].append(member)
    return found_sequences, not_found_sequences