def fasta_read(name):
    """
    Read a FASTA file, convert the protein sequence into an HP model, and return the HP sequence.

    Args:
    - name (str): The name of the FASTA file to read.

    Returns:
    - str: The HP sequence representing the protein.

    """
    SEQ = []
    # read the file and select only the sequence
    with open(name, "r") as filout:
        for lines in filout:
            if not lines.startswith(">"):
                SEQ.extend(lines.strip())

    # Change the sequence from amino acid to HP model
    final_seq = []
    for char in SEQ:
        if char in "VIFLMCWGPAC":
            final_seq.append("H")
        else:
            final_seq.append("P")

    # return string of the final sequence in the HP model
    return ''.join(final_seq)
