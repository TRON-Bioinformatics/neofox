import tempfile


def create_temp_file(prefix=None, suffix=None, dir=None):
    temp_file = tempfile.NamedTemporaryFile(prefix=prefix, suffix=suffix, dir=dir, delete=False)
    return temp_file.name


def create_temp_fasta(sequences, prefix=None, comment_prefix='seq'):
    """
    Writes seqs given in seqs list into fasta file
    """
    fasta_temp_file = create_temp_file(prefix=prefix, suffix='.fasta')
    counter = 1
    with open(fasta_temp_file, "w") as f:
        for seq in sequences:
            _id = ">{comment_prefix}{index}".format(comment_prefix=comment_prefix, index=counter)
            f.write(_id + "\n")
            f.write(seq + "\n")
            counter += 1
    return fasta_temp_file
