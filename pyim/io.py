import skbio


def read_fasta(file_path):
    with file_path.open('r') as file_:
        for seq in skbio.io.read(file_, format='fasta'):
            yield seq
