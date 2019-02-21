import os
import pickle
from sys import argv

# FASTQ files have metadata on odd lines and sequences on even lines
# ignore all errors
def pickle_input(in_dir):
    seqs = {}
    for dirname, _, filenames in os.walk(in_dir):
        for filename in filenames:
            base, ext = os.path.splitext(filename)
            if ext == '.fa':
                with open(os.path.join(dirname, filename)) as f:
                    contigs = []
                    for line in f:
                        if line[0] == '<':
                            pass
                        elif line[0] in ['A','T','C','G']:
                            contigs.append(line)
                    seqs[base] = contigs


if __name__ == '__main__':
    if len(argv) != 2:
        raise ValueError('Usage: python3 pickle_input.py <directory containing FASTQ files>')
    in_dir = str(argv[1])
    if not os.path.isdir(in_dir):
        raise ValueError('Usage: python3 pickle_input.py <directory containing FASTQ files>')    

    pickle_input(in_dir)
