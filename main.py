# Eli Pandolfo

from gene_pool import Gene_Pool, Variant
from gwas import Analyzer

# number of threads to use when simulating data and analyzing kmers
THREADS = 4

# kmer size
K = 30

# number of genomes to analyze
NUM_SAMPLES = 300

# number of base pairs per genome
GENOME_SIZE = 3000000

# the number of predefined genes that can be added to a genome during mutation simulation
WILD_GENES = 100

PICKLE_FILE = 'gene_pool.pickle'
SEQS_FILE = 'identified_sequences.pickle'

def instantiate_pool(filename):
    pool = Gene_Pool(NUM_SAMPLES, GENOME_SIZE, WILD_GENES, THREADS)
    pool.create()
    pool.dump(filename)

if __name__ == '__main__':

    instantiate_pool(PICKLE_FILE)

    #analyzer = Analyzer(THREADS, K, NUM_SAMPLES, PICKLE_FILE)
    #raw = analyzer.load_pool()
    #seqs = analyzer.find_seqs(raw['variants'], SEQS_FILE)
    #seqs = analyzer.load(SEQS_FILE)
    #analyzer.evaluate_accuracy(seqs, raw)
