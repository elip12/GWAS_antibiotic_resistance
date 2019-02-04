# Eli Pandolfo
# simulates a GWAS on microbial data

import pickle
from collections import OrderedDict

class Analyzer:

    def __init__(self, threads, k, num_samples, gene_pool_file):
        self.THREADS = threads
        self.K = k
        self.NUM_SAMPLES = num_samples
        self.gene_pool_file = gene_pool_file

    # for reusability
    def load(self, filename):
        with open(filename, 'rb') as f:
            raw = pickle.load(f)
        return raw

    # load the gene pool file
    def load_pool(self):
        print('Initializing gene pool...')
        return self.load(self.gene_pool_file)

    # there will be a bunch of kmers frame shifted 1 base from each other, all centered around
    # a single gene or SNP. This coalesces those kmers into a single sequence containing a gene or SNP.
    # Note: this only works if `for key in dict` returns the keys in the order in which they were added.
    def remove_duplicates(self, kmers):
        new_resistant_kmers = {}
        temp = {}
        i = 0

        prev_kmer = 'initial'

        for kmer in kmers: # assumes sorted
            if prev_kmer[1:] != kmer[:-1]: # if this is a totally new kmer
                i += 1
                temp[i] = (kmer, kmers[kmer])
            else: # if this is a previous kmer frame shifted 1
                temp[i] = (temp[i][0] + kmer[-1], temp[i][1])
            prev_kmer = kmer

        for k in temp:
            new_resistant_kmers[temp[k][0]] = temp[k][1]

        return new_resistant_kmers


    # 1. Creates a dict of all kmers in the gene pool, associated with the genomes that
    #    have that kmer
    # 2. Identifies kmers associated with antibiotic resistance by correlating the
    #    known resistance of a genome to that genome having a particular kmer
    # 3. Coalesces adjacent kmers into a longer sequence of nucleotides containing
    #    a SNP or gene
    # 4. Returns a dictionary with kmer strings as keys and lists of genome ids as values
    def find_seqs(self, variants, filename=''):
        kmers = OrderedDict()
        resistant_kmers = {}

        print('\nCreating kmer database...')
        for v_id, variant in variants.items():
            for i in range(len(variant.sequence) - self.K + 1):
                kmer = variant.sequence[i:i + self.K]
                if kmer in kmers:
                    kmers[kmer].append(v_id)
                else:
                    kmers[kmer] = [v_id]

        # consider parallelizing this
        print('\nIdentifying resistant kmers...')
        for kmer in kmers:
            num_resistant = 0
            for v_id in kmers[kmer]:
                if variants[v_id].resistance == 1:
                    num_resistant += 1
            p_resistant = num_resistant / len(kmers[kmer])
            if p_resistant > 0.95 and num_resistant / self.NUM_SAMPLES >= .01:
                resistant_kmers[kmer] = kmers[kmer]
        
        print('\nConsolidating adjacent kmers...')
        resistant_kmers = self.remove_duplicates(resistant_kmers)

        print(f'\n{len(resistant_kmers)} resistant sequences detected')
        if filename:
            print(f'\nDumping sequences to {filename}...')
            with open(filename, 'wb') as dump_file:
                pickle.dump(resistant_kmers, dump_file)
        return resistant_kmers

    def evaluate_accuracy(self, seqs, raw):
        variants = raw['variants']
        meta = raw['meta']
        gene_size = meta['GENOME_SIZE'] / 1000

        identified_genes = {gene: seqs[gene] for gene in seqs if len(gene) > 0.5 * gene_size}
        beneficial_gene = meta['genes'][0] # hacky but i dont want to make new data rn
        # print('\nBeneficial gene:\n', beneficial_gene)

        # print('\nIdentified genes:\n')
        for gene in identified_genes:
            if len(gene) > len(beneficial_gene):
                if beneficial_gene in gene:
                    print('Beneficial gene recognized!')
            else:
                if gene in beneficial_gene:
                    print('Beneficial gene recognized!')

        


        


