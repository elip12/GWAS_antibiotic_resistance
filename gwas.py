# Eli Pandolfo

import numpy as np
import random
import pickle

K = 30
NUM_SAMPLES = 300
GENOME_SIZE = 3000000
WILD_GENES = 10

RESISTANT_SNP_LOCS = [random.randrange(GENOME_SIZE) for loc in range(int(GENOME_SIZE / 100))]

def simulate_genes(size):
    gene_size = size / 1000
    genes = []
    for i in range(WILD_GENES): # 10 genes exist in gene pool that mutants can acquire
        gene = ''
        for j in range(random.randint(int(.7 * gene_size), int(1.3 * gene_size))):
            gene += random.choice('ACTG')
        genes.append(gene)
    return genes

def simulate_genome(size):
    s = ''
    for i in range(size):
        s += random.choice('ACTG')
    return s

def simulate_variant(wild_type, size, genes):
    mutant = ''
    num_snps = 0
    num_beneficial_snps = 0
    num_genes = 0
    num_beneficial_genes = 0
    resistance = 0

    for i,base in enumerate(wild_type):
        mut_chance = random.randint(0, size * 2)
        if mut_chance / (size * 2) > .9999: # .0001 chance of gene addition
            num_genes += 1
            gene_added = random.randint(0, WILD_GENES - 1)
            if gene_added <= WILD_GENES / 5:
                resistance = 1
                num_beneficial_genes += 1
            mutant += genes[gene_added]    
        elif mut_chance / (size * 2) > .99: # 1 percent chance of snp
            num_snps += 1
            mutation_options = 'ACTG'.replace(base, '')
            mutation_choice = random.choice(mutation_options)
            mutant += mutation_choice
            if mutation_choice == 'A' and i in RESISTANT_SNP_LOCS:
                resistance = 1
                num_beneficial_snps += 1
        else:
            mutant += base

    return (mutant, resistance, num_snps, num_beneficial_snps, num_genes, num_beneficial_genes)

def simulate_gene_pool(num, size):
    wild_type = simulate_genome(size)
    genes = simulate_genes(size)
    return [simulate_variant(wild_type, size, genes) for i in range(num)]

def sample_variants(variants, n):
    # compare mutants to wild-type
    sample = random.sample(variants, n)
    for genome in sample:
        print('SNPs: ', genome[2], 'Beneficial SNPS:', genome[3],
            'Genes added:', genome[4], 'Beneficial genes added:', genome[5])

def remove_duplicates(kmers):
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

def create_gene_pool():
    raw = simulate_gene_pool(NUM_SAMPLES, GENOME_SIZE)

    with open('gene_pool', 'wb') as gene_pool_file:
        pickle.dump(raw, gene_pool_file)

def identify_resistant_kmers():
    print('Initializing gene pool...')
    with open('gene_pool', 'rb') as gene_pool_file:
        raw = pickle.load(gene_pool_file)

    print('\nSampling variants...')
    sample_variants(raw, 10)
    
    kmers = {}
    resistant_kmers = {}

    print('\nCreating kmer database...')
    for index, genome in enumerate(raw):
        for i in range(len(genome[0]) - K + 1):
            kmer = genome[0][i:i+K]
            if kmer in kmers:
                kmers[kmer].append(index)
            else:
                kmers[kmer] = [index]

    print('\nIdentifying resistant kmers...')
    for kmer in kmers:
        num_resistant = 0
        for genome in kmers[kmer]:
            if raw[genome][1] == 1:
                num_resistant += 1
        p_resistant = num_resistant / len(kmers[kmer])
        if p_resistant > 0.95 and num_resistant / NUM_SAMPLES >= .01:
            resistant_kmers[kmer] = kmers[kmer]
            # print(kmer, kmers[kmer], len(kmers[kmer]), p_resistant)
    
    print('\nRemoving duplicates...')
    resistant_kmers = remove_duplicates(resistant_kmers)

    print('\nResistant kmers:', len(resistant_kmers))
    for kmer in resistant_kmers:
        print(kmer, resistant_kmers[kmer])


create_gene_pool()
#identify_resistant_kmers()
