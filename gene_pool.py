# Eli Pandolfo
# simulates a GWAS on microbial data

import random
import pickle
from multiprocessing import Process, Manager

class Gene_Pool:

    def __init__(self, num_samples, genome_size, wild_genes, threads):
        
        # passed parameters
        self.NUM_SAMPLES = num_samples
        self.GENOME_SIZE = genome_size
        self.WILD_GENES = wild_genes
        self.NUM_THREADS = threads
        
        # global vars
        self.threads = {}
        self.resistant_snp_loci = None
        self.genes = None
        self.variants = {}
        self.wild_type = None
            

    # create some random loci in the genome, such that at those loci, a SNP of value 'A'
    # confers resistance
    def instantiate_resistant_snp_loci(self):
        self.resistant_snp_loci = [random.randrange(self.GENOME_SIZE)
            for loc in range(int(self.GENOME_SIZE / 1000))]

    # simulate some random genes of roughly accurate nucleotide length for E. coli
    def simulate_genes(self):
        gene_size = self.GENOME_SIZE / 1000
        genes = []
        for i in range(self.WILD_GENES):
            gene = ''
            for j in range(random.randint(int(.7 * gene_size), int(1.3 * gene_size))):
                gene += random.choice('ACTG')
            genes.append(gene)
        self.genes = genes

    # creates a random wild-type genome
    def simulate_wild_type(self):
        self.wild_type = [random.choice('ACTG') for i in range(self.GENOME_SIZE)]

    def thread_simulate_variants(self, dic, tid, n, *args):
        variants = {tid * n + i: Variant(tid * n + i, *args) for i in range(n)}
        dic[tid] = variants

    # simulates an entire gene pool, creating a bunch of sample genomes
    # multithreaded for speed
    def create(self):
        print('Instantiating resistant SNP loci...')
        self.instantiate_resistant_snp_loci()

        print('Creating wild type genome...')
        self.simulate_wild_type()

        print('Creating wild genes...')
        self.simulate_genes()

        print(f'Creating variants using {self.NUM_THREADS} threads...')
        if self.NUM_THREADS == 1:
            self.variants = [Variant(
                                i,
                                self.wild_type,
                                self.GENOME_SIZE,
                                self.genes,
                                self.WILD_GENES,
                                self.resistant_snp_loci
                                ) for i in range(self.NUM_SAMPLES)]

        else:
            dic = Manager().dict()
            for i in range(self.NUM_THREADS):
                self.threads[i] = Process(target=self.thread_simulate_variants, args=(
                    dic,
                    i,
                    int(self.NUM_SAMPLES / self.NUM_THREADS)
                        + (self.NUM_SAMPLES % self.NUM_THREADS if i == self.NUM_THREADS - 1 else 0),
                    self.wild_type,
                    self.GENOME_SIZE,
                    self.genes,
                    self.WILD_GENES,
                    self.resistant_snp_loci))
                self.threads[i].start()

            for i in range(self.NUM_THREADS):
                self.threads[i].join()
                if self.threads[i].is_alive():
                    raise TimeoutError(f'Thread {i} timed out while joining')
                self.variants.update(dic[i])

        print(f'Gene pool created with {len(self.variants)} variants.')

    # take a random sample from the gene pool and output some information about the samples
    def sample(self, n):
        # compare mutants to wild-type
        sample = random.sample(self.variants, n)
        for v in sample:
            print('Resistant:', v.resistance, 
                'SNPs:', v.num_snps, 'Beneficial SNPS:', v.num_beneficial_snps,
                'Genes added:', v.num_genes, 'Beneficial genes added:', v.num_beneficial_genes)


    # simulates a gene pool, serializes it, and writes it to a file
    def dump(self, filename):
        print('Serializing and dumping to file...')
        meta = {
                'NUM_SAMPLES': self.NUM_SAMPLES,
                'GENOME_SIZE': self.GENOME_SIZE,
                'WILD_GENES': self.WILD_GENES,
                'NUM_THREADS': self.NUM_THREADS,
                'resistant_snp_loci': self.resistant_snp_loci,
                'genes': self.genes,
                'wild_type': self.wild_type
            }
        relevant_data = {
            'variants': self.variants,
            'meta': meta,
        }
        with open(filename, 'wb') as gene_pool_file:
            pickle.dump(relevant_data, gene_pool_file)

class Variant:

    P_GENE_ADDITION = .000005
    P_GENE_RESISTANCE = .01
    P_SNP = .0005
    BENEFICIAL_SNP_BASE = 'A'

    # creates a variant from a wild-type genome with some SNPs and possibly some genes inserted
    # there is a small chance that one or more of these mutations confer antibiotic resistance
    def __init__(self, variant_id, wild_type, size, genes, num_wild_genes, resistant_snp_loci):
        self.id = variant_id
        self.size = size
        self.num_snps = 0
        self.num_beneficial_snps = 0
        self.num_genes = 0
        self.num_beneficial_genes = 0
        self.resistance = 0
        self.sequence = ''
        self.seq_list = []

        for i,base in enumerate(wild_type):
            mut_chance = random.randint(0, self.size * 2)
            if mut_chance / (self.size * 2) < self.P_GENE_ADDITION: # .00001 chance of gene addition
                self.num_genes += 1
                gene_added = random.randint(0, num_wild_genes - 1)
                if gene_added < num_wild_genes * self.P_GENE_RESISTANCE: # given gene addition, .01 chance of resistance
                    self.resistance = 1
                    self.num_beneficial_genes += 1
                self.seq_list.append(genes[gene_added])
            elif mut_chance / (self.size * 2) < self.P_SNP: # .0005 chance of snp
                self.num_snps += 1
                mutation_options = 'ACTG'.replace(base, '')
                mutation_choice = random.choice(mutation_options)
                self.seq_list.append(mutation_choice)
                # for simplicity, all beneficial
                # SNPS are 'A', and must occur in one of the resistance snp locations
                if mutation_choice == self.BENEFICIAL_SNP_BASE and i in resistant_snp_loci:
                    self.resistance = 1
                    self.num_beneficial_snps += 1
            else:
                self.seq_list.append(base)

        self.sequence = ''.join(self.seq_list)

    @staticmethod
    def estimate_metrics(size, num_wild_genes):
        print('SNPS: ~', Variant.P_SNP * size)
        print('Beneficial SNPS: ~', round(0.25 * Variant.P_SNP * size / 1000, 3))
        print('Genes added: ~', round(Variant.P_GENE_ADDITION * size, 3))
        print('Beneficial genes added: ~', 
            round(Variant.P_GENE_ADDITION * Variant.P_GENE_RESISTANCE * size, 3))







    