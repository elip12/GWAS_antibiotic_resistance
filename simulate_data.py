# Eli Pandolfo
# simulates variant E. coli genomes 

import random
import pickle
import argparse

PICKLE_RAW = 'raw.pickle'
PICKLE_SIM = 'sim.pickle'
PERCENT_RESISTANT = .2
DEBUG = False

def printd(*args):
    if DEBUG == True:
        print(*args)

def load_raw():
    printd('Loading raw data...')
    with open(PICKLE_RAW, 'rb') as f:
        raw = pickle.load(f)
    return raw


# 1. Create a database of kmers such that:
#   - almost every genome has this exact kmer
#   - almost no genomes have more than 1 of this exact kmer
# 2. Sample any kmer from that database (we use the first)
# 3. Identify the value of a sampled location within that kmer
# 4. Return the relative location of the SNP within every genome, and the value to be changed
def identify_variant_location(raw):

    K = 50
    kmers = {
        # of the form:
        # 'ATGCCCCTC...': {'id_0': (0,234), 'id_1': (5,6625)}
        # where keys are the kmer sequence itself,
        # and values are a dict with
        # keys being the id of a genome,
        # and values being a tuple with
        # the contig index this kmer is present in,
        # and the string index within that contig where this kmer starts
    }
    printd('Creating kmer database...')
    for count, (raw_id, seq) in enumerate(raw.items()):
        for c_id, contig in enumerate(seq):
            l = len(contig)
            if l >= K: # ensure this contig is long enough to sample
                for i in range(l - K + 1):
                    kmer = contig[i:i + K] # sample the kmer itself
                    # store the id, contig number, and string index
                    if kmer in kmers:
                        kmers[kmer][raw_id] = (c_id, i)
                    elif count == 0:
                        kmers[kmer] = {raw_id: (c_id, i)}
        if (count < 10 and count > 0) or count % 15 == 1:
            # filter kmers to reduce mem usage
            kmers = {kmer:ids for kmer, ids in kmers.items() if len(ids) >= count * .83}
            printd('Usable kmers identified in first', count + 1, 'genomes:\t', len(kmers))
            if len(kmers) == 0:
                exit()
    
    kmers = {kmer:ids for kmer, ids in kmers.items() if len(ids) >= count * .83}
    printd('Usable kmers identified:', len(kmers))

    # pick first (arbitrary) kmer
    kmer, locs = next(iter(kmers.items()))
    for raw_id, tup in locs.items(): # iterate thru all ids in dict
        # store a unique index in every genome for the location of a SNP (to be added)
        locs[raw_id] = (tup[0], tup[1] + int(K / 2))
    
    val = kmer[int(K / 2)] # extract the base that will be changed to create the SNP
    return locs, val

# Adds a SNP to a random sample of genomes.
# The SNP should be at the same relative location in each genome
def create_resistant_variants(raw, locs, val):
    printd('Creating resistant variants...')
    sims = {}
    snp_choices = 'ACTG'.replace(val, '')
    snp = random.choice(snp_choices) # random base change to create the SNP
    
    for raw_id, seq in raw.items():
        resistance = 0
        if random.random() < PERCENT_RESISTANT and raw_id in locs:
            loc = locs[raw_id]
            contig = seq[loc[0]]
            seq[loc[0]] = contig[:loc[1]] + snp + contig[loc[1] + 1:] # add the SNP to this genome 
            resistance = 1
        sims[raw_id] = {'seq': seq, 'resistance': resistance}

    return sims

def dump_sims(sims):
    printd('Dumping sims data...')
    with open(PICKLE_SIM, 'wb') as f:
        pickle.dump(sims, f)

def create_sims():
    raw = load_raw()
    locs, val = identify_variant_location(raw)
    sims = create_resistant_variants(raw, locs, val)
    dump_sims(sims)
    printd('Done.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Simulates SNP in microbial genomes')
    parser.add_argument('-d', action='store_true',
        help='turn on debug mode')
    args = parser.parse_args()
    DEBUG = args.d
    create_sims()





