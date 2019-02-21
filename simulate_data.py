# Eli Pandolfo
# simulates variant E. coli genomes 

import random
import pickle

PICKLE_RAW = 'raw.pickle'
PICKLE_SIM = 'sim.pickle'
PERCENT_RESISTANT = .2

def load_raw():
    with open(PICKLE_RAW, 'rb') as f:
        raw = pickle.load(f)
    return raw

def identify_variant_location(raw):
    # create a database of kmers of length s.t. all samples have exactly 1 of this kmer
    # identify a location relative to the start of that kmer
    # return the relative location in all genomes (returns a list or tuple since location is different
    # in every genome)

    K = 500
    kmers = {}
    count = 0
    print('starting...')
    for raw_id, seq in raw.items():
        if count > 20:
            break
        count += 1
        for contig in seq:
            if len(contig) < K:
                break
            for i in range(len(contig) - K + 1):
                kmer = contig[i:i + K]
                if kmer in kmers:
                    kmers[kmer].append(raw_id)
                else:
                    kmers[kmer] = [raw_id]
    print('next step...')
    for kmer in kmers:
        if len(kmers[kmer]) == count:
            print('kmer identified')

    return (1,1)

def create_resistant_variants(raw, locs, val):
    sims = {}
    snp_choices = 'ACTG'.replace(val, '')
    snp = random.choice(snp_choices)
    
    for genome in raw:
        resistance = 0
        if random.random() < PERCENT_RESISTANT:
            raw[genome][locs[genome][0]][locs[genome[1]]] = snp
            resistance = 1
        sims[genome] = {'seq': raw[genome], 'resistance': resistance}

    return sims

def dump_sims(sims):
    with open(PICKLE_SIM, 'wb') as f:
        pickle.dump(sims, f)

def create_sims():
    raw = load_raw()
    locs, val = identify_variant_location(raw)
    #sims = create_resistant_variants(raw, locs, val)
    #dump_sims(sims)

if __name__ == '__main__':
    create_sims()





