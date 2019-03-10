import multiprocessing as mp, os, pickle, argparse

DEBUG = False
K = 30

### AUXILLIARY FUNCTIONS ###
def printd(*args, **kwargs):
    if DEBUG == True:
        print(*args, **kwargs)

def load_raw():
    printd('Loading data...')
    with open('pickle/sim.pickle', 'rb') as f:
        raw = pickle.load(f)
    return raw

# complements a single base
def complement(base):
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'G':
        return 'C'
    else:
        return 'G'

# breaks input file into chunks to minimize reads
def chunkify(fname,size=512):
    fileEnd = os.path.getsize(fname)
    with open(fname,'rb') as f:
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(size,1)
            f.readline()
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                break

### K-MER DATABASE FUNCTIONS ###

# writes a kmer doct to a file
# all in 1 single write to hopefully minimize data races
def write_kmers(kmers):
    printd('Writing kmers to file...')
    with open('output.txt', 'a+') as f:
        s = '\n'.join([f'{k}{v}' for k,v in kmers.items()]) + '\n'
        f.write(s)

# assigns locations to kmers in kmer dict by iterating through all
# genomes
def process(kmers):
    printd('Extracting kmer locations...')
    for count, (raw_id, raw_dict) in enumerate(raw.items()):
        seq = raw_dict['seq']
        for c_id, contig in enumerate(seq):
            l = len(contig)
            if l >= K: # ensure this contig is long enough to sample
                for i in range(l - K + 1):
                    kmer = contig[i:i + K]
                    if kmer in kmers:
                        kmers[kmer] += f' {raw_id},{c_id},{i}'
        printd(f'\tProcessed genome {count + 1}', end='\r')
    
    kmers = {kmer:val for kmer,val in kmers.items() if (len(val.split(' ')) > 3 and len(val.split(' ')) < 351)}
    write_kmers(kmers)
    printd('Done.')

# creates kmer dict from input chunk
def process_wrapper(chunkStart, chunkSize): 
    kmers = {}
    with open('input.txt') as f:
        f.seek(chunkStart)
        lines = f.read(chunkSize).splitlines()
        for line in lines:
            kmer = line.split(' ')[0]
            comp = ''.join(map(complement, list(kmer)))
            kmers[kmer] = ''
            kmers[comp] = ''
    process(kmers)
    # print(len(kmers))


### Classify by population structure ###

# 1. DSK for kmer input; we are going to use kmers of abundance > 2
def read_full_input(fname):
    kmers = {}
    with open(fname, 'r') as f:
        for line in f:
            kmer = line.split(' ')[0]
            comp = ''.join(map(complement, list(kmer)))
            kmers[kmer] = ''
            kmers[comp] = ''
    return kmers

# 2. Sample 1%

# 3. Distance Matrix
# call distance matrix on this 
# 4. Cluster from matrix




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runs GWAS on microbial genomes')
    parser.add_argument('-d', action='store_true',
        help='turn on debug mode')
    args = parser.parse_args()
    DEBUG = args.d

    #initialize raw data, multiprocessing pool, and jobs queue
    raw = load_raw()
    pool = mp.Pool(4)
    jobs = []

    # create jobs
    # n = 0
    for chunkStart,chunkSize in chunkify('input.txt'):
        # n += 1
        # if n > 4:
            # break
        # printd(f'Starting chunk {n}...')
        jobs.append(pool.apply_async(process_wrapper,(chunkStart,chunkSize)))

    # wait for all jobs to finish
    # n = 0
    for job in jobs:
        job.get()
        # n += 1
        # printd(f'Finished chunk {n}...')

    pool.close()