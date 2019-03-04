# Eli Pandolfo
# simulates a GWAS on microbial data

'''
The kmer dict gets too large to store in memory. We need to break it up into chunks of ~8gb,
pickle those chunks, and write them to disc. Then, when identifying the resistant k-mers,
we read one chunk at a time, we extract a single kmer's info from each one, and do the
association test (X^2 test) once we get all the info for a single kmer.

maybe worth using a real database?
'''

import pickle
import argparse
import io
import multiprocessing as mp

PICKLE_RAW = 'pickle/sim.pickle' #'pickle/raw.pickle'
DSK_OUTPUT = 'output.txt'
DATABASE_OUTPUT = 'database.txt'
DEBUG = False
K = 30

def printd(*args):
    if DEBUG == True:
        print(*args)

# load the gene pool file
def load_raw():
    printd('Loading data...')
    with open(PICKLE_RAW, 'rb') as f:
        raw = pickle.load(f)
    return raw

def read_infile(infile):
    printd('Reading DSK output...')
    kmers = []
    with open(infile, 'r') as f:
        for line in f:
            kmer = line.split(' ')[0]
            kmers.append(kmer)
    return kmers

def search_for_kmer(raw, search_kmer):
    tups = []
    print('start')
    for raw_id, raw_dict in raw.items():
        seq = raw_dict['seq']
        for c_id, contig in enumerate(seq):
            l = len(contig)
            if l >= K: # ensure this contig is long enough to sample
                for i in range(l - K + 1):
                    kmer = contig[i:i + K]
                    if kmer == search_kmer:
                        tups.append(str((raw_id, str(c_id), str(i))))
        print('end genome')
    return tups

def create_database(raw, infile, outfile):
    kmers = read_infile(infile)
    size = len(kmers)
    mod = size / 20
    
    printd('Retrieving kmer locations...')
    with open(outfile, 'w') as f:
        for count, kmer in enumerate(kmers):
            # if count % mod == 0:
            #     printd(f'\t{count / size * 100}%')
            printd(f'\tprocessing kmer: {count}')
            val = search_for_kmer(raw, kmer)
            s = ' '.join(val)
            outfile.write(kmer)
            outfile.write(s)
            outfile.write('\n')

def run_gwas():
    raw = load_raw()
    seqs = create_database(raw, DSK_OUTPUT, DATABASE_OUTPUT)
    #dump_seqs(seqs)
    printd('Done.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runs GWAS on microbial genomes')
    parser.add_argument('-d', action='store_true',
        help='turn on debug mode')
    args = parser.parse_args()
    DEBUG = args.d

    run_gwas()




# def create_database(raw, conn):

#     K = 30
#     count = 1
#     printd('Creating kmer database...')
#     cur = conn.cursor()
#     for raw_id, raw_dict in raw.items():
        
#         data = []
#         seq = raw_dict['seq']
#         for c_id, contig in enumerate(seq):
#             l = len(contig)
#             if l >= K: # ensure this contig is long enough to sample
#                 for i in range(l - K + 1):
#                     kmer = contig[i:i + K] # sample the kmer itself
#                     data.append(f'{kmer}\t{str(raw_id)}\t{str(c_id)}\t{str(i)}')

#         data = io.StringIO('\n'.join(data))
#         cur.copy_from(data, 'kmer')
#         conn.commit()


#         printd('Processed genome', count)
#         count += 1
#         #if count > 1: break
#     cur.close()
#     return 1


# def find_seqs(raw, conn):

#     K = 30
#     num_samples = len(raw)
#     count = 1
#     printd('Creating kmer database...')
#     cur = conn.cursor()
#     for raw_id, raw_dict in raw.items():
        
#         data = []
#         seq = raw_dict['seq']
#         for c_id, contig in enumerate(seq):
#             l = len(contig)
#             if l >= K: # ensure this contig is long enough to sample
#                 for i in range(l - K + 1):
#                     kmer = contig[i:i + K] # sample the kmer itself
#                     data.append((kmer, str(raw_id), str(c_id), str(i)))

#                     # store the id, contig number, and string index
#                     # if kmer in kmers:
#                     #     kmers[kmer][raw_id] = (c_id, i)
#                     # else:
#                     #     kmers[kmer] = {raw_id: (c_id, i)}
#         print('created list')
#         data = io.StringIO('\n'.join(['\t'.join(tup) for tup in data]))
#         print('created string')
#         #print(data.readlines())
#         cur.copy_from(data, 'kmer')
#         #upsert_kmer(cur, data)
#         conn.commit()


#         print('Processed genome', count)
#         count += 1
#         if count > 1: break
#     cur.close()

#     print(len(kmers))
#     return 1

#     printd('Identifying resistant kmers...')
#     for kmer, kmer_dict in kmers.items():
        
#         # if this kmer is present in >99% or <1% of samples, ignore it
#         p_kmer = len(kmer_dict) / num_samples
#         if p_kmer > .99 or p_kmer < .01:
#             break

#         # if <95% of the samples that have this kmer are resistant, ignore it
#         num_resistant = 0
#         for raw_id, loc  in kmer_dict.items():
#             if raw[raw_id]['resistance'] == 1:
#                 num_resistant += 1
#         p_resistant = num_resistant / len(kmer_dict)
#         if p_resistant > 0.95:
#             resistant_kmers[kmer] = kmers[kmer]
    
    # print('Consolidating adjacent kmers...')
    # resistant_kmers = self.remove_duplicates(resistant_kmers)

    # printd(f'{len(resistant_kmers)} resistant sequences detected')
    # return resistant_kmers