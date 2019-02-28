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
import psycopg2
import psycopg2.extras
import io

PICKLE_RAW = 'pickle/sim.pickle' #'pickle/raw.pickle'
PICKLE_KMERS = 'pickle/kmers.pickle'
PICKLE_SEQS = 'pickle/seqs.pickle'
DEBUG = False

def printd(*args):
    if DEBUG == True:
        print(*args)

# load the gene pool file
def load_raw():
    printd('Loading data...')
    with open(PICKLE_RAW, 'rb') as f:
        raw = pickle.load(f)
    return raw

def create_database(raw, conn):

    K = 30
    count = 1
    printd('Creating kmer database...')
    cur = conn.cursor()
    for raw_id, raw_dict in raw.items():
        
        data = []
        seq = raw_dict['seq']
        for c_id, contig in enumerate(seq):
            l = len(contig)
            if l >= K: # ensure this contig is long enough to sample
                for i in range(l - K + 1):
                    kmer = contig[i:i + K] # sample the kmer itself
                    data.append(f'{kmer}\t{str(raw_id)}\t{str(c_id)}\t{str(i)}')

        data = io.StringIO('\n'.join(data))
        cur.copy_from(data, 'kmer')
        conn.commit()


        printd('Processed genome', count)
        count += 1
        if count > 1: break
    cur.close()
    return 1


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


def dump_seqs(seqs):
    printd('Dumping seqs data...')
    with open(PICKLE_SEQS, 'wb') as f:
        pickle.dump(seqs, f)

def run_gwas():
    raw = load_raw()

    # need to be created AFTER THE FORK
    conn = psycopg2.connect('dbname=kmer user=localuser')
    
    seqs = create_database(raw, conn)
    conn.close()
    #dump_seqs(seqs)
    printd('Done.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runs GWAS on microbial genomes')
    parser.add_argument('-d', action='store_true',
        help='turn on debug mode')
    args = parser.parse_args()
    DEBUG = args.d

    run_gwas()



