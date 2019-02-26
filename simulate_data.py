import random
import pickle

PICKLE_RAW = 'raw.pickle'
PICKLE_SIM = 'sim.pickle'
PERCENT_RESISTANT = .2

def load_raw():
	print('Loading raw data...')
	with open(PICKLE_RAW, 'rb') as f:
		raw = pickle.load(f)
	return raw
def identify_variant_location(raw):
	# create a database of kmers of length s.t. all samples have exactly 1 of this kmer
	# identify a location relative to the start of that kmer
	# return the relative location in all genomes (returns a list or tuple since location is different
	# in every genome)

	K = 750
	kmers = {}
	count = 0
	print('Creating kmer database...')
