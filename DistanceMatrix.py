import random
import math
import numpy as np
import pickle

#Helper function for getting the difference 
#Input: 2 stings to compare
#Output: Percent difference between 
def diff(s1, s2):
    if len(s1) < len(s2):
        return diff(s2, s1)
    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 
            deletions = current_row[j] + 1 
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    
    return previous_row[-1]

def distanceMatrix(kmerDict):
	KmerOfOneSize = {}
	#list of the 1 percents
	#print len(kmerDict) #checks len of kmer TEST
	#onePercent = int(math.ceil(0.01*len(kmerDict))) #get 1 percent of len(kmerDict)
	for key, value in kmerDict.items():
		if len(key) == 13:
			temp = {key:value}
			KmerOfOneSize.update(temp)
	#print(KmerOfOneSize)
	onePercent = int(math.ceil(0.18*len(KmerOfOneSize))) #get 1 percent of len(kmerDict)
	#print onePercent #TEST
	p = random.sample(KmerOfOneSize.items(), onePercent)#add random 1 percent of kmerDict to a new dict for the matrix
	print( "Getting random sample...")
	#print(p)
	print(diff("AGGTDS", "AGHTDS")/6)
	x,y = zip(*p)
	print("this is x")
	print(x)
	my_array = np.zeros(shape=((len(x),len(x))))
	#time = time()
	print ("Building distance matrix...")
	for i in range(0, len(x)):
		for j in range(i+1, len(x)):
			print(y[i])
			print(" ")
			print(y[j])
			print(diff(y[i],y[j]))
			my_array[i][j] = str(round(diff(y[j], y[i])/len(y[j], 3))
	for i in range(0, len(x)):
		for j in range(0, len(x)):
			if i ==j:
				my_array[i][j] = 0
			elif i  > j:
				my_array[i][j] = my_array[j][i]
	#t1 = time()
	print (my_array)

with open('raw.pickle', 'rb') as f:
	raw = pickle.load(f)

distanceMatrix(raw)


