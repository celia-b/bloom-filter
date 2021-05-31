'''
This program runs a DNA sequence through the bloom filter
and returns the proportion of found kmers that are also present
in the human genome, each with a probability of 99%.
'''

#!/usr/bin/env python3
import sys, time, re, hashlib

## ---- Bloomfilter parameters ---- ##
n = 3000000000
m = 34359738368
k = 4 
p = 0.0076 
bitfieldsize = 35

## ---- Bloomfilter functions ---- ##
def hashit(kmer):
    myhashobj = hashlib.sha256()
    myhashobj.update(kmer)
    return int.from_bytes(myhashobj.digest(), byteorder=sys.byteorder)

def nextposition(hashnumber):
    position = hashnumber & (2**bitfieldsize - 1)
    byteposition = position >> 3
    bitposition = position & 7
    hashnumber >>= bitfieldsize # discarding the used hash slice
    return(hashnumber, byteposition, bitposition)

def isbitset(bloomfilter, byteposition, bitposition):
        return (bloomfilter[byteposition] & (1 << bitposition)) != 0
        #except: 
        #    return False


## ---- Read in the bloomfilter ---- ##
infile = open('/home/projects/pr_course/people/celbur/week12/bloomfilter.txt', 'rb')
bloomfilter = infile.read()
infile.close()


## ---- Index the fasta file ---- ##
def indexfasta(filename):
    try:
        infile = open(filename, 'rb')
    except IOError as err:
        print("Cant open file:", str(err));
        sys.exit(1)
    chunksize = 1024*1024
    filepos = 0
    headstart = list()
    headend = list()
    while True:
        content = infile.read(chunksize)
        if len(content) == 0:
            break
        # find headers
        chunkpos = 0
        while chunkpos != -1:
            chunkpos = content.find(b'>', chunkpos)
            if chunkpos != -1:
                headstart.append(chunkpos + filepos)
                chunkpos += 1
        # find corresponding headend
        for i in range(len(headend), len(headstart)):
            chunkpos = max(0, headstart[i] - filepos)
            chunkpos = content.find(b'\n', chunkpos)
            if chunkpos != -1:
                headend.append(chunkpos + filepos)
        filepos += len(content)
    infile.close()
    # Eliminating wrong headers due to extra > in header line
    for i in range(len(headstart)-1, 0, -1):
        if headend[i] == headend[i-1]:
            del headstart[i]
            del headend[i]
    headstart.append(filepos)
    fastaindex = list()
    for i in range(len(headend)):
        fastaindex.append((headstart[i], headend[i], headend[i]+1, headstart[i+1] - 1))
    return fastaindex

## ---- Index the sequence ---- ##
def indexsequence(seq):
    pointer = 0
    seqindex = list()
    while len(seq) > pointer:
        ## Find start of seq
        potenstart = [seq.find(b'A', pointer), seq.find(b'T', pointer), seq.find(b'G', pointer), seq.find(b'C', pointer)]
        realstart = min(potenstart)                             # It starts where the first A/C/T/G is
        if realstart == -1:                                     # If some of the bases are not in the seq
            potenstart = [ i for i in potenstart if i > -1 ]        # Get the ones that are
            if len(potenstart) == 0:                                # If there are none
                break                                                   # We are done with this seq
            realstart = min(potenstart)                             # Take the first
        ## Find end of seq
        realend = seq.find(b'N', realstart)                     # It ends with the first N
        if realend == -1:                                       # If there are no Ns in the seq
            realend = len(seq)                                      # Take it all
        seqindex.append((realstart, realend))
        pointer = realend
    return seqindex

## ---- Find proportion of human kmers in sample ---- ##
def humansample(filename, fastaindex, kmer):
    # Translation table so everything not ATGG becomes N.
    toN = bytes.maketrans(b'MRYKVHDBacgtmrykvhdbxnsw',b'NNNNNNNNACGTNNNNNNNNNNNN')
    
    # Read file
    try:
        infile = open(filename, 'rb')
    except IOError as err:
        print("Failed opening file:", str(err))
        sys.exit(1)
    
    # Open outfile for results:
    outfile = open('/home/projects/pr_course/people/celbur/week12/isthishuman.txt', 'wb')

    # Inspect each sample
    for idx in fastaindex:
        # Header
        infile.seek(idx[0])
        header = infile.read(idx[2]-idx[0]-1)
        # Sequence
        total_kmers = 0
        human_kmers = 0
        infile.seek(idx[2])
        seq = infile.read(idx[3]-idx[2]+1).translate(toN, b'\r\n\t ')
        # Index sequence
        seqindex = indexsequence(seq)
        # Look for kmers in filter
        for start, stop in seqindex:
            for i in range(start, stop-kmer+1):
                total_kmers += 1
                mer = seq[i:i+kmer]
                hashnumber = hashit(mer)
                there = True # We start by assuming the kmer is in the filter
                for i in range(0, k): # Check the k positions in the filter (all must be 1)
                    hashnumber, byteposition, bitposition = nextposition(hashnumber)
                    if isbitset(bloomfilter, byteposition, bitposition) == False: 
                        there = False # If it's not, we change the assumption
                if there == True:
                    human_kmers += 1
        
        # Save results
        outfile.write(header + b'\n')
        outfile.write(b'Proportion of 30-mers in common with humans: ' + str(human_kmers/total_kmers).encode('utf-8') + b'\n\n')

    infile.close()
    outfile.close()


## ---- Main program ---- ##
if __name__ == '__main__':
    # What file to find Kmers in
    filename = '/home/projects/pr_course/mixeddna.fsa'
    
    # What Kmers to find
    kmer = 30
    
    # Index fasta file
    fastaindex = indexfasta(filename)
    
    # Proportion of human kmers per sample?
    humansample(filename, fastaindex, kmer)

    