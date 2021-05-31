'''
This program constructs a bloom filter for all 30-mers in the human genome.
Together with the next program (query.py) it can be used to find whether a 
DNA sequence is human. 


# ------- Bloom filter parameters computation ------- #

- Desired false positive rate: p <= 1% = 0.01
- Size of the input: n = 3*10^9 = 3000000000 (number of 30-mers present in human.fsa)
- Size of the filter: m = (-n*ln(p)) / (ln(2)^2) = (-3*10^9*ln(0.01)) / (ln(2)^2) 
                        = 28755175132.1 ~ 34359738368 bits 
                        --> log2(m) = 35 bits needed from the hash function
- Number of hash functions = ceiling(-ln(0.01)/ln(2)) = 7
- New false positive rate: p = (1-e^(-7*3000000000/34359738368))^7 
                             = 0.00418160394 (smaller than what we want)
- New number of hash functions: 0.01 >= (1-e^(-k*3000000000/34359738368))^k = (1-e^(-k*0.08731149137))^k 
                                     --> k = 4: (1-e^(-4*0.08731149137))^4 = 0.00755082013 < 0.01
                                     --> k = 3: (1-e^(-3*0.08731149137))^3 = 0.01223673095 > 0.01

Therefore, our parameters are: 
n = 3000000000 elements
m = 34359738368 bits
k = 4 hash functions
p = 0.76% false positive rate
bitfield = 35
'''

import sys, time, re, hashlib


## ---- Bloom filter and parameters ---- ##
n = 3000000000
m = 34359738368
k = 4 
p = 0.0076 
bitfieldsize = 35
bloomfilter = bytearray(round(m/8)) # m bits long, packed into m/8 bytes

## ---- Hash function choice ---- ##
# We need a hash function that provides a hash that is at least
# k * bitfield size long (35*4). sha256 is a good hash function 
# since it returns a 256 bit hash

## ---- Approach ---- ##
# We run each kmer through the hash function. For each k, we take a slice 
# from the hash of the size of the bitfieldsize. The 3 most right values of 
# the slice are used to generate the bitposition, which will be a number 
# between 0 and 7. The other values are used to generate the byteposition, 
# which will be a number between 0 and m/8.


## ---- Bloomfilter functions ---- ##
def hashit(kmer):
    ''' 
    This function takes a k-mer (as a bytes object) and gives back an integer 
    representing a 32 byte (256 bit) integer number - the hash-number.
    '''
    myhashobj = hashlib.sha256()
    myhashobj.update(kmer)
    return int.from_bytes(myhashobj.digest(), byteorder=sys.byteorder)

def nextposition(hashnumber):
    '''
    This function finds the next position of the bit to set, given a hash-number.
    '''
    position = hashnumber & (2**bitfieldsize - 1)
    byteposition = position >> 3
    bitposition = position & 7
    hashnumber >>= bitfieldsize # discarding the used hash slice
    return(hashnumber, byteposition, bitposition)

def setbit(bloomfilter, byteposition, bitposition):
    '''
    This function sets a bit in a bytearray given positions.
    '''
    bloomfilter[byteposition] |= 1 << bitposition


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

## ---- Time taking functions ---- ##
def timepoint(comment, timelist):
    now = time.time()
    timelist.append((comment, now))

def showtimepoints(timelist):
    m = max([ len(timelist[i][0]) for i in range(1, len(timelist)) ])
    formatstring = "{:" + str(m+1) +"} {:.02f} seconds"
    for i in range(1, len(timelist)):
        print(formatstring.format(timelist[i][0]+':', timelist[i][1]-timelist[i-1][1]))
    print(formatstring.format('Total:', timelist[-1][1]-timelist[0][1]))

## ---- Build Bloom Filter ---- ##
def addtobloomfilter(filename, fastaindex, kmer):
    # Translation table so everything not ATGG becomes N.
    toN = bytes.maketrans(b'MRYKVHDBacgtmrykvhdbxnsw',b'NNNNNNNNACGTNNNNNNNNNNNN')
    # Read file
    try:
        infile = open(filename, 'rb')
    except IOError as err:
        print("Failed opening file:", str(err))
        sys.exit(1)
    # Find kmers and add them to the filter
    for idx in fastaindex:
        infile.seek(idx[2])
        seq = infile.read(idx[3]-idx[2]+1).translate(toN, b'\r\n\t ')
        # Index sequence
        seqindex = indexsequence(seq)
        # Add kmers to filter
        for start, stop in seqindex:
            for i in range(start, stop-kmer+1):
                mer = seq[i:i+kmer]
                hashnumber = hashit(mer)
                for i in range(0, k): # k hash slices
                    hashnumber, byteposition, bitposition = nextposition(hashnumber)
                    setbit(bloomfilter, byteposition, bitposition)
    infile.close()


## ---- Main program ---- ##
if __name__ == '__main__':
    # Start program
    timelist = list()
    timepoint('Start', timelist)
    
    # What file to find Kmers in
    filename = '/home/projects/pr_course/human.fsa'
    
    # What Kmers to find
    kmer = 30
    
    # Index fasta file
    fastaindex = indexfasta(filename)
    timepoint('Indexing file', timelist)
    
    # Bloom filter construction
    addtobloomfilter(filename, fastaindex, kmer)
    timepoint('Bloomfilter construction', timelist)

    # Save in files
    output = open("/home/projects/pr_course/people/celbur/week12/bloomfilter.txt", 'wb')
    output.write(bloomfilter)
    output.close()

    # Display timing
    showtimepoints(timelist)

