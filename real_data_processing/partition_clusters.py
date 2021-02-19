import editdistance
import matplotlib.pyplot as plt
import sys
import timeit
import distance 
import numpy as np

EDIT_DIST_THRESHOLD = 25
PRIMER_LEN = 25
REFS_FILE = sys.argv[1]
SEQS_FILE = sys.argv[2]

def read_sequences(file_name):
    '''read sequences from fastq file'''
    sequences = []
    with open(file_name) as fh:
        while True:
            fh.readline() 
            seq = fh.readline().rstrip() # read sequence
            fh.readline()
            fh.readline()
            if len(seq) == 0:
                break
            if 'N' not in seq:
                sequences.append(seq)
    return list(set(sequences))


def read_refs(file_name):
    '''reading references payloads'''
    refs = []
    with open(file_name) as fh:
        while True:
            idx = fh.readline() # index and primers
            seq = fh.readline().rstrip()[PRIMER_LEN:-PRIMER_LEN] # sequence
            if len(seq) == 0:
                break
            refs.append(seq)
    return list(set(refs))


def write_refs_clusters(refs, clusters):
    '''write references and clusters to file'''
    with open('data/reference.txt', 'w') as fr:
        for i in range(len(refs)):
            r = refs[i]
            fr.write(r + "\n")
            with open('data/clusters/' + str(i) + '.txt', 'w') as fc:
                for c in clusters[r]:
                    fc.write(c + "\n")
                        

def reverse_complement(seq):
    base_pairs = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(base_pairs.get(base, base) for base in reversed(seq))


def prune_sequences(seqs, threshold = 0.1):
    '''prune outlier sequences based on length '''
    hist = get_length_histogram(seqs)
    peak = np.argmax(hist)
    pruned_seqs = []
    for s in seqs:
        if len(s) >= peak * (1-threshold) and len(s) <= peak * (1+threshold):
            pruned_seqs.append(s)
    return pruned_seqs
        
        
def get_length_histogram(seqs, max_len = 1200):
    hist = [0]*max_len
    for s in seqs:
        hist[len(s)] += 1 
    return hist

def get_cluster_distribution(clusters, max_len = 150):
    hist = [0]*max_len
    for c in clusters:
        hist[len(clusters[c])] += 1 
    return hist 


def approximate_match(p, t, threshold=EDIT_DIST_THRESHOLD):
    '''approximate matching with edit distance'''
    pos = -1
    min_dist = len(p)
    for i in range(len(t) - len(p) + 1):
        ed = editdistance.eval(t[i:i + len(p)], p)
        if ed <= min_dist:
            min_dist = ed
            if min_dist <= threshold:
                pos = i  
    return pos 


def extract_noisy_copy(seq, ref, pos):
    '''finding the most likely noisy copy of ref within sequence seq'''
    # TODO: this is very naive, should be optimized
    l = len(ref)
    p = pos
    min_dist = editdistance.eval(seq[p:p + l], ref)
    #print(seq[p:p + l])
    #print(min_dist)
    
    for i in range(-10, 10):
        for j in range(-10, 10):
            ed = editdistance.eval(seq[pos+i:pos+l+j], ref)
            if ed < min_dist:
                min_dist = ed
                l = len(seq[pos+i:pos+l+j])
                p = pos+i
    #print(seq[p:p + l])
    #print(min_dist)
    return seq[p:p + l]    
    
    
def find_reference(seq, refs):
    '''find the reference strand for sequence seq'''
    # TODO: how should we handle reverse complement refs? should they form one cluster?
    #print('--------')
    for r in refs:
        s_rev = reverse_complement(seq)
        pos = approximate_match(r, seq)
        if pos != -1:
            #print(pos)
            #print("origin", r)
            #print("strand", seq[pos:pos+len(r)])
            return seq, r, pos
        else:
            pos = approximate_match(r, s_rev)
            if pos != -1:
                #print(pos)
                #print("origin", r)
                #print("strand", s_rev[pos:pos+len(r)])
                return s_rev, r, pos
    return None


def extract_clusters(seqs, refs):
    clusters = dict()
    for r in refs:
        clusters[r] = [] 
     
    for seq in seqs[:300]:
        try:
            s, r, pos = find_reference(seq, refs)
        except:
            continue
        noisy_copy = extract_noisy_copy(s, r, pos) 
        clusters[r].append(noisy_copy)
        #print(r)
        #print(noisy_copy)
        #print('---')
    
    return clusters

if __name__ == "__main__": 
    seqs = read_sequences(SEQS_FILE)
    refs = read_refs(REFS_FILE) 
    print(len(refs))
    print(len(seqs))
    print('finished reading data')
    start = timeit.default_timer()
    clusters = extract_clusters(seqs, refs)
    hist = get_cluster_distribution(clusters)
    print(hist)
    print("processing time: ", (stop - start))
    write_refs_clusters(refs, clusters)
    stop = timeit.default_timer()
