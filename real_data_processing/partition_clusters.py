import editdistance
import matplotlib.pyplot as plt
import sys
import timeit
import distance 
import numpy as np

EDIT_DIST_THRESHOLD = 25
PRIMER_LEN = 25

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
    return sequences


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
    return refs


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


def approximate_match(p, t, threshold=EDIT_DIST_THRESHOLD):
    '''approximate matching with edit distance (finding p in t)'''
    pos = -1
    min_dist = len(p)
    for i in range(len(t) - len(p) + 1):
        if editdistance.eval(t[i:i + len(p)], p) <= min_dist:
            min_dist = editdistance.eval(t[i:i + len(p)], p)
            if min_dist <= threshold:
                pos = i  
    return pos 


def extract_noisy_copy(seq, ref, pos):
    # TODO: this is very naive, should be optimized
    l = len(ref)
    p = pos
    min_dist = editdistance.eval(seq[p:p + l], ref)
    print(seq[p:p + l])
    print(min_dist)
    
    # removing junks at the beginning
    while True:
        ed = editdistance.eval(seq[p+1:p+l], ref) 
        if ed < min_dist:
            print(ed)
            p += 1
            l -= 1
            min_dist = ed
        else:
            break
    print(seq[p:p + l])

    # extending to cover potential omitted neucleotides
    while True:
        ed = editdistance.eval(seq[p:p+l+1], ref) 
        if ed < min_dist:
            print(ed)
            l += 1
            min_dist = ed
        else:
            break
    print(seq[p:p + l])
    
    # removing junks at the end
    while True:
        ed = editdistance.eval(seq[p:p+l-1], ref)
        if ed < min_dist:
            print(ed)
            l -= 1
            min_dist = ed
        else:
            break     
    print(seq[p:p + l])
    
    
    
def find_reference(seq, refs):
    '''find the reference strand for sequence seq'''
    print('--------')
    for r in refs:
        r_rev = reverse_complement(r)
        pos = approximate_match(r, seq)
        if pos != -1:
            print(pos)
            print("origin", r)
            print("strand", seq[pos:pos+len(r)])
            return pos, r
        else:
            pos = approximate_match(r_rev, seq)
            if pos != -1:
                print(pos)
                print("origin", r_rev)
                print("strand", seq[pos:pos+len(r)])
                return pos, r
    return None



        
    
