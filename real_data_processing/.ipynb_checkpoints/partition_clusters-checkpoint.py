import editdistance
import matplotlib.pyplot as plt
import sys
import timeit
import distance 
import numpy as np

EDIT_DIST_THRESHOLD = 20
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
    refs = []
    with open(file_name) as fh:
        while True:
            idx = fh.readline().rstrip() # index and primers
            if len(idx) == 0:
                break
            idx = idx.split('_')
            seq = fh.readline().rstrip()[PRIMER_LEN:-PRIMER_LEN] # sequence
            #p1 = idx[2]
            #p2 = idx[3]
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
    '''approximate matching with edit distance'''
    pos = -1
    min_dist = len(p)
    for i in range(len(t) - len(p) + 1):
        if editdistance.eval(t[i:i + len(p)], p) <= min_dist:
            min_dist = editdistance.eval(t[i:i + len(p)], p)
            if min_dist <= threshold:
                pos = i  
    return pos 


def extract_noisy_copy(seq, ref, pos):
    l = len(ref)
    p = pos
    min_dist = editdistance.eval(seq[p:p + l], ref)
    
    # removing junks at the beginning
    while True:
        p += 1
        # if seq[p:p + l] == ref:
            
        
    # removing junks at the end
    
    # extending to cover potential omitted neucleotides


def find_reference(seq, refs):
    '''find the reference strand for sequence seq'''
    print('--------')
    for [r, _, _] in refs:
        r_rev = reverse_complement(r)
        pos = approximate_match(r, seq)
        if pos != -1:
            print(pos)
            print("origin", r)
            print("strand", seq[pos:pos+len(r)])
            return
        else:
            pos = approximate_match(r_rev, seq)
            if pos != -1:
                print(pos)
                print("origin", r_rev)
                print("strand", seq[pos:pos+len(r)])
                return



        
    
