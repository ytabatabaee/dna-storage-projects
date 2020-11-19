import editdistance
import matplotlib.pyplot as plt
import sys
import timeit
from multiprocessing import Manager, Pool, cpu_count

# config
EDIT_DIST_THRESHOLD = 3
READS_FILE = sys.argv[1]
PRIMERS_FILE = sys.argv[2]
STRAND_LEN = 110
DIR_FILE = '5pto3p.txt'
REV_FILE = '3pto5p.txt'

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


def read_primers(file_name):
    '''read primers from reference'''
    with open(file_name) as fh:
        seq = fh.readline().rstrip()
        primer1 = seq[3:23]
        primer2 = seq[-23:-3]
    return primer1, primer2


def reverse_complement(s):
    base_pairs = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    t = ''
    for base in s:
        t = base_pairs[base] + t
    return t


def approximate_match(pl, pr, t, threshold=EDIT_DIST_THRESHOLD):
    '''approximate matching with edit distance'''
    p1, p2 = -1, -1
    min_dist = len(pl)
    for i in range(len(t) - len(pl) + 1):
        if editdistance.eval(t[i:i + len(pl)], pl) <= min_dist:
            min_dist = editdistance.eval(t[i:i + len(pl)], pl)
            if min_dist <= threshold:
                p1 = i  
            if min_dist == 0:
                break
    if p1 == -1:
        return p1, p2
    
    min_dist = len(pr)
    for i in reversed(range(len(t) - len(pr) + 1)):
        if editdistance.eval(t[i:i + len(pr)], pr) <= min_dist:
            min_dist = editdistance.eval(t[i:i + len(pr)], pr)
            if min_dist <= threshold:
                p2 = i  
            if min_dist == 0:
                break
    return p1, p2


def process_read(r):
    global direct_list
    global rev_list
    p1, p2 = approximate_match(primer1, primer2, r)
    if p2 != -1 and p2 != -1:
        payload = r[p1+len(primer1):p2] 
        if len(payload) >= 0.9 * STRAND_LEN and len(payload) <= 1.1 * STRAND_LEN:
            direct_list.append(payload)
    else:
        pr2, pr1 = approximate_match(primer2_rev, primer1_rev, r)
        if pr1 != -1 and pr2 != -1:
            payload = r[pr2+len(primer2):pr1]  
            if len(payload) >= 0.9 * STRAND_LEN and len(payload) <= 1.1 * STRAND_LEN:
                rev_list.append(payload) 


def extract_payloads(primer1, primer2, reads_file_name):
    global primer1_rev, primer2_rev
    reads = read_sequences(reads_file_name)
    primer1_rev = reverse_complement(primer1)
    primer2_rev = reverse_complement(primer2)

    pool = Pool(processes=cpu_count())
    pool.map_async(process_read, reads)
    pool.close() 
    pool.join()

    with open(DIR_FILE, "w") as f1:
        f1.write("\n".join(direct_list))  
    with open(REV_FILE, "w") as f2:
        f2.write("\n".join(rev_list))  
    
    # just for test
    print('5->3 strands count', len(direct_list))
    print('3->5 strands count', len(rev_list))
    print('percentage of excluded strands:', (1 - (len(direct_list) + len(rev_list)) / len(reads))* 100)
    return len(reads)

if __name__ == "__main__": 
    primer1_rev = ""
    primer2_rev = ""
    manager1 = Manager()
    direct_list = manager1.list()
    manager2 = Manager()
    rev_list = manager2.list()
    primer1, primer2 = read_primers(PRIMERS_FILE)  
    start = timeit.default_timer()
    len_reads = extract_payloads(primer1, primer2, READS_FILE)
    stop = timeit.default_timer()
    print('strands/sec: ', len_reads/(stop - start))
