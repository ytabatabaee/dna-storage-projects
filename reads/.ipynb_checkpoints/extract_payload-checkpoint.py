import editdistance
import matplotlib.pyplot as plt
import sys
import timeit
import pickle
from multiprocessing import Manager, Pool

#primer1 = ""
#primer2 = ""
primer1_rev = ""
primer2_rev = ""
direct_count = 0
rev_count = 0

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


def approximate_match(p, t, threshold=EDIT_DIST_THRESHOLD):
    '''approximate matching with edit distance'''
    pos = -1
    min_dist = len(p)
    for i in range(len(t) - len(p) + 1):
        if editdistance.eval(t[i:i + len(p)], p) <= min_dist:
            min_dist = editdistance.eval(t[i:i + len(p)], p)
            if min_dist <= threshold:
                pos = i  
            if min_dist == 0:
                break
    return pos


def process_read(r):
    #global direct_list
    #global rev_list
    #global direct_count
    #global rev_count
    p1 = approximate_match(primer1, r)
    if p1 != -1:
        p2 = approximate_match(primer2, r)
        if p2 != -1:
            payload = r[p1+len(primer1):p2] 
            if len(payload) >= 0.9 * STRAND_LEN and len(payload) <= 1.1 * STRAND_LEN:
                return payload
                #direct_list.append(payload)
                # f1.write(payload + '\n')
                #direct_count += 1
                hist[len(payload)] += 1
    else:
        pr1 = approximate_match(primer1_rev, r)
        if pr1 != -1:
            pr2 = approximate_match(primer2_rev, r)
            if pr2 != -1:
                payload = r[pr2+len(primer2):pr1]  
                if len(payload) >= 0.9 * STRAND_LEN and len(payload) <= 1.1 * STRAND_LEN:
                    return payload
                    #rev_list.append(payload)
                    # f2.write(payload + '\n')
                    #rev_count += 1
                    hist[len(payload)] += 1

                        
def extract_payloads(reads_file_name):
    global primer1_rev, primer2_rev
    global direct_list
    global rev_list
    reads = read_sequences(reads_file_name)
    primer1_rev = reverse_complement(primer1)
    primer2_rev = reverse_complement(primer2)
    hist = [0] * 170
    pool = Pool(processes=2)
    direct_list = pool.map_async(process_read, reads)
    pool.close()
    pool.join()
    rev_list = direct_list
    print(len(direct_list))
    with open(DIR_FILE, "w") as f1:
        f1.write("\n".join(direct_list))  
    with open(REV_FILE, "w") as f2:
        f2.write("\n".join(rev_list))      

    # just for test
    print('5->3 strands count', direct_count)
    print('3->5 strands count', rev_count)
    print('percentage of excluded strands:', (1 - (direct_count + rev_count) / len(reads))* 100)
    # print(hist)
    #plt.plot(range(int(0.9 * STRAND_LEN), int(1.1 * STRAND_LEN)), hist[int(0.9 * STRAND_LEN):int(1.1 * STRAND_LEN)])
    #plt.title('distribution of strands lengths')
    #plt.show()

if __name__ == "__main__": 
    # global primer1, primer2
    manager1 = Manager()
    direct_list = manager1.list()
    manager2 = Manager()
    rev_list = manager2.list()
    primer1, primer2 = read_primers(PRIMERS_FILE)  
    start = timeit.default_timer()
    extract_payloads(READS_FILE)
    stop = timeit.default_timer()
    print('Time (sec): ', stop - start)  
