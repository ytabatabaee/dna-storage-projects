import editdistance
import matplotlib.pyplot as plt
import sys
import timeit

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


def extract_payloads(primer1, primer2, reads_file_name):
    reads = read_sequences(reads_file_name)
    direct_list = []
    rev_list = []
    primer1_rev = reverse_complement(primer1)
    primer2_rev = reverse_complement(primer2)
    direct_count = 0
    rev_count = 0
    hist = [0] * 170
    for r in reads:
        p1, p2 = approximate_match(primer1, primer2, r)
        if p2 != -1 and p2 != -1:
            payload = r[p1+len(primer1):p2] 
            if len(payload) >= 0.9 * STRAND_LEN and len(payload) <= 1.1 * STRAND_LEN:
                direct_list.append(payload)
                direct_count += 1
                hist[len(payload)] += 1
        else:
            pr2, pr1 = approximate_match(primer2_rev, primer1_rev, r)
            if pr1 != -1 and pr2 != -1:
                payload = r[pr2+len(primer2):pr1]  
                if len(payload) >= 0.9 * STRAND_LEN and len(payload) <= 1.1 * STRAND_LEN:
                    rev_list.append(payload)
                    rev_count += 1
                    hist[len(payload)] += 1
            else:
                continue           
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
    primer1, primer2 = read_primers(PRIMERS_FILE)  
    start = timeit.default_timer()
    extract_payloads(primer1, primer2, READS_FILE)
    stop = timeit.default_timer()
    print('Time (sec): ', stop - start)