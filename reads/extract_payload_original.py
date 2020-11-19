import editdistance
import matplotlib.pyplot as plt
import sys
import timeit

# config
EDIT_DIST_THRESHOLD = 3
READS_FILE = sys.argv[1]
PRIMERS_FILE = sys.argv[2]
STRAND_LEN = 110

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
    return pos


def extract_payloads(primer1, primer2, reads_file_name):
    reads = read_sequences(reads_file_name)
    f1 = open("5pto3p.txt", "a")
    f2 = open("3pto5p.txt", "a")
    direct_count = 0
    rev_count = 0
    hist = [0] * 170
    for r in reads:
        p1 = approximate_match(primer1, r)
        p2 = approximate_match(primer2, r)
        pr1 = approximate_match(reverse_complement(primer1), r)
        pr2 = approximate_match(reverse_complement(primer2), r)
        # print(p1, p2, pr1, pr2)
        if p1 != -1 and p2 != -1:
            payload = r[p1+len(primer1):p2] 
            if len(payload) >= 0.9 * STRAND_LEN and len(payload) <= 1.1 * STRAND_LEN:
                f1.write(payload + '\n')
                direct_count += 1
                hist[len(payload)] += 1
        elif pr1 != -1 and pr2 != -1:
            payload = r[pr2+len(primer2):pr1]  
            if len(payload) >= 0.9 * STRAND_LEN and len(payload) <= 1.1 * STRAND_LEN:
                f2.write(payload + '\n')
                rev_count += 1
                hist[len(payload)] += 1
        else:
            # print('Failed\n')
            continue
    f1.close()
    f2.close()
    
    # just for test
    print('5->3 strands count', direct_count)
    print('3->5 strands count', rev_count)
    print('percentage of excluded strands:', (1 - (direct_count + rev_count) / len(reads))* 100)
    # print(hist)
    plt.plot(range(int(0.9 * STRAND_LEN), int(1.1 * STRAND_LEN)), hist[int(0.9 * STRAND_LEN):int(1.1 * STRAND_LEN)])
    plt.title('distribution of strands lengths')
    plt.show()

if __name__ == "__main__": 
    primer1, primer2 = read_primers(PRIMERS_FILE)  
    start = timeit.default_timer()
    extract_payloads(primer1, primer2, READS_FILE)
    stop = timeit.default_timer()
    print('Time (sec): ', stop - start) 