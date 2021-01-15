import configparser
import subprocess
import random
import sys
import os
import time
import editdistance
import distance

config_location = sys.argv[1]
test_index = sys.argv[2]
config = configparser.ConfigParser()
config.read(config_location)
parameters = config['parameters']
file_locations = config['file_locations']
flags = config['flags']

strand_num = int(parameters['strand_num'])
strand_length = int(parameters['strand_length'])
error_rate = float(parameters['error_rate'])
window_size = int(parameters['window_size'])
sample_cases = int(parameters['sample_cases'])

bases = ['A', 'C', 'G', 'T']

file_params = "w" + str(window_size) + "n" + str(strand_num) + "l" + str(strand_length) + "t" + str(
    sample_cases) + "p" + str(int(error_rate * 100)) + "_" + str(test_index)


def majority(clu, i):
    ems = []
    for ele in clu:
        if len(ele) > i:
            ems.append(ele[i])
    if len(ems) == 0:
        ems.append(random.choice(bases))
    return max(set(ems), key=ems.count)


def recover_strand_dp(cluster, strand_len):
    ans = ''

    i = 0
    while i < strand_len:
        ch, s_aligned = run_dp(cluster, i)
        first_col = [s_aligned[j][0] for j in range(len(cluster))]
        if first_col == ['-'] * strand_num:
            i += 1
            continue

        for j in range(len(cluster)):
            if ch == '-' and s_aligned[j][0] != '-':  # insertion error
                cluster[j] = cluster[j][:i] + cluster[j][i + 1:]

            elif s_aligned[j][0] != ch and s_aligned[j][0] == '-':  # deletion error
                cluster[j] = cluster[j][:i] + ch + cluster[j][i:]

            elif s_aligned[j][0] != ch:  # substitution error
                cluster[j] = cluster[j][:i] + ch + cluster[j][i + 1:]

        if ch != '-':
            i += 1
            ans += ch

    while len(ans) < strand_len:
        ans += random.choice(bases)
    return ans


def run_dp(cluster, i):
    filename = "test" + file_params + ".txt"
    f = open(filename, "w+")
    for strand in cluster:
        f.write(strand[i:min(i + window_size, strand_length)] + "\n")
    f.close()

    subprocess.run(["./dp" + str(strand_num), file_params], check=True)

    resname = "res" + file_params + ".txt"
    s_aligned = []
    if os.path.exists(resname):
        with open(resname) as g:
            s_aligned = g.readlines()
        s_aligned = [x.strip() for x in s_aligned]
        g.close()

    while s_aligned != [''] * strand_num:
        first_col = [s_aligned[j][0] for j in range(len(cluster))]
        if first_col == ['-'] * strand_num:
            for j in range(len(cluster)):
                s_aligned[j] = s_aligned[j][1:]
        else:
            break

    ch = majority(s_aligned, 0)
    if s_aligned == [''] * strand_num:
        s_aligned = [ch] * strand_num

    return ch, s_aligned


def recover_strand(cluster, strand_len):
    ans = ''
    recovered = ''

    for i in range(0, strand_len - 1):
        ch = majority(cluster, i)

        for j in range(len(cluster)):

            if len(cluster[j]) == i:
                cluster[j] += ch

            if cluster[j][i] != ch:

                ch2 = majority(cluster, i + 1)

                ch3_flag = -1
                if i + 2 < strand_len:
                    ch3_flag = 1
                    ch3 = majority(cluster, i + 2)

                ch4_flag = -1
                if i + 3 < strand_len:
                    ch4_flag = 1
                    ch4 = majority(cluster, i + 3)

                ch5_flag = -1
                if i + 4 < strand_len:
                    ch5_flag = 1
                    ch5 = majority(cluster, i + 4)

                if len(cluster[j]) > i + 2:
                    if cluster[j][i] == ch2 and (ch3_flag == -1 or cluster[j][i + 1] == ch3):  # erasure error
                        cluster[j] = cluster[j][:i] + ch + cluster[j][i:]

                    elif cluster[j][i + 1] == ch and cluster[j][i + 2] == ch2:  # insertion error
                        cluster[j] = cluster[j][:i] + cluster[j][i + 1:]

                    elif cluster[j][i + 1] == ch2 and (ch3_flag == -1 or cluster[j][i + 2] == ch3):  # subs
                        cluster[j] = cluster[j][:i] + ch + cluster[j][i + 1:]

                    elif cluster[j][i + 1] != ch2:

                        if cluster[j][i] == ch3 and (ch4_flag == -1 or cluster[j][i + 1] == ch4):  # erasure
                            cluster[j] = cluster[j][:i] + ch + ch2 + cluster[j][i:]

                        elif len(cluster[j]) > i + 3:
                            if cluster[j][i + 2] == ch3 and (ch4_flag == -1 or cluster[j][i + 3] == ch4):  # subs
                                cluster[j] = cluster[j][:i] + ch + ch2 + cluster[j][i + 1:]

                            elif cluster[j][i + 2] == ch and cluster[j][i + 3] == ch2:  # insertion
                                cluster[j] = cluster[j][:i] + cluster[j][i + 2:]

                            elif cluster[j][i + 1] == ch3 and (ch4_flag == -1 or cluster[j][i + 2] == ch4):
                                cluster[j] = cluster[j][:i] + ch + ch2 + cluster[j][i + 1:]

                            elif cluster[j][i + 1] == ch3 and (ch4_flag == -1 or cluster[j][i + 3] == ch4):
                                cluster[j] = cluster[j][:i] + ch + ch2 + cluster[j][i + 1:]

                            elif cluster[j][i + 2] == ch2 and cluster[j][i + 3] == ch3:
                                cluster[j] = cluster[j][:i] + ch + cluster[j][i + 2:]

                            elif cluster[j][i] == ch3 and (ch4_flag == -1 or cluster[j][i + 3] == ch4):
                                cluster[j] = cluster[j][:i] + ch + ch2 + ch3 + cluster[j][i + 3:]

                            elif cluster[j][i + 2] != ch3:

                                if cluster[j][i] == ch4 and (ch5_flag == -1 or cluster[j][i + 1] == ch5):  # erasure
                                    cluster[j] = cluster[j][:i] + ch + ch2 + ch3 + cluster[j][i:]

                                elif len(cluster[j]) > i + 4:
                                    if cluster[j][i + 3] == ch4 and (
                                            ch5_flag == -1 or cluster[j][i + 4] == ch5):  # subs
                                        cluster[j] = cluster[j][:i] + ch + ch2 + ch3 + cluster[j][i + 1:]

                                    elif cluster[j][i + 3] == ch and cluster[j][i + 4] == ch2:  # insertion
                                        cluster[j] = cluster[j][:i] + cluster[j][i + 3:]

                elif len(cluster[j]) == i + 2:
                    if cluster[j][i] == ch2:  # erasure error
                        cluster[j] = cluster[j][:i] + ch + cluster[j][i:]

                    elif cluster[j][i + 1] == ch2:  # subs
                        cluster[j] = cluster[j][:i] + ch + cluster[j][i + 1:]

                    elif cluster[j][i + 1] == ch:  # insertion error
                        cluster[j] = cluster[j][:i] + cluster[j][i + 1:]

                    else:
                        cluster[j] = cluster[j][:i] + ch

        recovered += ch
        ans = ans[:0] + recovered

    last_ch = majority(cluster, strand_len - 1)
    ans += last_ch

    return ans


def clean_up():
    if os.path.exists("res" + file_params + ".txt"):
        os.remove("res" + file_params + ".txt")
    if os.path.exists("test" + file_params + ".txt"):
        os.remove("test" + file_params + ".txt")


def main():
    input_file = "tests/test" + str(test_index) + ".txt"

    cluster = []
    rev_cluster = []
    prob = [0] * strand_length
    original = ""

    if os.path.exists(input_file):
        with open(input_file) as g:
            cluster = g.readlines()
        cluster = [x.strip() for x in cluster]
        original = cluster[0]
        # print(original)
        cluster = cluster[1::]
        # print(len(cluster))
        g.close()

    # print("noisy strands:")
    start_time = time.time()

    for i in range(0, len(cluster)):
        # print(cluster[i])
        rev_cluster.append(cluster[i][::-1])

    if flags.getboolean('non_dp'):
        mj = recover_strand(cluster, int(strand_length / 2))
        rev_mj = recover_strand(rev_cluster, int(strand_length / 2))

    else:
        mj = recover_strand_dp(cluster, strand_length / 2)
        rev_mj = recover_strand_dp(rev_cluster, (strand_length / 2))

    # print(mj)
    # print(rev_mj)

    rev_rev_mj = rev_mj[::-1]
    mj = mj[0:int(strand_length / 2)] + rev_rev_mj[0:int(strand_length / 2)]

    # print(len(original))
    # print(len(mj))

    for j in range(strand_length):
        if original[j] != mj[j]:
            prob[j] = 1

    f = open("outputs/out" + str(test_index) + ".txt", "w+")
    f.write(str(prob) + '\n')
    f.write(str(editdistance.eval(original, mj)) + '\n')
    f.write(str(distance.hamming(original, mj)))
    f.close()
    
    #print("reconstructed strand: ", mj)
    #print("--- %s seconds ---" % (time.time() - start_time))
    clean_up()


if __name__ == "__main__":
    main()

