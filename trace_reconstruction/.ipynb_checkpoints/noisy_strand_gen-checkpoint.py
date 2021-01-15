import random


def random_string(alphabet, length):
    return ''.join(random.choices(alphabet, k=length))


def weighted_choice(choices):
    r = random.uniform(0, sum(w for c, w in choices))
    upto = 0
    for c, w in choices:
        if upto + w >= r:
            return c
        upto += w
    assert False, "Shouldn't get here"


def add_noise(strand, choices, bases):
    out = []
    for c in strand:
        error = weighted_choice(choices)
        if 'sub' == error:
            out.append(random.choice(bases))
        elif 'in' == error:
            out.append(random.choice(bases))
            out.append(c)
        elif 'skip' == error:
            out.append(c)
    return ''.join(out)


def generate_strand_error(original, strand_num, strand_len, choices, bases):
    cluster = [original] * strand_num
    rev_cluster = []
    for i in range(0, len(cluster)):
        cluster[i] = add_noise(cluster[i], choices, bases) + '-' * (strand_len - len(cluster[i]))
        rev_cluster.append(cluster[i][::-1])

    return cluster, rev_cluster

