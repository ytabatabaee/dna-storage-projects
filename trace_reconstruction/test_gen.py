import dp_gen
import noisy_strand_gen
import configparser
import subprocess
import sys

config_location = sys.argv[1]
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
choices = [('sub', 100 * error_rate),
           ('in', 100 * error_rate),
           ('del', 100 * error_rate),
           ('skip', 300 - 300 * error_rate)]

file_params = "w" + str(window_size) + "n" + str(strand_num) + "l" + str(strand_length) + "t" + str(
    sample_cases) + "p" + str(int(error_rate * 100))


def compile_dp():
    file_name = "dp" + str(strand_num) + ".cpp"
    f = open(file_name, "w+")
    f.write(dp_gen.generate_msa(strand_num, strand_length))
    f.close()
    cmd = 'g++ dp' + str(strand_num) + '.cpp -o dp' + str(strand_num)
    subprocess.check_call(cmd, shell=True)


def main():
    if not flags.getboolean('non_dp'):
        compile_dp()
    for i in range(sample_cases):
        print("sample ", i)
        original = noisy_strand_gen.random_string(bases, strand_length)
        print("original strand: ", original)
        cluster, rev_cluster = noisy_strand_gen.generate_strand_error(original, strand_num, strand_length, choices,
                                                                      bases)

        f = open("tests/test" + str(i) + ".txt", "w+")
        f.write(original + "\n")

        for strand in cluster:
            f.write(strand + "\n")
        f.close()


if __name__ == "__main__":
    main()

