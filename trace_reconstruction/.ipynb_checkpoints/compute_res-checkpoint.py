import matplotlib.pyplot as plt
import configparser
import sys
import ast
import os

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

file_params = "w" + str(window_size) + "n" + str(strand_num) + "l" + str(strand_length) + "t" + str(
    sample_cases) + "p" + str(int(error_rate * 100))   


def show_results(prob):
    error = [(x) / sample_cases for x in prob]
    print(error)

    file_name = "results/error" + file_params + ".txt"
    if not flags.getboolean('non_dp'):
        file_name = "results/error" + file_params + "_dp" + ".txt"
    f = open(file_name, "w+")
    f.write(str(error))
    f.close()

    temp = list(range(1, strand_length + 1))
    fig, ax = plt.subplots()
    ax.plot(temp, error,
            label="P = " + str((int)(error_rate * 100)) + "%, N=" + str(strand_num) + ",tests=" + str(
                sample_cases) + ",max=" + str(max(error) * 100))
    legend = ax.legend()

    plt.xlabel("Position (1-L)")
    plt.ylabel("Probability of incorrect base")
    if flags.getboolean('non_dp'):
        plt.savefig("results/plot" + file_params + ".png")
    else:
        plt.savefig("results/plot" + file_params + "_dp" + ".png")


def clean_up():
    if flags.getboolean('cleanup_error'):
        if os.path.exists("results/error" + file_params + ".txt"):
            os.remove("results/error" + file_params + ".txt")
        if os.path.exists("results/error" + file_params + "_dp" + ".txt"):
            os.remove("results/error" + file_params + "_dp" + ".txt")

    if flags.getboolean('cleanup_codes'):
        if os.path.exists("dp" + str(strand_num) + ".cpp"):
            os.remove("dp" + str(strand_num) + ".cpp")
        if os.path.exists("dp" + str(strand_num)):
            os.remove("dp" + str(strand_num))
        if os.path.exists("dp" + str(strand_num) + ".exe"):
            os.remove("dp" + str(strand_num) + ".exe")
            
            
def plot_distance_histogram(dtype, hist):
    print(dtype)
    print(hist)
    plt.xlabel(dtype + ' distance')
    plt.ylabel('test count')
    plt.plot(range(len(hist)), hist)
    plt.title('distribution of ' + dtype + ' distances')
    if flags.getboolean('non_dp'):
        plt.savefig("results/plot_" + dtype + file_params + ".png")
    else:
        plt.savefig("results/plot_" + dtype + file_params + "_dp" + ".png")

        
def main():
    prob = [0] * strand_length
    avg_ed = 0
    avg_hd = 0
    ed_hist = [0] * strand_length
    hd_hist = [0] * strand_length
    for i in range(sample_cases):
        with open("outputs/out" + str(i) + ".txt", 'r') as f:
            error = ast.literal_eval(f.readline())
            prob = [a + b for a, b in zip(prob, error)]
            ed = int(f.readline())
            ed_hist[ed] += 1
            hd = int(f.readline())
            hd_hist[hd] += 1
            avg_ed += ed / sample_cases
            avg_hd += hd / sample_cases
            # print(prob)
            # print(len(prob))

    show_results(prob)
    plot_distance_histogram('edit', ed_hist)
    plot_distance_histogram('hamming', hd_hist)
    print(avg_ed)
    print(avg_hd)
    clean_up()


if __name__ == "__main__":
    main()

