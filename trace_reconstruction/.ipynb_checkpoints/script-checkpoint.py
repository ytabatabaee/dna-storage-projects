import subprocess
import configparser
import sys

config_location = sys.argv[1]
config = configparser.ConfigParser()
config.read(config_location)
parameters = config['parameters']

sample_cases = int(parameters['sample_cases'])


def main():
    subprocess.check_call(["python3.6", "test_gen.py", sys.argv[1]])

    processes = []
    for j in range(0, 10):
        for i in range(0, int(sample_cases/10)):
            processes.append(subprocess.Popen(["python3.6", "trace_reconstruction.py", sys.argv[1], str(j*int(sample_cases/10) + i)]))

        for proc in processes:
            proc.communicate()

    subprocess.check_call(["python3.6", "compute_res.py", sys.argv[1]])


if __name__ == "__main__":
    main()

