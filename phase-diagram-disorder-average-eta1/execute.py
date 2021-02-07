import numpy as np
import sys
import subprocess

idx = int(sys.argv[1])
MAX = int(sys.argv[2])
# idx = 0
# MAX = 1

mus = np.loadtxt('missing_mus.txt')
etas = np.loadtxt('missing_etas.txt')
seeds = np.loadtxt('missing_seeds.txt')

etas = etas[idx::MAX]
mus = mus[idx::MAX]
seeds = seeds[idx::MAX]

for i in range(len(etas)):
        eta = etas[i]
        mu = mus[i]
        seed = seeds[i]

        print(eta)
        bashCommand = f"./ed.exe {eta} {mu} {seed}"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        # print(output)
