import numpy as np
import sys
import subprocess

idx = int(sys.argv[1])
MAX = int(sys.argv[2])
# idx = 0

etas = np.array([1])
mus = np.linspace(0,0.2,80)
etas = np.repeat(etas,len(mus))
mus = np.tile(mus,len(etas))
seeds = np.arange(10)

etas = etas[idx::MAX]
mus = mus[idx::MAX]

for i in range(len(etas)):
        eta = etas[i]
        mu = mus[i]

        print(eta)
        for seed in seeds:
                bashCommand = f"./ed.exe {eta} {mu} {seed}"
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()
                print(output)
        
