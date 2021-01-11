import numpy as np
import sys
import subprocess

idx = int(sys.argv[1])
MAX = int(sys.argv[2])
# idx = 0

etas = np.linspace(0,2,100)
mus = np.array([0.2,0.3])
etas = np.repeat(etas,len(mus))
mus = np.tile(mus,len(etas))

etas = etas[idx::MAX]
mus = mus[idx::MAX]

for i in range(len(etas)):
        eta = etas[i]
        mu = mus[i]

        print(eta)
        bashCommand = f"./ed.exe {eta} {mu}"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        print(output)
