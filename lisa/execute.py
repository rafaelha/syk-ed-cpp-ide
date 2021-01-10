import numpy as np
import sys
import subprocess

idx = int(sys.argv[1])
# idx = 0

etas = np.linspace(0,2,100)
eta = etas[idx]

eta = 0

mus = np.linspace(0,1,30)
mu = mus[idx]

print(eta)
bashCommand = f"./ed.exe {eta} {mu}"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()
print(output)
