#coding: UTF-8

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys

import plotstyle

try:
    run_numbers = sys.argv[1:]
except:
    print("Give run number as command line argument.")
    exit(1)


def get_data(filename, variables):
    df = pd.read_csv(filename,\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df
    #using pandas to read the data files

corrzlist = []
corrpmlist = []

run_number = run_numbers[0]

if float(run_number) < 10:
    run_number = "00" + run_number
elif float(run_number) < 100:
    run_number = "0" + run_number
elif float(run_number) < 1000:
    run_number = run_number
else:
    print("Run number too big")
    exit(1)


#data = get_data("Run" + run_number + "/Partition.txt", ['Tinv' 'Z'])

Zfile = open("Run" + run_number + "/Partition.txt", "r")

Zbeta = []
Z = []

for line in Zfile:
    words = line.split()
    Zbeta.append(float(words[0]))
    Z.append(float(words[1]))

Zbeta = np.array(Zbeta)
Z = np.array(Z)

Zfile.close()

paramsfile = open("Run" + run_number + "/parameters.txt", "r")

for line in paramsfile:
    words = line.split()
    if len(words) > 0:
        if words[0] == "L":
            L = int(words[-1])
            break

print("L = %d" % L)

paramsfile.close()

infile = open("Run" + run_number + "/Corr.txt", "r")

beta = []
time = []

for line in infile:
    words = line.split()
    beta.append(float(words[0]))
    time.append(float(words[1]))

beta = np.unique(np.array(beta)) #May be sorting the elements?
time = np.unique(np.array(time)) #May be sorting the elements?

Nb = len(beta)
Nt = len(time)

infile.close()

corrzreal = np.zeros((2, 2*L, Nb, Nt))
corrzimag = np.zeros((2, 2*L, Nb, Nt))
corrpmreal = np.zeros((2, 2*L, Nb, Nt))
corrpmimag = np.zeros((2, 2*L, Nb, Nt))


for j in range(len(run_numbers)):

    run_number = run_numbers[j]

    if float(run_number) < 10:
        run_number = "00" + run_number
    elif float(run_number) < 100:
        run_number = "0" + run_number
    elif float(run_number) < 1000:
        run_number = run_number
    else:
        print("Run number too big")
        exit(1)

    infile = open("Run" + run_number + "/Corr.txt", "r")

    for line in infile:
        words = line.split()
        b = float(words[0])
        t = float(words[1])

        bind = np.where(beta == b)[0]
        tind = np.where(time == t)[0]

        words = words[2:]
        for i in range(2*L):
            cz = [float(word) for word in words[i][1:-1].split(",")]
            cpm = [float(word) for word in words[2*L+i][1:-1].split(",")]

            corrzreal[j,i,bind,tind] = cz[0]
            corrzimag[j,i,bind,tind] = cz[1]
            corrpmreal[j,i,bind,tind] = cpm[0]
            corrpmimag[j,i,bind,tind] = cpm[1]
        corrzlist.append(corrzreal)
        corrpmlist.append(corrpmreal)

plt.figure()
for b in range(len(beta)):
    plt.plot(range(2*L), 0.5*(corrpmreal[0,:,b,0]-corrpmreal[1,:,b,0]), label=r"$\beta = %.2f$" % beta[b])
plt.xlabel(r"site $j$")
plt.ylabel(r"$\Delta\frac{1}{2}\langle S^x_{0}S^x_{j} + S^y_{0}S^y_{j} \rangle$")
plt.title("One hole, Heisenberg AFM, difference between t = 0 and t = 1")
plt.legend()

plt.figure()
for b in range(len(beta)):
    plt.plot(range(2*L), corrzreal[0,:,b,0]-corrzreal[1][:,b,0], label=r"$\beta = %.2f$" % beta[b])
plt.xlabel(r"site $j$")
plt.ylabel(r"$\Delta\langle S^z_{0}S^z_{j} \rangle$")
plt.title("One hole, Heisenberg AFM, difference between t = 0 and t = 1")
plt.legend()

plt.show()
