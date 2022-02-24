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

for run_number in run_numbers:
    if float(run_number) < 10:
        run_number = "00" + run_number
    elif float(run_number) < 100:
        run_number = "0" + run_number
    elif float(run_number) < 1000:
        run_number = run_number
    else:
        print("Run number too big")
        exit(1)


    #data = get_data("Run" + run_number + "/eigvals.txt", ["", "OP1", "OP2", "OP1norm", "OP2norm"])

    eigvals = []

    infile = open("Run" + run_number + "/GSHoleDensity.txt", "r")

    GS = 1e9

    for line in infile:
        words = line.split()
        if float(words[1]) < GS:
            GS = float(words[1])
            HoleDens = [float(i) for i in words[2:]]

    plt.plot(HoleDens, label = "Run " + run_number)
    plt.xlabel(r"site $j$", fontsize=14)
    plt.ylabel(r"Hole density", fontsize=14)
plt.legend()
plt.show()
