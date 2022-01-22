import os
import numpy as np

mainproject = "OneHole"  #Set to zero if only one project.
project = "TwoSites"
description = "Testing."
jobname = "myjob"
time = "0:10:00"
runmin = 0
runmax = 0
runsame = 0
nruns = (runmax-runmin) + 1
NICE = 11

#LIBS#
#BOOST = 0       #Higher precision in Eigen-calculations. Time-consuming. Not implemented for now.

#LATTICE#
L        = 2*np.ones(nruns, int)
nruns = len(L)
runmax = runmin + (nruns-1)

Nh = 2*np.ones(nruns, int);

#EXCHANGE#
tl     = 1*np.ones(nruns)
tr     = 1*np.ones(nruns)
Jzl    = -1*np.ones(nruns)
Jzr    = Jzl
Jpml   = -1*np.ones(nruns) #-1*np.logspace(0, np.log10(2), nruns)+np.ones(nruns)
Jpmr   = Jpml

EIGVECS = 1         #Compute eigenvectors?
CORR = 1            #Compute correlations?

RESETOLDFILES = 1



if mainproject:
    os.system("mkdir " + "Data/" + mainproject)
    totalproject = mainproject + "/" + project
else:
    totalproject = project

os.system("mkdir " + "Data/" + totalproject)

descfile = open("Data/" + totalproject + "/description.txt", 'w')
descfile.write(description)

runsub = open("Data/" + totalproject + "/run.sub", "w")

runsub.write("#!/bin/sh\n")
runsub.write("#SBATCH --account=nn4563k\n")
runsub.write("#SBATCH --time=" + time + "\n")
runsub.write("#SBATCH --mem-per-cpu=6999M\n")
runsub.write("#SBATCH --ntasks=1\n")
runsub.write("#SBATCH --signal=B:USR1@60\n")

if runsame:
    runsub.write("#SBATCH --array={}-{}%{}\n".format(runmin, runmax, runsame))
else:
    runsub.write("#SBATCH --array={}-{}\n".format(runmin, runmax))
runsub.write("#SBATCH --job-name="+jobname+"\n")
runsub.write("#SBATCH --error=stderr.dat\n")
runsub.write("#SBATCH --output=stdout.dat\n")
runsub.write("#set -o errexit\n")
runsub.write("module purge\n")
runsub.write("module load intel/2019b\n")
runsub.write('#cleanup "rsync -av $SCRATCH/ $SUBMITDIR/ --exclude=stdout.dat --exclude=stderr.dat"\n')
runsub.write("#cd $SUBMITDIR\n")
runsub.write("#rsync -av $SUBMITDIR/ $SCRATCH/ --exclude=rundir\n")
runsub.write("#cd $SCRATCH\n")
runsub.write("echo Running program.....\n")
runsub.write("$HOME/Documents/ZRHeisenberg/Code/program " + totalproject + " $SLURM_ARRAY_TASK_ID\n")





for run in range(runmin, nruns + runmin):
    #delta = np.ones(num)*(np.logspace(dmin, dmax, num)[run])
    if run < 10:
        os.system("mkdir " + "Data/" + totalproject + "/Run00" + str(run))
        outfile = open("Data/" + totalproject + "/Run00" + str(run) + "/parameters.txt", 'w')
    elif run < 100:
        os.system("mkdir " + "Data/" + totalproject + "/Run0" + str(run))
        outfile = open("Data/" + totalproject + "/Run0" + str(run) + "/parameters.txt", 'w')
    elif run < 1000:
        os.system("mkdir " + "Data/" + totalproject + "/Run" + str(run))
        outfile = open("Data/" + totalproject + "/Run" + str(run) + "/parameters.txt", 'w')
    else:
        print("generate project: Run number is bigger than 999")
        exit(1)

    outfile.write("L = ")
    outfile.write(str(L[run-runmin]))
    outfile.write("\n")
    outfile.write("Nh = ")
    outfile.write(str(Nh[run-runmin]))
    outfile.write("\n")
    outfile.write("\n")

    outfile.write("tl = ")
    outfile.write(str(tl[run-runmin]))
    outfile.write("\n")
    outfile.write("tr = ")
    outfile.write(str(tr[run-runmin]))
    outfile.write("\n")
    outfile.write("Jzl = ")
    outfile.write(str(Jzl[run-runmin]))
    outfile.write("\n")
    outfile.write("Jzr = ")
    outfile.write(str(Jzr[run-runmin]))
    outfile.write("\n")
    outfile.write("Jpml = ")
    outfile.write(str(Jpml[run-runmin]))
    outfile.write("\n")
    outfile.write("Jpmr = ")
    outfile.write(str(Jpmr[run-runmin]))
    outfile.write("\n")
    outfile.write("\n")

    outfile.write("EIGVECS = ")
    outfile.write(str(EIGVECS))
    outfile.write("\n")

    outfile.write("CORR = ")
    outfile.write(str(CORR))
    outfile.write("\n")

    outfile.write("RESETOLDFILES = ")
    outfile.write(str(RESETOLDFILES))
    outfile.write("\n")

    #outfile.write(str(BOOST))


for i in range(runmin, runmax+1):
    if mainproject:
        shellfile = open('Jobs/' + mainproject + project + str(i) + '.sh', 'w')
    else:
        shellfile = open('Jobs/' + project + str(i) + '.sh', 'w')
    shellfile.write('#!/bin/bash\n')
    shellfile.write('nice -' + str(NICE) + ' ~/Documents/ZRHeisenberg/Code/program ' + totalproject + ' ' + str(i) + ' &\n')
    if mainproject:
        shellfile.write('rm ' + '~/Documents/ZRHeisenberg/Jobs/' + mainproject + project + str(i) + '.sh')
    else:
        shellfile.write('rm ' + '~/Documents/ZRHeisenberg/Jobs/' + project + str(i) + '.sh')
    shellfile.close()
    if mainproject:
        os.system('chmod +x Jobs/' + mainproject + project + str(i) + '.sh')
    else:
        os.system('chmod +x Jobs/' + project + str(i) + '.sh')
