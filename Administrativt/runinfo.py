import os

print("RUN                L       Nh      tl       tr       Jzl       Jzr       Jpml       Jplr")

directories = [d for d in os.listdir(os.getcwd()) if (d[0:3] == 'Run' and os.path.isdir(d))]
N = len(directories)

folder = 0
folders = 0
while folders < N:
    if folder >= 0 and folder < 10:
        if os.path.exists("Run00" + str(folder)):
            outfile = open("Run00" + str(folder) + "/parameters.txt", 'r')
            words = outfile.readline().split()
            print("Run00" + str(folder) + "     %15s    %4s    %4s    %4s    %+2.2f    %+2.2f    %+2.2f    %+2.2f"\
             % (words[0], words[1], words[2], words[3], float(words[4]), float(words[5]), float(words[6]), float(words[7])))
            folder += 1
            folders += 1
        else:
            folder += 1

    elif folder > 9 and folder < 100:
        if os.path.exists("Run0" + str(folder)):
            outfile = open("Run0" + str(folder) + "/parameters.txt", 'r')
            words = outfile.readline().split()
            print("Run0" + str(folder) + "     %15s    %4s    %4s    %4s    %+2.2f    %+2.2f    %+2.2f    %+2.2f"\
             % (words[0], words[1], words[2], words[3], float(words[4]), float(words[5]), float(words[6]), float(words[7])))
            folder += 1
            folders += 1
        else:
            folder += 1


    elif folder < 1000:
        if os.path.exists("Run" + str(folder)):
            outfile = open("Run" + str(folder) + "/parameters.txt", 'r')
            words = outfile.readline().split()
            print("Run" + str(folder) + "     %15s    %4s    %4s    %4s    %+2.2f    %+2.2f    %+2.2f    %+2.2f"\
             % (words[0], words[1], words[2], words[3], float(words[4]), float(words[5]), float(words[6]), float(words[7])))
            folder += 1
            folders += 1
        else:
            folder += 1



"""
if folders >= 0 and folders < 10:
    while os.path.exists("Run00" + str(folders)):
        outfile = open("Run00" + str(folders) + "/parameters.txt", 'r')
        words = outfile.readline().split()
        print("Run00" + str(folders) + "     %15s    %4s    %4s    %4s    %+2.2f    %+2.2f    %+2.2f    %+2.2f"\
         % (words[0], words[1], words[2], words[3], float(words[4]), float(words[5]), float(words[6]), float(words[7])))
        folders += 1
if folders > 9 and folders < 100:
    while os.path.exists("Run0" + str(folders)):
        outfile = open("Run0" + str(folders) + "/parameters.txt", 'r')
        words = outfile.readline().split()
        print("Run0" + str(folders) + "     %15s    %4s    %4s    %4s    %+2.2f    %+2.2f    %+2.2f    %+2.2f"\
         % (words[0], words[1], words[2], words[3], float(words[4]), float(words[5]), float(words[6]), float(words[7])))
        folders+=1
if folders < 1000:
    while os.path.exists("Run" + str(folders)):
        outfile = open("Run" + str(folders) + "/parameters.txt", 'r')
        words = outfile.readline().split()
        print("Run" + str(folders) + "     %15s    %4s    %4s    %4s    %+2.2f    %+2.2f    %+2.2f    %+2.2f"\
         % (words[0], words[1], words[2], words[3], float(words[4]), float(words[5]), float(words[6]), float(words[7])))
        folders+=1
"""
