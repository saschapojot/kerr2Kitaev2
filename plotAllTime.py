import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

R = 160  # small time step number

# Q = 1000#large time step number
dt = 0.000125  # small time step
ds = R * dt  # large time step
N = 96000
mu0 = 0
t0 = 1
d0 = -1
mu1 = -3

t1 = t0
d1 = d0
lmd = 0
folder = "/home/disk2/Documents/cppCode/kerr2Kitaev2/testkNum25/kNum"
dirOut = "/home/disk2/Documents/cppCode/kerr2Kitaev2/testkNum25/example"+str(N)+"/kNum"
inGName = folder + str(N) + "G" + "mu0" + str(mu0) + "t0" + str(
    t0) + "d0" + str(d0) + "mu1" + str(mu1) + "t1" + str(t1) + "d1" + str(d1) + "lmd" + str(lmd) + ".csv"
inGMat = pd.read_csv(inGName, header=None)
rowN, colN = inGMat.shape

dk = 2 * np.pi / (rowN - 1)
kvals = [k * dk for k in range(0, rowN - 1)]
for q in range(0, colN):
    gThisCol = inGMat.iloc[:, q]
    diffThis = [gThisCol[k + 1] - gThisCol[k] for k in range(0, rowN - 1)]
    plt.figure()
    plt.plot(kvals, diffThis, color="black")
    plt.xlabel("k")
    plt.ylabel("Phase difference")
    titleStr = "$\mu_{0}=$" + str(mu0) + ", $t_{0}=$" + str(t0) + ", $\Delta_{0}$=" + str(d0) + ", $\mu_{1}=$" + str(
        mu1) + ", $t_{1}=$" + str(t1) + ", $\Delta_{1}=$" + str(d1) + ", $\lambda=$" + str(lmd) + ", $N=$" + str(
        N) + ", Time = " + str(q * ds)
    plt.title(titleStr)

    outPicFileName = dirOut + str(N) + "diff" + "mu0" + str(
        mu0) + "t0" + str(t0) + "d0" + str(d0) + "mu1" + str(mu1) + "t1" + str(t1) + "d1" + str(d1) + "lmd" + str(
        lmd) + "q" + str(q) + ".png"
    plt.savefig(outPicFileName)
    plt.close()
