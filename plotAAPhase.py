import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
R = 160  # small time step number
#maybe incr is > 360
# Q = 1000#large time step number
dt = 0.000125  # small time step
ds = R * dt  # large time step
timeAxisParts = 10
N = 96000
mu0 = 1
t0 = 1
d0 = -1
mu1 = 1

t1 = t0
d1 = d0
lmd = 2.51
folder="/home/disk2/Documents/cppCode/kerr2Kitaev2/quench13/kNum"
inBetaName = folder+ str(N) + "beta" + "mu0" + str(
    mu0) + "t0" + str(t0) + "d0" + str(d0) + "mu1" + str(mu1) + "t1" + str(t1) + "d1" + str(d1) + "lmd" + str(
    lmd) + ".csv"
inBetaMat = pd.read_csv(inBetaName, header=None)
rowN = len(inBetaMat)

colN = len(inBetaMat.columns)

sPoints = [t * ds for t in range(0, rowN)]
#W = []
# for i in range(0, rowN):
#     tmp = 0
#     for j in range(0, int(colN / 2)):
#         tmp += inBetaMat.iloc[i, j]
#     W.append(tmp / (2 * np.pi))

halfInBetaTab=inBetaMat.iloc[:,0: int(colN / 2)]
W=[sum(halfInBetaTab.iloc[k,:])/(2*np.pi) for k in range(0, rowN)]
maxW = max(W)
minW = min(W)
wticks = []
for i in range(int(np.floor(minW)), int(np.ceil(max(W)))):
    wticks.append(i)

plt.figure()
if maxW - minW > 1:
    plt.yticks(wticks)

tTickStep = (rowN - 1) * ds / timeAxisParts
xTickVals = []
for i in range(0, timeAxisParts + 1):
    xTickVals.append(i * tTickStep)
plt.xticks(xTickVals)
plt.plot(sPoints, W, color="black")
plt.xlabel("time")
plt.ylabel("winding number")
titleStr = "$\mu_{0}=$" + str(mu0) + ", $t_{0}=$" + str(t0) + ", $\Delta_{0}$=" + str(d0) + ", $\mu_{1}=$" + str(
    mu1) + ", $t_{1}=$" + str(t1) + ", $\Delta_{1}=$" + str(d1) + ", $\lambda=$" + str(lmd) + ", $N=$" + str(N)
plt.title(titleStr)
outPicFileName = folder + str(N) + "wn" + "mu0" + str(
    mu0) + "t0" + str(t0) + "d0" + str(d0) + "mu1" + str(mu1) + "t1" + str(t1) + "d1" + str(d1) + "lmd" + str(
    lmd) + ".png"

plt.savefig(outPicFileName)
plt.close()
