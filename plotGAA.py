import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from  datetime import  datetime


startTime=datetime.now()
R = 160  # small time step number

# Q = 1000#large time step number
dt = 0.000125  # small time step
ds = R * dt  # large time step
timeAxisParts = 10
N = 3000
mu0 = -6
t0 = 1
d0 = -1
mu1 = -6

t1 = t0
d1 = d0
lmd = 6
folder = "/home/disk2/Documents/cppCode/kerr2Kitaev2/testkNum11/kNum"
inGName = folder + str(N) + "G" + "mu0" + str(
    mu0) + "t0" + str(t0) + "d0" + str(d0) + "mu1" + str(mu1) + "t1" + str(t1) + "d1" + str(d1) + "lmd" + str(
    lmd) + ".csv"

inGMat = pd.read_csv(inGName, header=None)
rowN, colN = inGMat.shape
Q = colN - 1

sPoints = [t * ds for t in range(0, colN)]

cutoff = 1


def jumpCut(incr):
    tmp = incr / np.pi
    if tmp > cutoff:
        return incr - 2 * np.pi
    elif tmp < -cutoff:
        return incr + 2 * np.pi
    else:
        return incr


betaMat = pd.DataFrame(index=np.arange(colN), columns=np.arange(rowN - 1))
for q in range(0, Q + 1):
    for k in range(0, N):
        betaMat.iloc[q, k] = jumpCut(inGMat.iloc[k + 1, q] - inGMat.iloc[k][q])



halfBeta=betaMat.iloc[:,0:int(N/2)]
W=[sum(halfBeta.iloc[k,:])/(2*np.pi) for k in range(0, colN)]

maxW=max(W)
minW=min(W)
wticks = []
for i in range(int(np.floor(minW)), int(np.ceil(max(W)))):
    wticks.append(i)


plt.figure()
if maxW - minW > 1:
    plt.yticks(wticks)


tTickStep = (colN - 1) * ds / timeAxisParts
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
outPicFileName = folder + str(N) + "Gwn" + "mu0" + str(
    mu0) + "t0" + str(t0) + "d0" + str(d0) + "mu1" + str(mu1) + "t1" + str(t1) + "d1" + str(d1) + "lmd" + str(
    lmd) + ".png"

plt.savefig(outPicFileName)
plt.close()

endTime=datetime.now()
print("Time: ", str(endTime-startTime))