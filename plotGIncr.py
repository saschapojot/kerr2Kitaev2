import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

R = 160  # small time step number

# Q = 1000#large time step number
dt = 0.000125  # small time step
ds = R * dt  # large time step
timeAxisParts = 10
N = 250
mu0 = -6
t0 = 1
d0 = -1
mu1 = -6

t1 = t0
d1 = d0
lmd = 6
folder="/home/disk2/Documents/cppCode/kerr2Kitaev2/testkNum251/kNum"

inGName = folder+ str(N) + "G" + "mu0" + str(
    mu0) + "t0" + str(t0) + "d0" + str(d0) + "mu1" + str(mu1) + "t1" + str(t1) + "d1" + str(d1) + "lmd" + str(
    lmd) + ".csv"

inGMat=pd.read_csv(inGName,header=None)

rowN,colN=inGMat.shape
sPoints=[t*ds for t in range(0, colN)]
W=[(inGMat.iloc[int(rowN/2)-1,t]-inGMat.iloc[0,t])/(2*np.pi) for t in range(0,colN)]
titleStr = "$\mu_{0}=$" + str(mu0) + ", $t_{0}=$" + str(t0) + ", $\Delta_{0}$=" + str(d0) + ", $\mu_{1}=$" + str(
    mu1) + ", $t_{1}=$" + str(t1) + ", $\Delta_{1}=$" + str(d1) + ", $\lambda=$" + str(lmd) + ", $N=$" + str(N)
plt.figure()
plt.plot(sPoints,W,color="black")
outPicFileName = folder +"noCorr"+ str(N) + "wn" + "mu0" + str(
    mu0) + "t0" + str(t0) + "d0" + str(d0) + "mu1" + str(mu1) + "t1" + str(t1) + "d1" + str(d1) + "lmd" + str(
    lmd) + ".png"
maxW = max(W)
minW = min(W)
wticks = []
for i in range(int(np.floor(minW)), int(np.ceil(max(W)))):
    wticks.append(i)

plt.figure()
if maxW - minW > 1:
    plt.yticks(wticks)
plt.plot(sPoints,W,color="black")

xTickVals=[]
tTickStep=(colN-1)*ds/timeAxisParts
for i in range(0, timeAxisParts + 1):
    xTickVals.append(i * tTickStep)
plt.xticks(xTickVals)
plt.xlabel("time")
plt.ylabel("winding number")
plt.title(titleStr)
plt.savefig(outPicFileName)
