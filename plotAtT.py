import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

R = 160  # small time step number
#maybe incr is > 360
# Q = 1000#large time step number
dt = 0.000125  # small time step
ds = R * dt  # large time step

timeAxisParts = 10
N1 = 500
N2=1000
mu0 = 1
t0 = 1
d0 = -1
mu1 = 1

t1 = t0
d1 = d0
lmd = 2.5
prefixFolder="/home/disk2/Documents/cppCode/kerr2Kitaev2/quench15/"
folder=prefixFolder+ "kNum"

inBetaName = folder+ str(N1) + "beta" + "mu0" + str(
    mu0) + "t0" + str(t0) + "d0" + str(d0) + "mu1" + str(mu1) + "t1" + str(t1) + "d1" + str(d1) + "lmd" + str(
    lmd) + ".csv"
inBetaMat = pd.read_csv(inBetaName, header=None)
rowN = len(inBetaMat)

colN = len(inBetaMat.columns)
inBetaName2=folder+ str(N2) + "beta" + "mu0" + str(
    mu0) + "t0" + str(t0) + "d0" + str(d0) + "mu1" + str(mu1) + "t1" + str(t1) + "d1" + str(d1) + "lmd" + str(
    lmd) + ".csv"
inBetaMat2=pd.read_csv(inBetaName2,header=None)


for tRow in range(0,rowN):

  time=tRow*ds
  row1BetaVals=inBetaMat.iloc[tRow,:]
  row2BetaVals=inBetaMat2.iloc[tRow,:]

  row1Num=len(row1BetaVals)
  row2Num=len(row2BetaVals)
  kVals1=[2*np.pi/N1* kTmp for kTmp in range(0,row1Num)]
  kVals2=[2*np.pi /N2* kTmp for kTmp in range(0,row2Num)]
  ph1=[sum(row1BetaVals[0:(k+1)]) for k in range(0,row1Num)]
  ph2=[sum(row2BetaVals[0:(k+1)]) for k in range(0,row2Num)]

  plt.figure()
  plt.plot(kVals1,ph1,color="red",label=str(N1))
  plt.plot(kVals2,ph2,color="black",label=str(N2))
  plt.legend(loc="best")
  plt.xlabel("momentum")
  plt.ylabel("AA phase")
  titleStr = "$\mu_{0}=$" + str(mu0) + ", $t_{0}=$" + str(t0) + ", $\Delta_{0}$=" + str(d0) + ", $\mu_{1}=$" + str(
    mu1) + ", $t_{1}=$" + str(t1) + ", $\Delta_{1}=$" + str(d1) + ", $\lambda=$" + str(lmd)+", time = "+str(time)
  plt.title(titleStr)
  outFolder=prefixFolder+"allPhases2/"
  outPicFileName =outFolder + "cmp" +str(N1)+"and"+str(N2) + "mu0" + str(
    mu0) + "t0" + str(t0) + "d0" + str(d0) + "mu1" + str(mu1) + "t1" + str(t1) + "d1" + str(d1) + "lmd" + str(
    lmd) +"row"+str(tRow)+ ".png"
  plt.savefig(outPicFileName)
  plt.close()