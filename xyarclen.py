# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 15:38:58 2019

@author: sseremet
"""
from numpy import *
import matplotlib.pyplot as plt

xlist=[]
ylist=[]
def xyarclen(r,s,ca):
    if s>=0 and s<=pi*r/2:
        x=r*(1-cos(s/r))
        y=r*sin(s/r)
    if s>pi*r/2:
        x=r+(s-pi*r/2)*cos(arctan(r/(ca-r)))
        y=r-(s-pi*r/2)*sin(arctan(r/(ca-r)))
    return x,y

#Only calculates positive y
#for s in range(0,pi*17.3/4+sqrt((48.4-17.3/2)**2+(17.3/2)**2)):
for s in range(0,54):
    a,b=xyarclen(17.3/2,s,0.484)
    xlist.append(a)
    ylist.append(b)
plt.plot(xlist,ylist)