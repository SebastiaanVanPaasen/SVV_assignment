# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 15:38:58 2019

@author: sseremet
"""
from numpy import *
def xyarclen(r,s,ca):
    if s>=0 and s<=pi*r/2:
        x=r*(1-cos(s/r))
        y=r*sin(s/r)
    if s>pi*r/2:
        x=r+(s-pi*r/2)*cos(arctan(r/(ca-r)))
        y=r-(s-pi*r/2)*sin(arctan(r/(ca-r)))
    return x,y

#Only calculates positive y