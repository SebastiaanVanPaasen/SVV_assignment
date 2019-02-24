# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:58:20 2019

@author: thiba
"""

from math import *
import numpy as np

'Needed parameters'
Vy = 37.9*10**3 #N
Izz = 8.40796188796e-05 #m4
boomarea = [9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-053, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05, 9.97962448e-05] #m2
boomy = [0.0, 0.3, 0.5, 0.4, 0.35, 0.3, 0.233, 0.166, 0.10, 0.05, -0.05, -0.1, -0.166, -0.233, -0.30, -0.35, -0.4, -0.5, -0.3]
AI = 2.0
AII = 8.0
h = 17.3*10**(-2) #m
Ca = 0.484 #m
tskin = 1.2*10**(-3) #m
tspar = 2.4*10**(-3) #m
G = 28.0 #GPa
x = 0.1

'Formulas'
LI = 2.0*pi*(h/2.0)*(1.0/8.0)
LII = sqrt((h/2.0)**2 + (Ca-h/2.0)**2)

f = -(Vy/Izz)

#QB
qb23   = 0.0
qb12   = qb23   + f*boomarea[1]*boomy[1]
qb119  = qb12   + f*boomarea[0]*boomy[0]
qb1819 = qb119  + f*boomarea[18]*boomy[18]
qb318  = qb1819 + f*boomarea[17]*boomy[17]

qb34   = 0.0
qb318  = qb34   + f*boomarea[2]*boomy[2]
qb1718 = qb318  + f*boomarea[17]*boomy[17]
qb1617 = qb1718 + f*boomarea[16]*boomy[16]
qb1516 = qb1617 + f*boomarea[15]*boomy[15]
qb1415 = qb1516 + f*boomarea[14]*boomy[14]
qb1314 = qb1415 + f*boomarea[13]*boomy[13]
qb1213 = qb1314 + f*boomarea[12]*boomy[12]
qb1112 = qb1213 + f*boomarea[11]*boomy[11]
qb1011 = qb1112 + f*boomarea[10]*boomy[10]
qb910  = qb1011 + f*boomarea[9]*boomy[9]
qb89   = qb910  + f*boomarea[8]*boomy[8]
qb78   = qb89   + f*boomarea[7]*boomy[7]
qb67   = qb78   + f*boomarea[6]*boomy[6]
qb56   = qb67   + f*boomarea[5]*boomy[5]
qb45   = qb56   + f*boomarea[4]*boomy[4]

#QS0
intqb1 = ((qb318*h)/(G*tspar)) + 2.0*((qb12*LI)/(G*tskin))
intqb2 = ((qb318*h)/(G*tspar)) + ((LII)/(G*tskin))*(2.0*qb45+2.0*qb56+2.0*qb67+2.0*qb78+2.0*qb89+2.0*qb910+qb1011)

mtx1 = np.array([[((LI/(G*tskin))+(h/(G*tspar))), -h/(G*tspar), -2.0*AI],
        [-h/(G*tspar), ((2*LII/(G*tskin))+(h/(G*tspar))), -2.0*AII],
        [2.0*AI, 2.0*AII, 0.0]])
        
mtx3 = np.array([[-intqb1],
        [-intqb2],
        [Vy*Ca - 2*qb12*LI*(Ca-x)]])

mtxsol = np.linalg.solve(mtx1, mtx3)

#Solve every part of airfoil
q12   = qb12   + mtxsol[0]
q23   = qb23   + mtxsol[0]
q34   = qb34   + mtxsol[1]
q45   = qb45   + mtxsol[1]
q56   = qb56   + mtxsol[1]
q67   = qb67   + mtxsol[1]
q78   = qb78   + mtxsol[1]
q89   = qb89   + mtxsol[1]
q910  = qb910  + mtxsol[1]
q1011 = qb1011 + mtxsol[1]
q1112 = qb1112 + mtxsol[1]
q1213 = qb1213 + mtxsol[1]
q1314 = qb1314 + mtxsol[1]
q1415 = qb1415 + mtxsol[1]
q1516 = qb1516 + mtxsol[1]
q1617 = qb1617 + mtxsol[1]
q1718 = qb1718 + mtxsol[1]
q1819 = qb1819 + mtxsol[0]
q119  = qb119  + mtxsol[0]
q318  = qb318  + mtxsol[0] - mtxsol[1]

print q318