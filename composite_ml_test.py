# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 13:46:54 2019

@author: shahk
"""

#Generate random stacking sequences
import random
import numpy as np
import subprocess

randtxt = np.loadtxt(fname = "C:/Users/shahk/Documents/ldata.csv", delimiter=',', unpack=False)
zero = np.zeros((18001,1))
mat = np.concatenate((randtxt,zero),axis=1)

for i in range(len(mat)):
    np.random.shuffle(mat[i, :])

# Calling CLT functions
subprocess.call("C:/Users/shahk/Documents/functions.py", shell=True)

# Lamina properties
E1 = 163e9;
E2 = 10e9;
nu_12 = 0.35;
G12 = 5.17e9;

# Applied Loading
Nx = 1 ; 

# Laminate properties
h = 0.125e-3;  #thickness of each layer
H = 8*h;
N=8
Z = np.zeros((N+1,1));
l=1
for k in range(0,N+1):
    Z[k] = -H/2 + (l-1)*h
    l+=1
    
    
[S, Q] = S_Q_mat(E1,E2,nu_12,G12);

j_array = np.zeros((len(layups),1));
Nt_array = np.zeros((len(layups),1));


    
Theta = np.array([0,45,-45,90,90,-45,45,0])


for j in range(len(layups)):
    Theta = np.array(layups[j]);
    N = len(Theta); # number of layers
    m = (min(abs(Theta)));
    L0 = np.where(Theta == np.absolute(Theta).min())
    

    T_array, Qbar_array = transMat(Theta,S);

    temp=0;
    Xt = 2720e6;
    Xc = 1690e6;
    Yt = 64e6;
    Yc = 200e6;
    shear = 120e6;
    while Xt > 0:
    
        ABD = ABDmat(Qbar_array,Z,N,h);
        det = np.linalg.det(ABD)
        if det != 0:
            abd = np.linalg.inv(ABD);
        else:
            Nt = 0;
            j_array[j] = j;
            Nt_array[j] = 0;
            break
        
    
        stress123_array, Lsig1,Lsig2,Lsig6 = stress123Cal(abd,Nx,T_array,Qbar_array,N);
        Sf_min, Lfail, Lstress, Xt, Xc, Yt, Yc, shear = maxStress(Xt,Xc,Yt,Yc,shear,Lsig1,Lsig2,Lsig6,N,L0);
        for i in range(len(Lfail)):
            Qbar_array[Lfail[i]]=[0];
    
        Nt = temp + Sf_min;
        temp = Nt;
    Nt_array[j] = Nt/H/(1e9);
    #print(Nt)
    #[Sf_min, Lfail, Lstress, Xt, Xc, Yt, Yc, shear] = TsaiWu(Xt,Xc,Yt,Yc,shear,F12_star,Lsig1,Lsig2,Lsig6,N,L0);

 #writematrix(data,'data.csv')

zero_val = np.where(Nt_array == 0)[0] 
mat = np.delete(mat, zero_val, axis=0)
Nt_array = np.delete(Nt_array, zero_val)

data = np.concatenate((mat,(Nt_array/H).reshape(len(Nt_array),1)), axis=1)