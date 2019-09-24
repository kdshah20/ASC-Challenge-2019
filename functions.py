# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 13:56:07 2019

@author: shahk
"""

import numpy as np
import math

def S_Q_mat(E1,E2, nu_12, G12):
    S11 = 1/E1;
    S12 = -nu_12/E1;
    S22 = 1/E2;
    S66 = 1/G12;

    # Reduced compliance and stiffness matrix
    S = np.matrix([[S11, S12, 0], [S12, S22, 0], [0, 0, S66]]);
    Q = S.I;
    return S, Q;

def transMat(Theta, S):
    N = len(Theta);
    T_array = np.empty((N,3,3));
    Qbar_array = np.empty((N,3,3));
    for i in range(N):
        m = math.cos(math.radians(Theta[i]));
        n = math.sin(math.radians(Theta[i]));
        # Transformation matrix
        T = np.matrix([[m**2, n**2, 2*m*n], [n**2, m**2, -2*m*n], [-m*n, m*n, (m**2 - n**2)]]);
        T[np.absolute(T)<1e-10] = 0;
        T_array[i] = T;
        Sbar = T.transpose()*S*T
        Qbar = Sbar.I;
        Qbar_array[i] = Qbar;
        
    return T_array, Qbar_array;


def ABDmat(Qbar_array, Z, N, h):
    
    A = np.zeros((3,3));
    B = np.zeros((3,3));
    D = np.zeros((3,3));
    
    for i in range(N):
        A = A + Qbar_array[i] * h
        B1 = (Qbar_array[i] * (Z[i+1]**2 - Z[i]**2))*0.5
        B = B + B1
        D1 = (Qbar_array[i] * (Z[i+1]**3  - Z[i]**3))*1/3
        D = D + D1
        
    AB = np.concatenate((A,B),axis=0);
    BD = np.concatenate((B, D),axis=0);
    ABD = np.concatenate((AB,BD),axis=1);
    ABD[np.absolute(ABD)<1e-8] = 0;
    
    return ABD;
        
    

def stress123Cal(abd, Nx, T_array, Qbar_array, N):
    
    epsi_x0 = abd[0,0] * Nx;
    epsi_y0 = abd[1,0] * Nx;
    gamma_xy0 = 0;
    strain_xyz = np.matrix([[epsi_x0], [epsi_y0] , [gamma_xy0]]);
    #strain_xyz.shape = (3,1);
    stress123_array = np.empty((N,3,1));
    Lsig1 = np.zeros((N,1));
    Lsig2 = np.zeros((N,1));
    Lsig6 = np.zeros((N,1));

    for i in range(N):
        stress_xyz = Qbar_array[i]*strain_xyz;
        stress_123 = T_array[i] * stress_xyz;
        stress123_array[i] = stress_123;
        Lsig1[i] = stress_123[0];
        Lsig2[i] = stress_123[1];
        Lsig6[i] = stress_123[2];    
        
    return stress123_array, Lsig1, Lsig2, Lsig6;
    

def maxStress(Xt,Xc,Yt,Yc,shear,Lsig1,Lsig2,Lsig6,N,L0):
    
    Rsig1 = np.zeros((N,1));
    Rsig2 = np.zeros((N,1));
    Rsig6 = np.zeros((N,1));
    
    for i in range(N):
        if Lsig1[i] == 0 and Lsig2[i] ==0 and Lsig6[i] == 0:
            continue   
        
        
        if Lsig1[i] > 0:
            Rsig1[i] = Xt/Lsig1[i];
        elif Lsig1[i] < 0:
            Rsig1[i] = Xc/Lsig1[i];
        
        if Lsig2[i] > 0:
            Rsig2[i] = Yt/Lsig2[i];
        elif Lsig1[i] < 0:
            Rsig2[i] = Yc/Lsig2[i];
         
        if Lsig6[i] == 0:
            Rsig6[i] = 0;
        else:
           Rsig6[i] = shear/Lsig6[i];
          
     
    Sf_array = np.concatenate((Rsig1, Rsig2, Rsig6),axis=1);
    Sf_array[Sf_array ==0]= np.nan;
    M = np.nanmin(np.absolute(Sf_array),axis=1);
    Sf_min = np.nanmin(np.absolute(M))
    #print(Sf_min)
    lfail = np.where((np.absolute(M) == Sf_min))
    
   #Find all layers which failed
    Lfail=lfail;
   
    Lstress = np.concatenate((Lsig1*Sf_min, Lsig2*Sf_min, Lsig6*Sf_min), axis=1);
    
    L0stress = Lstress[L0]
    #print(L0stress);
    
    if L0stress[0,0] > 0:
        Xt = Xt - L0stress[0,0];
    elif L0stress[0,0] < 0:
        Xc = Xc + L0[0,0];
    
    
    if L0stress[0,1] > 0:
        Yt = Yt - L0stress[0,1];
    elif L0stress[0,1] < 0:
        Yc = Yc + L0stress[0,1];
    
    if L0stress[0,2] > 0:
        shear = shear - L0stress[0,2];
    elif L0stress[0,2] < 0:
        shear = shear + L0stress[0,2];
    
    return Sf_min, Lfail, Lstress, Xt, Xc, Yt, Yc, shear;
        
def ranker(popsize,load):
    fit=[0]*(popsize+1)
    fitness=[0]*popsize
    fitn=0
    for i in range(0,popsize):
        fitness[i]=load[i]/sum(load)    #normalize the loads   
        fitn+=fitness[i]
        fit[i+1]=fitn*100               #setting fitness range for selection of parents
    return fit
    
def new_generation(popsize,fit,totallayup,N,layups):
    newlayup = [];
    #selection of parents, random chance, higher fitness has better chance to be selected
    for k in range(0,popsize):
        p1c=(random.randrange(0, 10000, 1))/100
        i=0;chk=0;
        while chk==0:
            if fit[i]<p1c<fit[i+1]:
                p1=i
                chk=1
            i+=1
        p2c=(random.randrange(0, 10000, 1))/100
        i=-1;chk=0;
        while chk==0:
            i+=1
            if fit[i]<p2c<fit[i+1]:
                p2=i
                #make sure parent 2 is not the same as parent 1
                if p2==p1:
                    p2c=(random.randrange(0, 10000, 1))/100
                    i=-1;
                else:
                    chk=1
    
#crossover and mutation
        muta=5 #5%chance
        c1=totallayup[1]
        while c1 in totallayup:
            #randomly select a crossover point
            cpoint=(random.randrange(0, N, 1))
            c1t=totallayup[p1];c2t=totallayup[p2];
            c1=c1t[0:4*cpoint]+c2t[4*cpoint:4*N]                #mix the 2 parents to create an offspring
            lup=[]
            for i in range(0,N):
                ang=int(c1[i*4:(i+1)*4],2)
                lup.append(-90+ang*15)
                if random.randrange(1, 100, 1)<=muta:           #chance to mutate a random ply to a different angle
                    mut=random.randrange(1, 12, 1)
                    lup[i]=-90+mut*15
                    ct=c1
                    c1=ct[0:4*(i)] + '{0:04b}'.format(mut) + ct[4*(i+1):4*N]
        layups.append(lup)
        #print(lup)
        totallayup.append(c1)
        newlayup.append(lup)
    return newlayup

def lstrength(mat):
    
    # Lamina properties
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
    
    [S, Q] = S_Q_mat(E1,E2,nu_12,G12);
    
    
    
    N=8
    Z = np.zeros((N+1,1));
    l=1
    for k in range(0,N+1):
        Z[k] = -H/2 + (l-1)*h
        l+=1

    Theta = np.array(mat);
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
        #print(det)
        if det != 0:
            abd = np.linalg.inv(ABD);
            
            stress123_array, Lsig1,Lsig2,Lsig6 = stress123Cal(abd,Nx,T_array,Qbar_array,N);
        #print(Lsig1, Lsig2, Lsig6)
            Sf_min, Lfail, Lstress, Xt, Xc, Yt, Yc, shear = maxStress(Xt,Xc,Yt,Yc,shear,Lsig1,Lsig2,Lsig6,N,L0);
            for i in range(len(Lfail)):
                Qbar_array[Lfail[i]]=[0];
            
            Nt = temp + Sf_min;
            temp = Nt;
            
            #print(abd)
        else:
            Xt=0;
        
        return temp/H/(1e9);
        

    


        