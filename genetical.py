# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 07:14:37 2019

@author: k
"""
import random
import numpy as np
popsize=10;  #size of GA parent population
N=8;        #number of plies
nIter=10   # number of iterations

mat=[0,45,-45,90,90,-45,45,0]


layups = [];
totallayup=[];
load=[];
sindex=[]*popsize;
fitness=[0]*popsize
#create random initial population
for i in range(0,popsize):
    j=0;total="";lup=[]
    while j<N:                              
        ang=random.randrange(1, 12, 1)
        stri='{0:04b}'.format(ang)
        total=total+stri
        if mat[j]+ang > 90:
            lup.append(mat[j]-ang)
        else:
            lup.append(mat[j]+ang)              #append lamina angle
        j=j+1
    layups.append(lup)                      #adds new random sequence to layup[]
    totallayup.append(total)
    l = lstrength(np.array(lup));                     #get peak load from the clt and damage analysis
    load.append(l)                          
print("\nBest loads and stacking sequences")    
for z in range(0,nIter):                                                       
    sindex=sorted(range(len(load)), key=lambda i: load[i])[-popsize:]           #sorted index of load[] for best 5 laminates
    sload=[];bestlayups=[];
    print("\nIteration Number"+str(z+1))
    print("LOAD      STACKING SEQUENCE")
    for x in sindex: 
        sload.append(load[x])
        bestlayups.append(layups[x])
        print("{:.1f}".format(load[x]), layups[x])
    fit = ranker(popsize,sload);                                    #evaluates fitness of top 5 laminates
    newlayup = new_generation(popsize,fit,totallayup,N,layups);     #get new generation of stacking sequences based on the top 5 laminates
    for i in range(0,popsize):
        lup2=newlayup[i]
        l2 = lstrength(np.array(lup2));     #get peak load from the clt and damage analysis for new laminates
        #print(l2,lup2)
        
        load.append(l2)
#Rank

