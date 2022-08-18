#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 16:55:06 2022

@author: BJP
"""

import numpy as np
import itertools
import random

#c[0]=1 S, c[1 ]=1 e1, c[2 ]=1 e2, c[3 ]=1 e3, c[4 ]=1 e4, c[5 ]=1 e5, c[6 ]=1 e6, c[7 ]=1 b1, c[8 ]=1 b2, c[9 ]=1 z1, c[10]=1 z2, 
#          c[11]=S e1, c[12]=S e2, c[13]=S e3, c[14]=S e4, c[15]=S e5, c[16]=S e6, c[17]=S b1, c[18]=S b2, c[19]=S z1, c[20]=S z2,  
#                      c[21]=e1e2, c[22]=e1e3, c[23]=e1e4, c[24]=e1e5, c[25]=e1e6, c[26]=e1b1, c[27]=e1b2, c[28]=e1z1, c[29]=e1z2,
#                                  c[30]=e2e3, c[31]=e2e4, c[32]=e2e5, c[33]=e2e6, c[34]=e2b1, c[35]=e2b2, c[36]=e2z1, c[37]=e2z2,
#                                              c[38]=e3e4, c[39]=e3e5, c[40]=e3e6, c[41]=e3b1, c[42]=e3b2, c[43]=e3z1, c[44]=e3z2, 
#                                                          c[45]=e4e5, c[46]=e4e6, c[47]=e4b1, c[48]=e4b2, c[49]=e4z1, c[50]=e4z2, 
#                                                                      c[51]=e5e6, c[52]=e5b1, c[53]=e5b2, c[54]=e5z1, c[55]=e5z2,
#                                                                                  c[56]=e6b1, c[57]=e6b2, c[58]=e6z1, c[59]=e6z2,  
#                                                                                              c[60]=b1b2, c[61]=b1z1, c[62]=b1z2,
#                                                                                                          c[63]=b2z1, c[64]=b2z2,                                                                                                                 
#                                                                                                                      c[65]=z1z2,

                                                                                      

def ChiralSecs(c):   
    #CHIRAL MASSLESS STATES
    
    TRU1=0
    TRU2=0
    TRU3=0
    #--------------------------------------------------------------------------------------------------------  
    ### x phases
    One_x=(1+c[0]+c[1]+c[2]+c[3]+c[4]+c[5]+c[6]+c[9]+c[10])%2
    S_x= (1+c[11]+c[12]+c[13]+c[14]+c[15]+c[16]+c[19]+c[20])%2
    e1_x=(c[11]+1+c[21]+c[22]+c[23]+c[24]+c[25]+c[28]+c[29])%2
    e2_x=(c[12]+1+c[21]+c[30]+c[31]+c[32]+c[33]+c[36]+c[37])%2
    e3_x=(c[13]+1+c[22]+c[30]+c[38]+c[39]+c[40]+c[43]+c[44])%2
    e4_x=(c[14]+1+c[23]+c[31]+c[38]+c[45]+c[46]+c[49]+c[50])%2
    e5_x=(c[15]+1+c[24]+c[32]+c[39]+c[45]+c[51]+c[54]+c[55])%2
    e6_x=(c[16]+1+c[25]+c[33]+c[40]+c[46]+c[51]+c[58]+c[59])%2
    b1_x=(1+c[7]+c[17]+c[26]+c[34]+c[41]+c[47]+c[52]+c[56]+c[61]+c[62])%2
    b2_x=(1+c[8]+c[18]+c[27]+c[35]+c[42]+c[48]+c[53]+c[57]+c[63]+c[64])%2
    z1_x=(c[19]+c[28]+c[36]+c[43]+c[49]+c[54]+c[58]+c[65])%2
    z2_x=(c[20]+c[29]+c[37]+c[44]+c[50]+c[55]+c[59]+c[65])%2
    
    #--------------------------------------------------------------------------------------------------------           
    #b3 phases
    e5_b3=(e5_x+c[52]+c[53])%2
    e6_b3=(e6_x+c[56]+c[57])%2
    z1_b3=(z1_x+c[61]+c[63])%2
    z2_b3=(z2_x+c[62]+c[64])%2
    b1_b3=(b1_x+c[7]+c[60])%2
    S_b3=(S_x+c[17]+c[18])%2
    e3_b3=(e3_x+c[41]+c[42])%2
    e4_b3=(e4_x+c[47]+c[48])%2
    b3_x=(1+b1_x+b2_x+1+One_x)%2 #leave the 1+1 in for now
    #--------------------------------------------------------------------------------------------------------  
    #SPINORIALS
    #--------------------------------------------------------------------------------------------------------  
    #F1_pqrs=S+b1+pe3+qe4+re5+se6
    # Projectors: e1,e2,z1,z2
    subsec=list(itertools.product([0, 1], repeat=4))
    for pqrs in subsec:
        p=pqrs[0]
        q=pqrs[1]
        r=pqrs[2]
        s=pqrs[3]
        
        F1e1=(c[11]+c[26]+p*c[22]+q*c[23]+r*c[24]+s*c[25])%2
        if F1e1==1:
            
            F1e2=(c[12]+c[34]+p*c[30]+q*c[31]+r*c[32]+s*c[33])%2
            if F1e2==1:
                F1z1=(c[19]+c[61]+p*c[43]+q*c[49]+r*c[54]+s*c[58])%2
                if F1z1==1:
                    F1z2=(c[20]+c[62]+p*c[44]+q*c[50]+r*c[55]+s*c[59])%2
                    if F1z2==1:
                        F1chi= (p+q+ c[0] + c[17] + p*c[13] + q*c[14] + c[15] + c[16] \
                        + c[18] + c[60] + p*c[42] + q*c[48] + r*c[53] + s*c[57] \
                        + S_x + b1_x + p*e3_x + q*e4_x + r*e5_x + s*e6_x \
                        + (1-r)*c[52] + (1-r)*p*c[39] + (1-r)*q*c[45] \
                        + (1-s)*c[56] + (1-s)*p*c[40] + (1-s)*q*c[46])%2
                        
                        if F1chi==1:
                            TRU1+=0.5
                        elif F1chi==0:
                            TRU1-=0.5
    #--------------------------------------------------------------------------------------------------------      
    #F2_pqrs=S+b2+pe1+qe2+re5+se6
    # Projectors: e3,e4,z1,z2
    for pqrs in subsec:
        p=pqrs[0]
        q=pqrs[1]
        r=pqrs[2]
        s=pqrs[3]
        
        F2e3=(c[13]+c[41]+p*c[22]+q*c[30]+r*c[39]+s*c[40])%2
        if F2e3==1:
            
            F2e4=(c[14]+c[47]+p*c[23]+q*c[31]+r*c[45]+s*c[46])%2
            if F2e4==1:
                F2z1=(c[19]+c[63]+p*c[28]+q*c[36]+r*c[54]+s*c[58])%2
                if F2z1==1:
                    F2z2=(c[20]+c[64]+p*c[29]+q*c[37]+r*c[55]+s*c[59])%2
                    if F2z2==1:
                        F2chi= (p+q+ c[0] + c[18] + p*c[11] + q*c[12] + c[15] + c[16] \
                        + c[17] + c[60] + p*c[26] + q*c[34] + r*c[52] + s*c[56] \
                        + S_x + b2_x + p*e1_x + q*e2_x + r*e5_x + s*e6_x \
                        + (1-r)*c[53] + (1-r)*p*c[24] + (1-r)*q*c[32] \
                        + (1-s)*c[57] + (1-s)*p*c[25] + (1-s)*q*c[33])%2
                        
                        if F2chi==1:
                            TRU2+=0.5
                        elif F2chi==0:
                            TRU2-=0.5
    #--------------------------------------------------------------------------------------------------------      
    #F3_pqrs=S+b3+pe1+qe2+re3+se4
    # Projectors: e5,e6,z1,z2
    for pqrs in subsec:
        p=pqrs[0]
        q=pqrs[1]
        r=pqrs[2]
        s=pqrs[3]
        
        F3e5=(c[15]+e5_b3+p*c[24]+q*c[32]+r*c[39]+s*c[45])%2
        if F3e5==1:
            
            F3e6=(c[16]+e6_b3+p*c[25]+q*c[33]+r*c[40]+s*c[46])%2
            if F3e6==1:
                F3z1=(c[19]+z1_b3+p*c[28]+q*c[36]+r*c[43]+s*c[49])%2
                if F3z1==1:
                    F3z2=(c[20]+z2_b3+p*c[29]+q*c[37]+r*c[44]+s*c[50])%2
                    if F3z2==1:
                        F3chi= (p+q+ c[0] + S_b3 + p*c[11] + q*c[12] + c[13] + c[14] \
                        + c[17] + b1_b3 + p*c[26] + q*c[34] + r*c[41] + s*c[47] \
                        + S_x + b3_x + p*e1_x + q*e2_x + r*e3_x + s*e4_x \
                        + (1-r)*e3_b3 + (1-r)*p*c[22] + (1-r)*q*c[30] \
                        + (1-s)*e4_b3 + (1-s)*p*c[23] + (1-s)*q*c[31])%2
                        
                        if F3chi==1:
                            TRU3+=0.5
                        elif F3chi==0:
                            TRU3-=0.5
    #--------------------------------------------------------------------------------------------------------      
    #F4_pqrs=S+b1+x+z1+pe3+qe4+re5+se6
    # Projectors: e1,e2,z2

    #--------------------------------------------------------------------------------------------------------      
    #F5_pqrs=S+b2+x+z1+pe1+qe2+re5+se6
    # Projectors: e3,e4,z2

    #--------------------------------------------------------------------------------------------------------      
    #F6_pqrs=S+b3+x+z1+pe1+qe2+re3+se4
    # Projectors: e5,e6,z2
    
    #--------------------------------------------------------------------------------------------------------      
    #F7_pqrs=S+b1+x+z2+pe3+qe4+re5+se6
    # Projectors: e1,e2,z1

    #--------------------------------------------------------------------------------------------------------      
    #F8_pqrs=S+b2+x+z2+pe1+qe2+re5+se6
    # Projectors: e3,e4,z1

    #--------------------------------------------------------------------------------------------------------      
    #F9_pqrs=S+b3+x+z2+pe1+qe2+re3+se4
    # Projectors: e5,e6,z1
    
    #BOSONIC SPARTNER SECTORS:... 
    
    #--------------------------------------------------------------------------------------------------------
    #VECTORIALS - 
    #--------------------------------------------------------------------------------------------------------
    #V1_pqrs=S+b1+x+pe3+qe4+re5+se6
    # Projectors: e1,e2,z1,z2

    #--------------------------------------------------------------------------------------------------------      
    #V2_pqrs=S+b2+x+pe1+qe2+re5+se6
    # Projectors: e3,e4,z1,z2

    #--------------------------------------------------------------------------------------------------------      
    #V3_pqrs=S+b3+x+pe1+qe2+re3+se4
    # Projectors: e5,e6,z1,z2
    
    
    return [TRU1,TRU2,TRU3]


import multiprocessing as mp
import csv
import os

# The following code will read the tachyon free models in from the file produced by TachCheckSO10.py and quickly computes the traces 

def matrixReader(): #this is an efficient reader of the tach free matrix file that doesn't store all matrices in memory 
    with open("TachFreeSO10Batch.csv", newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        i = 0
        acc = []
        for line in reader:
            acc.append([float(x) for x in line])
            #acc=np.array(acc)
            i+=1
            if (i % 12 == 0):
                yield np.array(acc)
                acc = []

def getTraces(GGSOs):#this converts the upper triangle entries into the form used here and applies the ChiralSecs fn to get the U(1) traces
    
    Gstrip=list(GGSOs[np.triu_indices(12, k=1)])
    GGSOs=[1 if phase==-1 else 0 for phase in Gstrip ]
    traces=ChiralSecs(GGSOs)
    
    return traces


if __name__=='__main__':
    
    cpu_count = os.cpu_count()
    #print("cpu count: ", cpu_count)


    with mp.Pool(cpu_count) as p:
        TRes1 = p.map(getTraces,matrixReader(),chunksize=100) # this is just some multiprocessing thing where the getTraces fn is mapped onto each of the matrices from the file one by one across the CPU cores...
        #for x in p.imap(consttermCheck, matrixes(), chunksize=100):
        for traces in TRes1:
            print(traces)
    
    
                
            
    
    
