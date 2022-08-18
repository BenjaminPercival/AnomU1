#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 17:44:38 2022

@author: viktor
"""


import numpy as np
#--------------------------------------------------------------------------------------------------------      
# TACHYON CHECKER- finds SO(10) tach free mats and saves to a file 

One = 0
S =  1
e1 = 2
e2 = 3
e3 = 4
e4 = 5
e5 = 6
e6 = 7
b1 = 8
b2 = 9
z1 = 10
z2 = 11

BP=np.array([
 (-12, 4,  0, 0,  0, 0,  0,  0, -4, -4,  -4, -4),
 (4, 4,  0, 0,  0, 0, 0, 0, 2, 2,  0,  0),
 (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
 (0, 0, 0, 0,  0, 0,  0,  0,  0,  0,  0, 0),
 (0, 0, 0, 0, 0, 0, 0, 0,  0, 0,  0,  0),
 (0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0,  0),
 (0, 0, 0, 0, 0,  0, 0, 0, 0, 0, 0, 0),
 (0, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0),
 (-4, 2, 0, 0, 0, 0, 0, 0, -4, -4, 0, 0),
 (-4, 2, 0, 0, 0, 0, 0, 0, -4, -4, 0, 0),
 (-4, 0, 0, 0,  0, 0, 0, 0,  0,  0, -4,  0),
 (-4, 0, 0, 0,  0, 0, 0, 0,  0,  0,  0, -4)])



def x(bvec,GSO):
    phase=(GSO[bvec][One]+GSO[bvec][S]+GSO[bvec][e1]+GSO[bvec][e2]+GSO[bvec][e3]+GSO[bvec][e4]+GSO[bvec][e5]+GSO[bvec][e6]+GSO[bvec][z1]+GSO[bvec][z2])%2
    return phase

from itertools import combinations 
#from operator import mul
#import timeit

def rSubset(arr, r): 
  
    return list(combinations(arr, r)) 


comb3=rSubset([2,3,4,5,6,7],3)
comb2=rSubset([2,3,4,5,6,7],2)
comb1=rSubset([2,3,4,5,6,7],1)

def tachyonCheckerNew3s(GSO):

    full=[2,3,4,5,6,7]
    TachyonPresent=False
    
    for comb in comb3:
        combLst=list(comb)
        row=np.zeros(8) # S, z1,z2,el,em,en,x IN THAT ORDER
        projes=[x for x in full if x not in combLst]
        projlst=[S,z1,z2]+projes
        #print(projlst)
        #print("for comb: ", comb)
        for proj in projlst:
            #print("projlst.index(proj):", projlst.index(proj))
            #print("for combLst: ", combLst)
            #print("with proj: ", proj)
            phase=[GSO[proj][i] for i in combLst] 
            #print(np.prod(phase))
            row[projlst.index(proj)]=np.prod(phase)
        
        
        row[-1]=np.prod([x(i,GSO) for i in combLst])
        #print(row)
        
        #check vectorial tachyons:
        if row[0]==1:
            if np.count_nonzero(row == -1)==1 and row[3]==1:
                TachyonPresent=True
                #print("vect tach: ", comb)
                #print("row is: ", row)
                break
            
            elif np.count_nonzero(row == -1)==2 and row[3]==row[-1]==-1: #psi45
                TachyonPresent=True
                #print("vect tach: ", comb)
                #print("row is: ", row)
                break
            
            elif np.count_nonzero(row == -1)==2 and row[3]==row[1]==-1: # phi12 
                TachyonPresent=True
                #print("vect tach: ", comb)
                #print("row is: ", row)
                break
            
            else:
                pass #still tach free
        else: 
            pass
        
        
        #check +z2 tachyons (3,7) - ignore z2  construct new 'row'
        
        z2row=np.zeros(8) # 
        projlstz2=projlst.copy()
        projlstz2.remove(z2) #z2
        for proj in projlstz2: 
            
            #print("projlst.index(proj):", projlst.index(proj))
            #print("for combLst: ", combLst)
            #print("with proj: ", proj)
            phase=[GSO[proj][i] for i in combLst] 
            phase.append(GSO[proj][z2])
            #print(np.prod(phase))
            z2row[projlst.index(proj)]=np.prod(phase)
        
        z2row[0]=-1*z2row[0]
        #print("z2row: ", z2row)
        
        if np.count_nonzero(z2row == -1)==0:
            TachyonPresent=True
            #print("z2 tach: ", comb)
            #print("row is: ", z2row)
            break
        else:
            pass #still tach free
            
        #check +z1 tachyons (3,7) - ignore z1 and Stilde- construct new 'row'
        #combLst.append(11)
        z1row=np.zeros(8) #so 0 and 1st elements should be zero
        projlstz1=projlst.copy()
        projlstz1.remove(z1)
        
        for proj in projlstz1:
            
            #print("projlst.index(proj):", projlst.index(proj))
            #print("for combLst: ", combLst)
            #print("with proj: ", proj)
            phase=[GSO[proj][i] for i in combLst] 
            #print(np.prod(phase))
            phase.append(GSO[proj][z1])
            z1row[projlst.index(proj)]=np.prod(phase)
        
        z1row[0]=-1*z1row[0]
        phasex=[x(i,GSO) for i in combLst]
        phasex.append(x(z1,GSO))
        z1row[-1]=np.prod(phasex) 
        #print("z1row: ", z1row)
        
        if np.count_nonzero(z1row == -1)==0:
            TachyonPresent=True
            #print("z1 tach: ", comb)
            #print("row is: ", z1row)
            break
        else:
            pass #still tach free
    
    return TachyonPresent
            
def tachyonCheckerNew2s(GSO):

    full=[2,3,4,5,6,7]
    TachyonPresent=False
    
    for comb in comb2:
        combLst=list(comb)
        row=np.zeros(10) # S, z1,z2,ek,el,em,en,x,b1/2/3 IN THAT ORDER
        projes=[x for x in full if x not in combLst]
        projlst=[S,z1,z2]+projes
        #print(projlst)
        #print("for comb: ", comb)
        for proj in projlst:
            #print("projlst.index(proj):", projlst.index(proj))
            #print("for combLst: ", combLst)
            #print("with proj: ", proj)
            phase=[GSO[proj][i] for i in combLst] 
            #print(np.prod(phase))
            row[projlst.index(proj)]=np.prod(phase)
        row[0]=-1*row[0]
        
        
        xphase=[x(i,GSO) for i in combLst]
        row[-2]=np.prod(xphase)
        
        if combLst==[6,7]:
            b3phase=[x(i,GSO)*GSO[b1][i]*GSO[b2][i] for i in combLst]
            row[-1]=np.prod(b3phase)
        elif combLst==[4,5]:
            b2phase=[GSO[b2][i] for i in combLst]
            row[-1]=np.prod(b2phase)
        elif combLst==[2,3]:
            b1phase=[GSO[b1][i] for i in combLst]
            row[-1]=np.prod(b1phase)
        else:
            pass
            
        #print("row: ", row)
        
        #check vectorial tachyons:
        if row[0]==1:
            if np.count_nonzero(row == -1)==1 and row[3]==1:
                TachyonPresent=True
                #print("vect tach 1: ", comb)
                #print("row is: ", row)
                break
            
            elif np.count_nonzero(row == -1)==2 and row[3]==row[1]==-1: #phi12
                TachyonPresent=True
                #print("vect tach 2: ", comb)
                #print("row is: ", row)
                break
            
            elif np.count_nonzero(row == -1)==2 and row[3]==row[-2]==-1: # psi^45 oscillator
                TachyonPresent=True
                #print("vect tach 3: ", comb)
                #print("row is: ", row)
                break
            
            elif np.count_nonzero(row == -1)==2 and row[-1]==row[-2]==-1: #psi123/eta1/2/3 in e1+e2,e3+e4,e5+e6 cases
                    
                    TachyonPresent=True
                    #print("vect tach 4: ", comb)
                    #print("row is: ", row)
                    break
            
            elif np.count_nonzero(row == -1)==2 and row[-1]==-1: 
                if row[4]==-1 or row[5]==-1 or row[6]==-1 or row[7]==-1: #i.e. one of y/w's in case of e1+e2,e3+e4,e5+e6
                    TachyonPresent=True
                    #print("vect tach 5: ", comb)
                    #print("row is: ", row)
                    break
            elif np.count_nonzero(row == -1)==3 and row[-1]==row[-2]==row[3]==-1: 
                
                TachyonPresent=True
                #print("vect tach 5: ", comb)
                #print("row is: ", row)
                break
            
            else:
                pass #still tach free
        else:
            pass
        
        
        #check +z2 tachyons (3,7) - ignore z2 
        
        z2row=np.zeros(10)
        projlstz2=projlst.copy()
        projlstz2.remove(z2) #z2
        for proj in projlstz2: #S,z1,+projes 
            
            #print("projlst.index(proj):", projlst.index(proj))
            #print("for combLst: ", combLst)
            #print("with proj: ", proj)
            phasez2=[GSO[proj][i] for i in combLst] 
            phasez2.append(GSO[proj][z2])
            #print(np.prod(phase))
            z2row[projlst.index(proj)]=np.prod(phasez2)
        
        z2row[-2]=x(z2,GSO)*np.prod([x(i,GSO) for i in combLst]) 
        
        if combLst==[6,7]:
            z2row[-1]=np.prod([x(i,GSO)*GSO[b1][i]*GSO[b2][i] for i in combLst])*x(z2,GSO)*GSO[z2][b1]*GSO[z2][b2]
            
        elif combLst==[4,5]:
            z2row[-1]=np.prod([GSO[b2][i]  for i in combLst])*GSO[b2][z2]
        elif combLst==[2,3]:
            z2row[-1]=np.prod([GSO[b1][i]  for i in combLst])*GSO[b1][z2]
        else:
            pass
        #print("z2row: ",z2row)
        if np.count_nonzero(z2row == -1)==0:
            TachyonPresent=True
            #print("z2 tach: ", comb)
            #print("row is: ", z2row)
            break
        else:
            pass #still tach free
        
        #check +z1 tachyons (3,7) - 
        #combLst.append(11)
        z1row=np.zeros(10)
        projlstz1=projlst.copy()
        projlstz1.remove(z1)
        
        for proj in projlstz1:
            
            #print("projlst.index(proj):", projlst.index(proj))
            #print("for combLst: ", combLst)
            #print("with proj: ", proj)
            phase=[GSO[proj][i] for i in combLst] 
            phasez1=np.prod(phase)*GSO[proj][z1]
            z1row[projlst.index(proj)]=np.prod(phasez1)
        
        
        z1row[-2]=x(z1,GSO)*np.prod([x(i,GSO) for i in combLst])
        
        if combLst==[6,7]:
            phase1=[x(i,GSO)*GSO[b1][i]*GSO[b2][i] for i in combLst]
            phase1.append(x(z1,GSO)*GSO[b1][z1]*GSO[b2][z1])
            z1row[-1]=np.prod(phase1)
            
        elif combLst==[4,5]:
            phase2=[GSO[b2][i] for i in combLst]
            phase2.append(GSO[b2][z1])
            z1row[-1]=np.prod(phase2)
            
        elif combLst==[2,3]:
            phase3=[GSO[b1][i] for i in combLst]
            phase3.append(GSO[b1][z1])
            z1row[-1]=np.prod(phase3)
        
        else:
            pass
        
        #print("z1row: ", z1row)
        
        if np.count_nonzero(z1row == -1)==0:
            TachyonPresent=True
            #print("z1 tach: ", comb)
            #print("row is: ", z1row)
            break
        else:
            pass #still tach free
    
    return TachyonPresent    
    
def tachyonCheckerNew1s(GSO): 
    
    
    
    TachyonPresent=False
    e1row=[GSO[e1][S],GSO[e1][z1],GSO[e1][z2],GSO[e1][e2],GSO[e1][e3],GSO[e1][e4],GSO[e1][e5],GSO[e1][e6],x(e1,GSO),GSO[e1][b1]]
    e2row=[GSO[e2][S],GSO[e2][z1],GSO[e2][z2],GSO[e2][e1],GSO[e2][e3],GSO[e2][e4],GSO[e2][e5],GSO[e2][e6],x(e2,GSO),GSO[e2][b1]]
    e3row=[GSO[e3][S],GSO[e3][z1],GSO[e3][z2],GSO[e3][e1],GSO[e3][e2],GSO[e3][e4],GSO[e3][e5],GSO[e3][e6],x(e3,GSO),GSO[e3][b2]]
    e4row=[GSO[e4][S],GSO[e4][z1],GSO[e4][z2],GSO[e4][e1],GSO[e4][e2],GSO[e4][e3],GSO[e4][e5],GSO[e4][e6],x(e4,GSO),GSO[e4][b2]]
    e5row=[GSO[e5][S],GSO[e5][z1],GSO[e5][z2],GSO[e5][e1],GSO[e5][e2],GSO[e5][e3],GSO[e5][e4],GSO[e5][e6],x(e5,GSO),GSO[e5][b1]*GSO[e5][b2]*x(e5,GSO)]
    e6row=[GSO[e6][S],GSO[e6][z1],GSO[e6][z2],GSO[e6][e1],GSO[e6][e2],GSO[e6][e3],GSO[e6][e4],GSO[e6][e5],x(e6,GSO),GSO[e6][b1]*GSO[e6][b2]*x(e6,GSO)]
    
    while TachyonPresent is False:
    
        #check vectorial tachyons:
        if e1row[0]==1:
            if e1row.count(-1)==1 and e1row[-1]==1: # covers some y/w cases , eta23,phis 
                TachyonPresent=True 
                break
                
            elif e1row.count(-1)==2 and e1row[-2]==e1row[-1]==-1: # psi12345,eta1
                TachyonPresent=True
                break
                
            elif e1row.count(-1)==2 and e1row[-1]==-1: # ybar3456 case 
                if e1row[3:8].count(-1)==1:
                    TachyonPresent=True
                    break
                
            else:
                pass #still tach free
        else:
            pass
        
        if e2row[0]==1:
            if e2row.count(-1)==1 and e2row[-1]==1: # covers some y/w cases , eta23,phis 
                TachyonPresent=True 
                break
                
            elif e2row.count(-1)==2 and e2row[-2]==e2row[-1]==-1: # psi12345,eta1
                TachyonPresent=True
                break
                
            elif e2row.count(-1)==2 and e2row[-1]==-1: # ybar3456 case 
                if e2row[3:8].count(-1)==1:
                    TachyonPresent=True
                    break
                
            else:
                pass #still tach free
        else:
            pass
        
        if e3row[0]==1:
            if e3row.count(-1)==1 and e3row[-1]==1: # covers some y/w cases , eta23,phis 
                TachyonPresent=True 
                break
                
            elif e3row.count(-1)==2 and e3row[-2]==e3row[-1]==-1: # psi12345,eta1
                TachyonPresent=True
                break
                
            elif e3row.count(-1)==2 and e3row[-1]==-1: # ybar3456 case 
                if e3row[3:8].count(-1)==1:
                    TachyonPresent=True
                    break
                
            else:
                pass #still tach free
        else:
            pass
        
        if e4row[0]==1:
            if e4row.count(-1)==1 and e4row[-1]==1: # covers some y/w cases , eta23,phis 
                TachyonPresent=True 
                break
                
            elif e4row.count(-1)==2 and e4row[-2]==e4row[-1]==-1: # psi12345,eta1
                TachyonPresent=True
                break
                
            elif e4row.count(-1)==2 and e4row[-1]==-1: # ybar3456 case 
                if e4row[3:8].count(-1)==1:
                    TachyonPresent=True
                    break
                
            else:
                pass #still tach free
        else:
            pass
        #print("e5row: ", e5row)
        if e5row[0]==1:
            if e5row.count(-1)==1 and e5row[-1]==1: # covers some y/w cases , eta23,phis 
                TachyonPresent=True 
                break
                
            elif e5row.count(-1)==2 and e5row[-2]==e5row[-1]==-1: # psi12345,eta1
                TachyonPresent=True
                break
                
            elif e5row.count(-1)==2 and e5row[-1]==-1: # ybar3456 case 
                if e5row[3:8].count(-1)==1:
                    TachyonPresent=True
                    break
                
            else:
                pass #still tach free
        else:
            pass
            
        if e6row[0]==1:
            if e6row.count(-1)==1 and e6row[-1]==1: # covers some y/w cases , eta23,phis 
                TachyonPresent=True 
                break
                
            elif e6row.count(-1)==2 and e6row[-2]==e6row[-1]==-1: # psi12345,eta1
                TachyonPresent=True
                break
                
            elif e6row.count(-1)==2 and e6row[-1]==-1: # ybar3456 case 
                if e6row[3:8].count(-1)==1:
                    TachyonPresent=True
                    break
                
            else:
                pass #still tach free
        else:
            pass
            
        #ei+z1
        
        e1z1row=np.array([-1*GSO[S][e1]*GSO[S][z1],GSO[e2][e1]*GSO[e2][z1],GSO[e3][e1]*GSO[e3][z1],GSO[e4][e1]*GSO[e4][z1],GSO[e5][e1]*GSO[e5][z1],GSO[e6][e1]*GSO[e6][z1],GSO[z2][e1]*GSO[z2][z1],x(e1,GSO)*x(z1,GSO),GSO[b1][e1]*GSO[b1][z1]])
        e2z1row=np.array([-1*GSO[S][e2]*GSO[S][z1],GSO[e1][e2]*GSO[e1][z1],GSO[e3][e2]*GSO[e3][z1],GSO[e4][e2]*GSO[e4][z1],GSO[e5][e2]*GSO[e5][z1],GSO[e6][e2]*GSO[e6][z1],GSO[z2][e2]*GSO[z2][z1],x(e2,GSO)*x(z1,GSO),GSO[b1][e2]*GSO[b1][z1]])
        e3z1row=np.array([-1*GSO[S][e3]*GSO[S][z1],GSO[e1][e3]*GSO[e1][z1],GSO[e2][e3]*GSO[e2][z1],GSO[e4][e3]*GSO[e4][z1],GSO[e5][e3]*GSO[e5][z1],GSO[e6][e3]*GSO[e6][z1],GSO[z2][e3]*GSO[z2][z1],x(e3,GSO)*x(z1,GSO),GSO[b2][e3]*GSO[b2][z1]])
        e4z1row=np.array([-1*GSO[S][e4]*GSO[S][z1],GSO[e1][e4]*GSO[e1][z1],GSO[e2][e4]*GSO[e2][z1],GSO[e3][e4]*GSO[e3][z1],GSO[e5][e4]*GSO[e5][z1],GSO[e6][e4]*GSO[e6][z1],GSO[z2][e4]*GSO[z2][z1],x(e4,GSO)*x(z1,GSO),GSO[b2][e4]*GSO[b2][z1]])
        e5z1row=np.array([-1*GSO[S][e5]*GSO[S][z1],GSO[e1][e5]*GSO[e1][z1],GSO[e2][e5]*GSO[e2][z1],GSO[e3][e5]*GSO[e3][z1],GSO[e4][e5]*GSO[e4][z1],GSO[e6][e5]*GSO[e6][z1],GSO[z2][e5]*GSO[z2][z1],x(e5,GSO)*x(z1,GSO),GSO[e5][b1]*GSO[e5][b2]*x(e5,GSO)*GSO[z1][b1]*GSO[z1][b2]*x(z1,GSO)])
        e6z1row=np.array([-1*GSO[S][e6]*GSO[S][z1],GSO[e1][e6]*GSO[e1][z1],GSO[e2][e6]*GSO[e2][z1],GSO[e3][e6]*GSO[e3][z1],GSO[e4][e6]*GSO[e4][z1],GSO[e5][e6]*GSO[e5][z1],GSO[z2][e6]*GSO[z2][z1],x(e6,GSO)*x(z1,GSO),GSO[e6][b1]*GSO[e6][b2]*x(e6,GSO)*GSO[z1][b1]*GSO[z1][b2]*x(z1,GSO)])
        eiz1=np.vstack((e1z1row,e2z1row,e3z1row,e4z1row,e5z1row,e6z1row))
        
        #ei+z2
        e1z2row=np.array([-1*GSO[S][e1]*GSO[S][z2],GSO[z1][e1]*GSO[z1][z2],GSO[e2][e1]*GSO[e2][z2],GSO[e3][e1]*GSO[e3][z2],GSO[e4][e1]*GSO[e4][z2],GSO[e5][e1]*GSO[e5][z2],GSO[e6][e1]*GSO[e6][z2],x(e1,GSO)*x(z2,GSO),GSO[b1][e1]*GSO[b1][z2]])
        e2z2row=np.array([-1*GSO[S][e2]*GSO[S][z2],GSO[z1][e2]*GSO[z1][z2],GSO[e1][e2]*GSO[e1][z2],GSO[e3][e2]*GSO[e3][z2],GSO[e4][e2]*GSO[e4][z2],GSO[e5][e2]*GSO[e5][z2],GSO[e6][e2]*GSO[e6][z2],x(e2,GSO)*x(z2,GSO),GSO[b1][e2]*GSO[b1][z2]])
        e3z2row=np.array([-1*GSO[S][e3]*GSO[S][z2],GSO[z1][e3]*GSO[z1][z2],GSO[e1][e3]*GSO[e1][z2],GSO[e2][e3]*GSO[e2][z2],GSO[e4][e3]*GSO[e4][z2],GSO[e5][e3]*GSO[e5][z2],GSO[e6][e3]*GSO[e6][z2],x(e3,GSO)*x(z2,GSO),GSO[b2][e3]*GSO[b2][z2]])
        e4z2row=np.array([-1*GSO[S][e4]*GSO[S][z2],GSO[z1][e4]*GSO[z1][z2],GSO[e1][e4]*GSO[e1][z2],GSO[e2][e4]*GSO[e2][z2],GSO[e3][e4]*GSO[e3][z2],GSO[e5][e4]*GSO[e5][z2],GSO[e6][e4]*GSO[e6][z2],x(e4,GSO)*x(z2,GSO),GSO[b2][e4]*GSO[b2][z2]])
        e5z2row=np.array([-1*GSO[S][e5]*GSO[S][z2],GSO[z1][e5]*GSO[z1][z2],GSO[e1][e5]*GSO[e1][z2],GSO[e2][e5]*GSO[e2][z2],GSO[e3][e5]*GSO[e3][z2],GSO[e4][e5]*GSO[e4][z2],GSO[e6][e5]*GSO[e6][z2],x(e5,GSO)*x(z2,GSO),GSO[e5][b1]*GSO[e5][b2]*x(e5,GSO)*GSO[z2][b1]*GSO[z2][b2]*x(z2,GSO)])
        e6z2row=np.array([-1*GSO[S][e6]*GSO[S][z2],GSO[z1][e6]*GSO[z1][z2],GSO[e1][e6]*GSO[e1][z2],GSO[e2][e6]*GSO[e2][z2],GSO[e3][e6]*GSO[e3][z2],GSO[e4][e6]*GSO[e4][z2],GSO[e5][e6]*GSO[e5][z2],x(e6,GSO)*x(z2,GSO),GSO[e6][b1]*GSO[e6][b2]*x(e6,GSO)*GSO[z2][b1]*GSO[z2][b2]*x(z2,GSO)])
        eiz2=np.vstack((e1z2row,e2z2row,e3z2row,e4z2row,e5z2row,e6z2row))
        
        def spinPhaseCheck(row):
            
            TachPresent=False
            if np.count_nonzero(row == -1)==0:
                TachPresent=True
                
            else:
                pass #still tach free
                
            return TachPresent
        
        spinTachs1=np.apply_along_axis( spinPhaseCheck, axis=1, arr=eiz1 )
        spinTachs2=np.apply_along_axis( spinPhaseCheck, axis=1, arr=eiz2 )
        #print(spinTachs1)
    
        for bool in spinTachs1:
            if bool==True:
                TachyonPresent=True
                #print("ei z1 tach")
                #print(spinTachs1)
                break
            else: 
                pass
            
        for bool in spinTachs2:
            if bool==True:
                TachyonPresent=True
                #print("ei z2 tach")
                #print(spinTachs2)
                break
            else: 
                pass
        
        z1tachs=np.array([GSO[z1][S],GSO[z1][e1],GSO[z1][e2],GSO[z1][e3],GSO[z1][e4],GSO[z1][e5],GSO[z1][e6],GSO[z1][b1],GSO[z1][b2],GSO[z1][z2]])
        z2tachs=np.array([GSO[z2][S],GSO[z2][e1],GSO[z2][e2],GSO[z2][e3],GSO[z2][e4],GSO[z2][e5],GSO[z2][e6],GSO[z2][b1],GSO[z2][b2],GSO[z2][z1]])
        #check z1 and z2 tachs
        if np.count_nonzero(z1tachs == -1)==0:
            TachyonPresent=True
            #print("tach from z1")
            break
        if np.count_nonzero(z2tachs == -1)==0:
            TachyonPresent=True
            #print("tach from z2")
        break
            
    
    return TachyonPresent

import random

def GSOmatrix():
    G=np.zeros((12,12))
    G[0][0] = -1 #still applies?

    #randomise upper triangle
    for i in range(0,12):
        for j in range(1,12):
            if (j>i):
                G[i][j]=random.choice([-1,1])
                G[j][i]=np.real(np.around(np.exp(1j*np.pi*BP[i][j]/2)))*G[i][j]
            else:
                pass

    #fix diagonal
    for i in range(0,12):
        G[i][i] = -np.real(np.around(np.exp(1j*np.pi*BP[i][i]/4))*G[i][0])

    #avoid susy case...

    if G[S][e1]==-1 and G[S][e2]==-1 and G[S][e3]==-1 and G[S][e4]==-1 and G[S][e5]==-1 and G[S][e6]==-1 and G[S][z1]==-1 and G[S][z2]==-1:
        G[S][e1]=1
        
    
    return G

from pandas import DataFrame
import csv

import multiprocessing as mp

import timeit

def RandomClassfn(N):

    
    GSOt=GSOmatrix()
    
    #avoid susy case:
    
    if tachyonCheckerNew3s(GSOt) is False:
    
        if tachyonCheckerNew2s(GSOt) is False:
                #print("No 2s")
            if tachyonCheckerNew1s(GSOt) is False:
            
                return GSOt
                            
            else: return 2
            
        else: return 2

    else: return 2


Workers = mp.cpu_count()

SampleSize=np.zeros(10000)

if __name__ == '__main__':

    print("Number of Parallel Processes: ", Workers)

    pool = mp.Pool(processes=Workers) #Initialise Multiprocessing Pool of Workers

    start = timeit.default_timer()
    TRes = pool.map(RandomClassfn,SampleSize)
    stop = timeit.default_timer()
    print("Time:", stop - start)

    #print("Number tachyonic:", TRes.count(2))
    #print("Number tachyon free but with chiral exotics: ", TRes.count(3))
    #print("Numb tach free, no chiral exots but not 3 gen: ", TRes.count(4))
    #print("Numb tach free, no chiral exots, 3 gen but no bidoub: ", TRes.count(5))
    #print("TRes = ",TRes)
    numModels=0
    for item in TRes:

        if type(item)==np.ndarray:
            numModels+=1
            #print("item is: ", item)
            fd = open("TachFreeSO10Batch.csv", "a")

            GSO_dframe=DataFrame(item)
            #print("DataFrame is: ", GSO_dframe)
            GSO_dframe.to_csv(fd,header=None,index=False)
            fd.close()

    print("Num tachfree models:", numModels)
    
    
    
    