# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 18:13:13 2019

@author: batuhan
"""
#Batuhan Özkan
# This program to evaluate and draw the axial force vs. moment
# interaction diagram for a reinforced concrete section with multi layer reinforcements.
# INTERACTION CURVE 

import numpy as np
import matplotlib.pyplot as plt
fc=300 # in kgf/cm2
fcd=(fc/1.5)*(1e-3)  # in t/cm2
k1=0.82
Ey=2e6 # kgf/cm2 
fy=4200 # kgf/cm2

eps_cup=0.003
eps_sy=((fy/1.15)/Ey)
dia=2.4 # 24fi (in 2cm is used)
area=((np.pi)*(dia**2))/4    # 1 reinforcement area 

h=35 # cm
bw=50#cm
clear_cover=4
St=np.array([[area*5,clear_cover,Ey,fy],
            [area*2,12,Ey,fy],
            [area*2,23,Ey,fy],
            [area*5,h-clear_cover,Ey,fy]])

#Calculation of cb 
cb=((St[3,1])/((eps_sy/eps_cup)+1))
print("cb(balance depth)= " ,cb)

# CREATING OF INTERACTION DIAGRAM :
def function_Interaction(cb,St,fc,bw,h,eps_cup,eps_sy,fy,k1):
    # First Step : PURE TENSION
    fyd=(fy/1.15)*(1e-3)
    # in t/cm2
    fcd=(fc/1.5)*(1e-3)

    counter=0
    for cp in np.arange(clear_cover,St[3,1]+0.1,0.1):
        counter=counter+1
    
    N_M=np.zeros((counter+2,2),dtype=float)   # FIRST COLUMN IS AXIAL FORCES , AND THE SECOND COLUMN THEIR MOMENTS
    
    # d. Maximum Tensile load
    N_M[0,0]=(-1)*10*((St[0,0]+St[1,0]+St[2,0]+St[3,0])*fyd)
    N_M[0,1]=0
    print("d. Maximum Tensile Load Capacity =",N_M[0,0],"kN")
        
    # a. Maximum Axial Load    
    N_M[counter+1,0]=10*(0.85*fcd*bw*h+((St[0,0]+St[1,0]+St[2,0]+St[3,0])*fyd))
    N_M[counter+1,1]=0
    print("a. Maximum Axial Load Capacity =",N_M[counter+1,0],"kN")
 
    counter_c=0
    for cp in np.arange(clear_cover,St[3,1]+0.1,0.1):
        
        sigma4p=fyd
        Fs4p=St[3,0]*sigma4p
      
        Fcp=0.85*k1*cp*bw*fcd
        
        eps1p=((eps_cup*(cp-St[0,1]))/cp)
        if eps1p>eps_sy :
            sigma1p=fyd
            Fs1p=St[0,0]*sigma1p
        else:
            sigma1p=(St[0,2]*eps1p)*(1e-3)
            Fs1p=St[0,0]*sigma1p
    
       
        eps2p=((eps_cup*(cp-St[1,1]))/cp)
        if eps2p<eps_sy:
            
            sigma2p=(St[1,2]*eps2p)*(1e-3)
            Fs2p=St[1,0]*sigma2p
        else:
            sigma2p=fyd
            Fs2p=St[1,0]*sigma2p
            
        
        eps3p=((eps_cup*(St[2,1]-cp))/cp)
        if eps3p>eps_sy :
            sigma3p=fyd
            Fs3p=St[2,0]*sigma3p
        else:
            sigma3p=(St[0,2]*eps3p)*(1e-3)
            Fs3p=St[2,0]*sigma1p
        
             
        N=Fcp+Fs1p+Fs2p-Fs3p-Fs4p
        Mp=Fs1p*((h/2)-St[0,1])+Fcp*((h/2)-((k1*cp)/2))+Fs2p*((h/2)-St[1,1]) + Fs3p*(St[2,1]-(h/2)) + Fs4p*(St[3,1]-(h/2))
        
        if N<0.01 and N>0:
            print("c. No axial force bending only = ",1e-1*Mp,"kN*m","cp :",cp,"cm")
        
        N_M[counter_c+1,0]=10*N      
        N_M[counter_c+1,1]=Mp
        
        counter_c=counter_c+1
    
    Mb=np.max(N_M)
    a=np.where(N_M==Mb)
       
    row_max=a[0][0]
    M_balanced=N_M[row_max,1]
    N_balanced=N_M[row_max,0]
    print("b. Balanced Axial force is =",N_balanced,"kN",",","Balanced Moment is =",1e-1*M_balanced,"kN.m")
    
    N_axial=(N_M[:,0])   # in kN
    M_bending=1e-1*(N_M[:,1])    # in kN*m
    
    print(N_M)
       
    eq=np.loadtxt("data2.txt",dtype="float")
    #print(eq[0][1])
    #print(N_M)
    
    plt.plot(M_bending,N_axial)
    for index_x in range(0,len(eq),1):
            x_value = eq[index_x,1]
            y_value = eq[index_x,0]
            plt.scatter(x_value, y_value, s=10)
    plt.plot(M_bending,N_axial,'b')
    plt.ylabel("Axial Force(kN)")
    plt.xlabel("Bending Moment(kN*m)")
    plt.title("Interaction Curve")
    
    plt.show()
    
   
function_Interaction(cb,St,fc,bw,h,eps_cup,eps_sy,fy,k1)        