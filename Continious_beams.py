# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 15:45:12 2023

@author: user
"""
import numpy as np


user_input = input("Enter length of elements in meter separated by spaces: ")
l = user_input.split()
l = [int(element) for element in l]
print("Length of elements are:", l)
N_element= len(l)
N_node = N_element+1
E = float(input("Enter the Elastic Modulus in Pascal :"))
I = float(input("Enter the Moment of Inertia m^4 :"))
EI = E*I
K=np.zeros([2*N_node,2*N_node])
for i in range(N_element):
    n = i
    f = n+1
    K[2*n,2*n]=K[2*n,2*n]+12/(l[i]**3)
    K[2*n,2*n+1]=K[2*n,2*n+1]+6/(l[i]**2)
    K[2*n,2*f]=K[2*n,2*f]-12/(l[i]**3)
    K[2*n,2*f+1]=K[2*n,2*f+1]+6/(l[i]**2)
    
    K[2*n+1,2*n]=K[2*n+1,2*n]+6/(l[i]**2)
    K[2*n+1,2*n+1]=K[2*n+1,2*n+1]+4/l[i]
    K[2*n+1,2*f]=K[2*n+1,2*f]-6/(l[i]**2)
    K[2*n+1,2*f+1]=K[2*n+1,2*f+1]+2/l[i]
    
    K[2*f,2*n]=K[2*f,2*n]-12/(l[i]**3)
    K[2*f,2*n+1]=K[2*f,2*n+1]-6/(l[i]**2)
    K[2*f,2*f]=K[2*f,2*f]+12/(l[i]**3)
    K[2*f,2*f+1]=K[2*f,2*f+1]-6/(l[i]**2)
    
    K[2*f+1,2*n]=K[2*f+1,2*n]+6/(l[i]**2)
    K[2*f+1,2*n+1]=K[2*f+1,2*n+1]+4/l[i]
    K[2*f+1,2*f]=K[2*f+1,2*f]-6/(l[i]**2)
    K[2*f+1,2*f+1]=K[2*f+1,2*f+1]+4/l[i]

print('Global stiffness matrix')    
print(K*E*I)

#******Loading ******#
F=np.zeros(2*N_node)

for j in range(N_element):
   
       t = int(input("Enter the loading pattern in element wise type 1 for UDL:"))
       
       if (t == 1):
           u= int(input("Enter the UDL value in N/m :"))
            
           F[2*j]= F[2*j] + u*l[j]/2
           F[2*j+1]= F[2*j+1] + u*l[j]**2/12
           F[2*j+2]= F[2*j+2] + u*l[j]/2
           F[2*j+3]= F[2*j+3] - u*l[j]**2/12
        
           
       else :
           p= int(input("Enter the Concentrated load value acting at midpoint in N:"))
          
           F[2*j]= F[2*j] + p/2
           F[2*j+1]= F[2*j+1] + p*l[j]/8
           F[2*j+2]= F[2*j+2] + p/2
           F[2*j+3]= F[2*j+3] - p*l[j]/8

print('Global load vector')
print(F)

#******Boundary Conditions ********#
user_input = input("Enter Degree of freedom where Displacement is zero sequentially : ")
d = user_input.split()
d = [int(element) for element in d]
N_d= len(d)

for i in range(N_d):
    s= d[i]-i
    
    K=np.delete(K, s, 0)  
    K=np.delete(K, s, 1)
    F=np.delete(F,s,0)


K_modified= K*EI
print('Stiffness Matrix after elimination')
print (K_modified)

#******Solutions******#
print('Global load vector after elimination')
print(F)
All_DOF= np.zeros(2*N_node)
Q= np.linalg.solve(K_modified, F)
Qg=np.linspace(0,2*N_node-1,2*N_node) # all the dof indices
All_DOF=np.linspace(0,0,2*N_node)# displacements at all the dofs
NZQ=np.setdiff1d(Qg,d)  # the array containing the indices of the free dofs
ZQ=np.setdiff1d(Qg,NZQ)  # the array containing the indices of the constrained dofs

nNZDOF= 2*N_node - N_d

for i in range(nNZDOF):
  j=int(NZQ[i])
  All_DOF[j]=Q[i]
print('Displacement at all dofs')
print(All_DOF)


#******End of code*******#




           
        
         
       
    
    
    
    
    
    
    
    