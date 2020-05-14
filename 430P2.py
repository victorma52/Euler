# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 13:56:45 2019

@author: Victor
"""

import numpy as np
from matplotlib import pyplot as plt
#forward euler step positive trial
t0=0;
#t0=28.5;
n0=1;
c0=(0.0075)/(0.08*6E-5);
tf=30;
dt=0.5;
length=int((tf-t0)/(dt));
a1=(0.0015-0.0075)/(6E-5)
a2=(-0.050025-0.0075)/(6E-5)
b=(0.0075/6E-5)
lam=0.08
#A inputs for functions are Apos or Aneg
Apos=np.array([[a1,lam],[b,-lam]])
Aneg=np.array([[a2,lam],[b,-lam]])
t=np.linspace(t0,tf,length);
f=np.zeros((2,length))
I=np.identity(2)
f[:,0]=np.array([n0,c0])

def FEuler(A,dt):
    length=int((tf-t0)/(dt));
    t=np.linspace(t0,tf,length);
    f=np.zeros((2,length))
    I=np.identity(2)
    f[:,0]=np.array([n0,c0])
    for i in range(1,length):
        Adt=A*dt
        mat=(I+Adt)
        f[:,i]=np.matmul(mat,f[:,i-1])
    N=f[0]
    C=f[1]
    
    plt.plot(t,N)
    plt.title("N over time(t) for Forward Euler")
    plt.show()
    plt.plot(t,C)
    plt.title("C over time(t) for Forward Euler")
    plt.show()

def BEuler(A,dt):
    length=int((tf-t0)/(dt));
    t=np.linspace(t0,tf,length);
    f=np.zeros((2,length))
    I=np.identity(2)
    f[:,0]=np.array([n0,c0])
    for i in range(1,length):
        Adt=A*dt
        mat=np.linalg.inv((I-Adt))
        f[:,i]=np.matmul(mat,f[:,i-1])
        
    N=f[0]
    C=f[1]
    
    plt.plot(t,N)
    plt.title("N over time(t) for Backwards Euler")
    plt.show()
    plt.plot(t,C)
    plt.title("C over time(t) for Backwards Euler")
    plt.show()

def CNick(A,dt):
    length=int((tf-t0)/(dt));
    t=np.linspace(t0,tf,length);
    f=np.zeros((2,length))
    I=np.identity(2)
    f[:,0]=np.array([n0,c0])
    for i in range(1,length):
        Adt=(A*dt)/2
        mat1=np.linalg.inv((I-Adt))
        mat2=(I+Adt)
        X=np.matmul(mat1,mat2)
        f[:,i]=np.matmul(X,f[:,i-1])
        
    N=f[0]
    C=f[1]
    
    plt.plot(t,N,label="N over time(t) for Crank-Nicholson")
    plt.title("N over time(t) for Crank-Nicholson")
    plt.show()
    plt.plot(t,C)
    plt.title("C over time(t) for Crank-Nicholson")
    plt.show()

def BDF2(A,dt):
    length=int((tf-t0)/(dt));
    t=np.linspace(t0,tf,length);
    f=np.zeros((2,length))
    I=np.identity(2)
    f[:,0]=np.array([n0,c0])
    for i in range(1,length):
        Adt=(A*dt)/2
        mat1=np.linalg.inv((I-(7/12)*Adt))
        mat2=(I+(5/12)*Adt)
        X=np.matmul(mat1,mat2)
        f[:,i]=np.matmul(X,f[:,i-1])
        
    N=f[0]
    C=f[1]
    #print(f)
    plt.plot(t,N)
    plt.title("N over time(t) for BDF2")
    plt.show()
    plt.plot(t,C)
    plt.title("C over time(t) for BDF2")
    plt.show()