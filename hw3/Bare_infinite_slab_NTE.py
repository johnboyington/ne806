# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 14:14:49 2018

@author: Daniel Nichols
"""
import mpmath
import numpy as np

T = 26.868*2                      # domain size
N = 10                          # number of nodes

def E(n,z):
  return mpmath.expint(n,z)

def y(x):
  return (T/2)*x+(T/2)


#x = np.linspace(0,T,N)
B = np.zeros((N,N))
X = np.array([-0.973906528517172,-0.865063366688985,-0.679409568299024,\
               -0.433395394129247,-0.148874338981631,0.148874338981631,\
               0.433395394129247,0.679409568299024,0.865063366688985,\
               0.973906528517172])
W = np.array([0.066671344308688,0.149451349150581,0.219086362515982,\
               0.269266719309996,0.295524224714753,0.295524224714753,\
               0.269266719309996,0.219086362515982,0.149451349150581,\
               0.066671344308688])
S=0

for i in range(N):
  for j in range(N):
    if i != j:
      #print int(abs(x[j]-x[i]))
      #print "i: {}, j: {}".format(str(i),str(j))
      #print W[j]*E(1,int(abs(y(X[j])-y(X[i]))))
      #print abs(y(X[j])-y(X[i]))
      #print "x_j: {}, x_i: {}".format(str(y(X[j])),str(y(X[i])))
      B[i,j] = W[j]*E(1,abs(y(X[j])-y(X[i])))
      S += W[j]*E(1,abs(y(X[j])-y(X[i])))
      print S 
    else:
      pass
    B[i,i] = float(2) - E(2,y(X[i]))- E(2,T-y(X[i]))-float(T/(2*S))
    S = 0
    
#for i in range(N):
#  for j in range(N):
#    if i == j:
#      B[i,j] = float(2) - E(2,y(X[i]))- E(2,T-y(X[i]))-float(S)
     
      
P = np.array([1]*10)    # Initial P vector of size (m,1) as starting guess
diff = 10
n = .0001          # This will be the tolerance at which the while loop quits
lam_1 = lam_2 = 0    # Initializes the first lambda as the trivial solution
while diff > n:
  lam_1 = lam_2                               # set lambda_(i-1) = lambda_(i) 
  x = P                                       # set f_(i-1) = f_(i) 
  P = np.matmul(B,P)                          # B*P = new f_(i)
  lam_2 = np.linalg.norm(P)/np.linalg.norm(x) # Find new lambda_(i)
  diff = lam_2-lam_1                          # diff(old l, new l)
  c= 2 /lam_2                                 # Calculate c_crit with max lamda
  #print c
