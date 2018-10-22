# -*- coding: utf-8 -*-
"""
Monte Carlo Code for calculating bare infinite slab characteristics

Authored by: Daniel Nichols 10/9/2018
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from datetime import datetime

'''
To run this program just modify the constant variables. Global variables
should be left alone. Note that the thickness specified is the half thickness.
'''

#%% Constant Variables
Sig_c = 0.011437;Sig_s = 0.05;Sig_f = 0.013         # [cm^-1] cross sections
Sig_t = Sig_c + Sig_s + Sig_f                       # [cm^-1] total xs
nu = 2.5                                            # Avg num of n produced
N = 1000000                                         # Num of particles
T = 78.000                                          # Initial thickness
n_b = 10                                            # Num of batches
N_b = N/n_b                                         # Num of histories/batch
bins =100                                           # Num of bins for histogram


#%% Global Variables
FDP = np.zeros(((bins+1),2))                        # Fission Density Profile
GENS = np.zeros(((bins+1),n_b+1))                   # FDP for n_b generations
GENS[:,0] = [x for x in np.linspace(0,2*T,bins+1)]  # Add x coordinates
k_a=k_f=k_c=k_t = 0.0                               # Initialize  k estimators
e_f=e_a=e_c=e_t = 0.0                               # Initialize error
# Initialize variables
S_f = S_c = S_a = S_t = 0.0                         # Init all scores to zero
S2_f = S2_c = S2_a = S2_t = 0.0                     # Init all std
N_L = N_f = N_c = N_s = 0                           # Init all counters
D = 0.0                                             # Init track length

#%% Time marker
def current_time():
    return datetime.now().strftime('%I:%M:%S %p')

#%%  Park and Miller Function for random number generation
def rand( seed = 1, A = 16807. , M = 2147483647.):
    minv = 1. / M
    yield minv * seed
    while True :
        seed = (A * seed ) % M
        yield minv * seed
#%% Functions for neutron tracking
def tracking(m,gen,CDF):
  ''' Overarching function to track single neutron'''
  x=0;delt=0
  a = start(m,gen,CDF)*2*T            # Record birth location
  x = a                               # Initialize particle start location
  D_i = 0                             # Tracking total displacement
  c = 0                               # Tracks number of scatters
  r=direction(m)                      # start wih no trajectory
  k=0                                 # no starting magnitude
  global D; global N_f; global N_c; global N_s; global N_L
  global S_f; global S_c; global S_a; global S_t
  global S2_f; global S2_c; global S2_a; global S2_t
  while x <= 2*T and x >= 0:          # while neutron is in slab
    r, k, i  = direction(m), displace(m), interaction(m)
    x += r *k                         # displace neutron location
    if x >=2*T or x <= 0:             # confirm while loop is during its job
      break
    else:
      pass
    D_i += abs(k)                     # total displacement
    if i == "fission":
      D +=D_i; N_f +=1                # Add to distance, Fission tally
      N_s+=c
      #S_f += nu; S2_f += nu**2        # increase average tally and std tally
      return a, x, i                  # birth, death, case
    elif i == "capture":
      #if c == 0:
        #S_c += nu*(Sig_f/Sig_t)
        #S2_c +=(nu*(Sig_f/Sig_t))**2
      #else:
        #S_c += nu*(Sig_f/Sig_t)*float(c+1)
        #S2_c +=(nu*(Sig_f/Sig_t)*float(c+1))**2
      N_s+=c
      D +=D_i; N_c +=1                # Add to distance, Capture tally
      return a, x, i                  # birth, death, case
    else:
      c += 1
  #print x
  i='leakage'
  q = x - r*k                         # Subtract last displacemnet
  if r < 0:                           # determine distance to nearest slab face
    delt = -q/r
  else:
    delt = (2*T-q)/r
  D +=(D_i+delt); N_s +=c; N_L +=1    # Add to distance, Scat and leak tally
  return a, x, i                     # birth, death, case
#%%
def start(r,gen,CDF):
  ''' Initialize particle birth location'''
  x=next(r)
  if gen == 1:
    return x
  else:
    return np.interp(x, CDF, FDP[:,0])/(2*float(T))
#%%
def direction(r):
  ''' function to determine direction of travel for particle'''
  x = next(r)
  if x <= 0.5:                       # Determine positive or negative direction
    return -1
  else:
    return 1
#%%
def displace(r):
  ''' Function to determine path length of neutron'''     
  d = -(1/Sig_t)*np.log(next(r))    # Sample path length
  return d
#%%
def interaction(r):
  ''' Function to determine interaction type'''                      
  x = next(r)
  if x <= Sig_f/Sig_t:              # Determine interaction type
    return "fission"
  elif x <= (Sig_f+Sig_c)/Sig_t and  x > Sig_f/Sig_t:
    return "capture"
  else:
    return "scatter"
#%%
def distr(FDP):
  ''' function used to find location of neutron birth using distribution'''
  CDF = [1/np.sum(FDP[:,1])*np.sum(FDP[0:i,1]) for i in range(len(FDP[:,1]))]
  return CDF
#%% 
def batch(T,bins,gen):
  global r; global N_b; global FDP
  case, fiss, birth = [],[],[]        # Int type, fiss loc, birth loc
  if gen != 1:
    CDF = distr(FDP)
  else:
    CDF = 0
  for i in range(N_b):
    [f,g,h] = tracking(r,gen,CDF)
    birth.append(f); case.append(h)   # Track birth and final interaction type
    if h == "fission":
      fiss.append(g)                  # Track fission locations
    if i % (N_b/5) == 0:
      print "Gen {}, {}% complete, time: {}"\
      .format(str(gen),str(100*i/N_b),str(current_time()))

  #%% build source from histogram
  '''Scans fission location values and builds FDPtribution list '''

  [x,y] = sp.histogram(fiss,bins,(0,2*T))
  x = np.append(x,(x[0]))                                # Apply sym conditions
  #x = np.append(x,(x[-1]+(x[-1]-x[-2])))                # Apply sym conditions
  [FDP[:,0],FDP[:,1]] = y,x

  #%% K_eff methods
  global k_f; global k_c; global k_a; global k_t
  N_col = N_c + N_s + N_f
  k_f = nu*N_f/float(N_b)                              # Fission estimator
  k_c = ((N_col)*((nu*Sig_f)/Sig_t))/float(N_b)  # Collison Estimator
  k_a = (N_c+N_f)/float(N_b)*(nu*Sig_f/(Sig_c+Sig_f))  # Absorption Estimator
  k_t = D*Sig_f*nu/float(N_b)                          # Track Length Estimator
  
  global e_f; global e_c; global e_a; global e_t
  e_f = np.sqrt((1/(float(N_b-1)))*abs(((float(N_f)*nu**2)/float(N_b))-k_f**2))
  
  e_c = np.sqrt((1/(float(N_b-1)))\
                *abs(((float(N_col)*(nu*Sig_f/Sig_t)**2)/float(N_b))-k_c**2))
  
  e_a = np.sqrt((1/(float(N_b-1)))\
                *abs(((float(N_f+N_c)\
                    *(nu*Sig_f/(Sig_c+Sig_f))**2)/float(N_b))-k_a**2))
  
  e_t = np.sqrt((1/(float(N_b-1)))\
                *abs(((float(D)*(Sig_f*nu)**2)/float(N_b))-k_t**2))
  
  return fiss,birth
#%%####################### PROGRAM START ######################################
r = rand()                                             # Begin rand num gen
events = np.zeros((14,n_b+1),dtype='a9')
names = ['fiss frac','col frac','scatters','leak frac','k_f','error_f',\
         'k_c','error_c','k_a','error_a','k_track','error_track','k_avg','error_avg']
events[:,0] = [i for i in names]
#global D; global N_s; global N_L
for gen in range(1,1+n_b,1):
  [fiss,birth] = batch(T,bins,gen)
  GENS[:,gen] = [FDP[u,1] for u in range(bins+1)]
  events[:,gen] = [float(N_f)/float(N_b),float(N_c)/float(N_b)\
        ,N_s,float(N_L)/float(N_b),k_f,e_f,k_c,e_c,k_a,e_a,k_t,e_t,\
        (k_c+k_f+k_a+k_t)/float(4),\
        np.sqrt(e_c**2+e_f**2+e_a**2+e_t**2)/float(4)]
  D = 0.0; N_s = N_L =N_f=N_c =0
  # Print out all estimators
  print"\nGen {}: \nK_f: {}\tK_c: {}\tK_a: {}\tK_t: {}"\
  .format(str(gen),str(round(k_f,7)),str(round(k_c,7)),\
          str(round(k_a,7)),str(round(k_t,7)))
  print"   +-{}\t   +-{}\t   +-{}\t   +-{}\n"\
  .format(str(round(e_f,7)),str(round(e_c,7)),\
          str(round(e_a,7)),str(round(e_t,7)))                                           

###############################################################################  
#%% Plotting Section
plt.figure(1)
plt.figure(figsize=(11,8.5))
plt.title('Fission Density Profile')
for i in range(1,1+n_b,1):
  plt.plot(GENS[:,0],GENS[:,i],'-', label='Gen {}'.format(str(i)))
plt.xlim(0,2*T), plt.xlabel('Location (cm)'),plt.ylabel('# of fissions')
plt.legend(),plt.grid('True'),plt.savefig("FDP.jpg"), plt.show()


del D,N,N_b,N_c,N_f,N_s,S2_a,S2_c,S2_f,S2_t,S_a,S_c,S_f,S_t,Sig_c,Sig_f,Sig_s
del Sig_t,T,fiss,i,k_a,k_c,k_f,k_t,n_b,nu,bins,birth,gen,u,x,e_c,e_a,e_f,e_t
del names