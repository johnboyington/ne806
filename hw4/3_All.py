import numpy as np
import matplotlib.pyplot as plt
import scipy.special as ss
#===========================================================================


#Select geometry by uncommenting desired option
#Leave remaining two selections commented
Geometry = "slab"
#Geometry = "cylinder"
#Geometry = "sphere"

#Define Diffusion Coefficient
D = 0.6 #[cm]
#Define Absorption Coefficient
SigA = 0.005 #[cm^-1]
#Define slab thickness or radius
T = 45.0 #[cm]
#Define volumetric source strength
So = 1.0E+8 #[n*cm^-1*s^-1]
#Define number of nodes
Nodes = 100
#Define differential thickness
h = T/Nodes
#Define matrix size
N=int(Nodes)
M=N
#Define Matrix A
A=np.zeros( (M,N) )
#Define Matrix b
#This is a 1x1 matrix
b = np.zeros((N,1))

#If the geometry is a slab, do the following
if Geometry == "slab":
    #Set up coefficeints for matrix rows 2 througn N-1
    a_0 = -D/(2.0*h**2.0)
    a_1 = -D/(2.0*h**2.0)
    a_2 = SigA + 2.0*D/(2.0*h**2.0)
    #Input above coefficients into Matrix A  
    for i in range(M - 2):
        A[i+1][i] = 1
        A[i+1][i+1] = -(2.0 + SigA/D*h**2)
        A[i+1][i+2] = 1
    #Input matrix coefficients for x=0 non-reentrant boundary
    A[0][0] = 0.25 + 0.25*D 
    A[0][1] = -D*0.25
    #Input matrix coefficients for x=T non-reentrant boundary    
    A[N-1][N-2] = +D*0.25
    A[N-1][N-1] = 0.25 - 0.25*D
    #Populate Matrix b with source values
    #Recall that the slab had an x^2 dependent source term
    for i in range(N-1):
        b[i][0] = -So/D*h**2.0*(float(i)/float(Nodes)*T)**2.0    
    #Note that the b[0][0] and b[N][0] are  0 due to the boundary conditions 
   
#    #Use Python's solve to determine the flux values
#    y = np.linalg.solve(A,b)
    #Create an x matrix for plotting purposes
    x = np.linspace(0,T,Nodes)
    
    #Create a new matrix, Fie, that is calculated using the TDMA I developed 
    Fie = np.zeros((N,1))
    #Implement the TDMA
    #Forward elimination
    for i in range(1,N):
        m = A[i][i-1]/A[i-1][i-1]
        A[i][i] = A[i][i] - m*A[i-1][i]
        b[i][0] = b[i][0] - m*b[i-1][0]    
    #Backward Substitution    
    Fie[N-1][0] = b[N-1][0]/A[N-1][N-1]
    for i in range(N-2, -1, -1):
        Fie[i][0] = (b[i][0] - A[i][i+1]*Fie[i+1][0]) / A[i][i]
    
    #Analytical Solution
    C1 = -7.02477e+11
    C2 = -3.76233e+12
    X_Analytical = np.linspace(0,45,1000)
    Y_Analytical = C1*np.exp(X_Analytical*(SigA/D)**0.5) + \
        C2*np.exp(-1.0*X_Analytical*(SigA/D)**0.5) + \
        So/D**1.0/( (SigA/D) )*X_Analytical**2 + \
        2.0*So/D*1.0/( (SigA/D)**2 )
   
    #Plot the results
    plt.plot(x, Fie, '--*', color = 'black', label = 'Finite Diff (Via TDMA)')
    plt.plot(X_Analytical, Y_Analytical, '-', color = 'red', label = 'Analytical')   
    plt.xlabel("Slab Thickness [cm]")
    plt.ylabel("Flux [cm^-2*s^-1]")
    plt.title("Slab Flux")
    plt.legend()
    plt.show() 
    
#If the geometry is either a cylinder or sphere
else:    
    T = 45./2.
    h = T/Nodes
    if Geometry == "cylinder":
        c = 1
        PlotTitle = "Cylinder Flux"  
        #Analytical Solution
        B = np.sqrt(SigA/D)
        X_Analytical = np.linspace(0,T,1000)
        Y_Analytical = So/SigA*(1.0 - ss.iv(0,B*X_Analytical)/( ss.iv(0,B*T) + 2.0*D*B*ss.iv(1,B*T) ) )

    if Geometry == "sphere":
        c = 2
        PlotTitle = "Sphere Flux"
        #Analytical Solution
        B = np.sqrt(SigA/D)
        C1 = -T*So/SigA*( 1.0/(np.sinh(B*T) + 2.0*D*B*np.cosh(B*T) - 2.0*D/T*np.sinh(B*T)) )
        X_Analytical = np.linspace(0,T,1000)
        Y_Analytical = C1*np.sinh(B*X_Analytical)/X_Analytical + So/SigA
    #Set up coefficeints for matrix rows 2 througn N-        
    for i in range(1,M - 1):
        A[i][i-1] = -D/(h**2)*( 1.0 - 1.0/(2.0*i) )**c
        A[i][i] = SigA + D/h**2*((1.0 - 1.0/(2.0*i))**c+(1.0 + 1.0/(2.0*i))**c)
        A[i][i+1] = -D/(h**2)*( 1.0 + 1.0/(2.0*i) )**c
    #Input matrix coefficients for 0 slope flux at radius = 0       
    A[0][0] = 1.0
    A[0][1] = -1.0   
    #Input matrix coefficients for x=T non-reentrant boundary       
    A[N-1][N-2] = D*0.25
    A[N-1][N-1] = 0.25 - 0.25*D
    #Populate Matrix b with source values
    #Recall that the slab had an x^2 dependent source term    
    for i in range(N-2):
        b[i+1][0] = So/D
#    #Use Python's solve to determine the flux values
#    y = np.linalg.solve(A,b)
    #Create an x matrix for plotting purposes
    x = np.linspace(0,T,Nodes)
    
    #Create a new matrix, Fie, that is calculated using the TDMA I developed 
    Fie = np.zeros((N,1))
    #Implement the TDMA
    #Forward elimination
    for i in range(1,N):
        m = A[i][i-1]/A[i-1][i-1]
        A[i][i] = A[i][i] - m*A[i-1][i]
        b[i][0] = b[i][0] - m*b[i-1][0]   
    #Backward Substitution    
    Fie[N-1][0] = b[N-1][0]/A[N-1][N-1]
    for i in range(N-2, -1, -1):
        Fie[i][0] = (b[i][0] - A[i][i+1]*Fie[i+1][0]) / A[i][i]    
    #Plot the results
    plt.plot(x, Fie, '--*', color = 'black', label = 'TDMA')
    plt.plot(X_Analytical, Y_Analytical, '-', color = 'red', label = 'Finite Diff (Via TDMA)') 
    plt.xlabel("Radial Distance [cm]")
    plt.ylabel("Flux [cm^-2*s^-1]")
    plt.title(PlotTitle)
    plt.legend()
    plt.show()    
#===========================================================================