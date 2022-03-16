from scipy import *
from scipy import interpolate

import scipy.linalg as lin
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from scipy.interpolate import lagrange
from scipy.interpolate import interp1d
 #Size of the List



Dep =[]  
Den =[]  
Nfq =[]  

for line in open('newalldata.txt', 'r'):
  values = [float(s) for s in line.split()]
  Nfq.append(values[0])
  Dep.append(values[1])
  Den.append(values[2])
  

plt.plot(Nfq,Dep, '-')
plt.show()
  
n =[]   
f=[] 


for i in range(len(Dep)-1):    
    n.append((6.67)*Den[i]*(Den[i+1]-Den[i])/(Dep[i+1]-Dep[i]))

   
    
newNfq=[]     
dDen=[]    
newDen=[]
newDep = []
newn = []

print('enter your data dimension : '),
k = input()
d=int(k)
print('enter your path of calculation : ')
l = input()
pa = int(l)
print('enter the number of mode that you would generate : ')
pl = input()
nb=int(pl)



for j in range(0,int(d),int(pa)):        
    newDep.append(Dep[j])
print(newDep)
for i in range(0,int(d),int(pa)):        
    newDen.append(Den[i])


for k in range(0,int(d),int(pa)):        
    newn.append(n[k])
    
for l in range(0,int(d),int(pa)):
    newNfq.append(Nfq[l])
print(newNfq)    
for i in range (0,int(d),int(pa)):
    dDen.append((Den[i+1]-Den[i])/(Dep[i+1]-Dep[i]))
    


#matrix filling    

def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

g=len(newNfq)-1
q=len(newNfq)



a = [-1]*g;
b = [2]*q;
c = [-1]*g;
A1 = tridiag(a, b, c)


a = [1]*g;
b = [0]*q;
c = [-1]*g;
B1 = tridiag(a, b, c)



c = [0]*g;
a = [0]*g;
b=[]
for i in range(0,len(newNfq)):
    b.append(abs(newNfq[i]*newDen[i]))
    
C1 = tridiag(a, b, c)


C=[]
for i in range (q):
    C.append([0]*q)
    for j in range (0,q):
        C[i][j]=newDen[i]
    
Af=1/pa**2*np.multiply(A1,C)

        
D=[]
for i in range (q):
    D.append([0]*q)
    for j in range (0,q):
        D[i][j]=dDen[i]
        
        
Bf=2/pa*np.multiply(B1,D)  




#Calculate the Eigenmodes
Efu=[]
Efu=lin.eigh(Af+Bf,C1,type=1)[1]



#Calculate the Eigenvalues
Eva=[]
Eva=lin.eigh(Af+Bf,C1,type=1)[0]


#plt.plot(newDep,newDen, "-")
#plt.show()
newEva=[]
for i in range(0,len(Eva)):
    newEva.append(np.sqrt(1/Eva[i]))

print(newEva)    
x1=np.linspace(1,len(newEva),len(newEva))

plt.plot(x1, newEva, 'o')
plt.show()
    

#Using the Chebychev Nodes
mnewDep=[]
for i in range(len(newDep)):
    mnewDep.append(-newDep[i])

for i in range(0,nb):
    plt.plot(mnewDep, 10*Efu[i], '-')
plt.show()




for i in range(0,nb):
    x=np.asarray(newDep)
    y=np.asarray(Efu[i])

#  Calculate the polynomial coefficients
    L = lagrange(x, y)
    L = np.polynomial.polynomial.Polynomial(L).coef


#  plot the polynomial
    X = np.linspace(1,len(y), 10000)
  

    x = np.arange(1,len(y)+1 )
    n = x.size
    x = np.cos((2*x-1) / (2 * n) * np.pi )



#  Polynomial coefficients using the new nodes
    L = lagrange(x, y)
    L = np.polynomial.polynomial.Polynomial(L).coef




    X = np.linspace(0, 1, 100)

    plt.plot(-X*260,np.polyval(L, X))
   
    plt.grid(True)
   
  
    coefs = [0,] * (len(y)+1)
    coefs[-1] = 1

#  Make the polynomial object.
    C = np.polynomial.chebyshev.Chebyshev(coefs)

#  Calculate the roots
    R = C.roots()

#  Make sure our nodes from the cell above are in ascending order
    x = np.sort(x)

#  Print both the roots and the nodes from above out together.
    for r, u in zip(R, x):
        print(r, u)

#  Are they the same?
    print( np.allclose(x, R)) 
plt.show()      






