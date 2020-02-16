import matplotlib.pyplot as plt
import numpy as np
from math import pi
#parameters
a = 0.2
b = 1.3
N = 100
#mode
n = 2
#variations d
d = np.zeros(N)
hd = 10

#variations delta
delta = np.zeros(N)
hdel = 0.1

#initialisation tableaux
for k in range(N) :
    d[k] = k*hd
    delta[k] = k*hdel+50

det = np.zeros((N,N))

def determinant(mode, dpar, deltapar) :
    a1 =((mode*pi)**4)*dpar*(a+b)
    b1 = ((mode*pi)**2)*(dpar*(a-b)+(a+b)**3)*deltapar
    c = (deltapar**2)*((a+b)**3)
    return (a1+b1+c)


#d_critique
def dc(n, delta) :
    beta = (n*pi)**2*delta*(a+b)**3+delta**2*(a+b)**3
    alpha = (n*pi)**4*(a+b)+(n*pi)**2*delta*(a-b)
    return (-beta/alpha)

def deltac(n) :
    return ((n*pi)**2)*(a+b)/(b-a)


#calcul determinant
for i in range(N) :
    #listes pour determinant positif ou negatif
    deltatrue = []
    deltafalse = []
    dtrue = []
    dfalse = []
    for j in range(N) :
        det[i][j] = determinant(n,d[i], delta[j])
        #si determinant positif on ajoute d et delta aux listes 'true'
        if (det[i][j]>0):
            dtrue.append(d[i])
            deltatrue.append(delta[j])
        #sinon aux listes 'false'
        else :
            dfalse.append(d[i])
            deltafalse.append(delta[j])
    plt.scatter(dtrue, deltatrue, c='lightblue')
    plt.scatter(dfalse, deltafalse, c='coral')
print(dtrue, deltatrue)
print("end")
plt.savefig('scatterplot.png')
plt.show()
