import matplotlib.pyplot as plt
import numpy as np
from math import pi
import matplotlib.patches as mpatches
#parameters
a = 0.2
b = 1.3
N = 100
#mode
n = 2
#variations d
d = np.zeros(N)
hd = 5

#variations delta
delta = np.zeros(N)
hdel = 4

#initialisation tableaux
for k in range(N) :
    d[k] = k*hd
    delta[k] = k*hdel

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

def tracer(n) :
    deltatrue = []
    deltafalse = []
    dtrue = []
    dfalse = []
#calcul determinant
    for i in range(N) :
        #listes pour determinant positif ou negatif
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
    plt.xlabel("d")
    plt.ylabel("delta")
    plt.scatter(dtrue, deltatrue, c='lightblue', label = 'stable')
    plt.scatter(dfalse, deltafalse, c='coral', label = 'instable')
    coral_dot = mpatches.Circle((0.5, 0.5), 0.25, facecolor="coral")
    blue_dot = mpatches.Circle((0.5, 0.5), 0.25, facecolor="lightblue")
    plt.legend([coral_dot, blue_dot],["instable", "stable"])
    plt.title("n=" + str(n))
    
    
def courbe(n) :
    dcrit = []
    deltaabs = []
    for i in range(N) :
        dcal = dc(n, delta[i])
        if dcal>0 and dcal<=N*hd :
            dcrit.append(dcal)
            deltaabs.append(delta[i])
    plt.plot(dcrit, deltaabs)

plt.subplot(1,2,1)
tracer(4)
plt.subplot(1,2,2)
tracer(5)

print("end")
plt.savefig('scatterplot.png')
plt.show()
