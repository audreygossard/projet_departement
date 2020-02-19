import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


a = 0.2
b = 1.3
d_crit = ((1+np.sqrt(2*b*(a+b)))*(a+b)/(a-b))**2

def stabilite(n,d,delta):
    B = np.zeros((2,2))
    pi2 = np.pi*np.pi
    n2 = n*n
    B[0][0] = -n2*pi2+delta*(b-a)/(a+b)
    B[0][1] = delta*(a+b)*(a+b)
    B[1][0] = -2*b*delta/(a+b)
    B[1][1] = -d*n2*pi2-delta*(a+b)*(a+b)
    if np.trace(B)<0 and np.linalg.det(B)>0:
        return True
    return False

couleurs = ['b','g','r','c','m','y']

def diagramme_turing(d1,d2,d_nb,delta1,delta2,delta_nb,n1,n2):
    plt.figure()
    list_d = np.linspace(d1,d2,d_nb)
    list_delta = np.linspace(delta1,delta2,delta_nb)
    for n in range(n1,n2+1):
        d_stables = []
        delta_stables = []
        for d in list_d:
            for delta in list_delta:
                if stabilite(n,d,delta):
                    d_stables.append(d)
                    delta_stables.append(delta)
                    # print("stable!")
                # else:
                    # print("instable!")
        print(couleurs[n])
        plt.plot(d_stables,delta_stables,color = couleurs[n],linestyle = 'None', marker = 'o',alpha = 0.2)
    plt.show()


def diagramme_turing_bis(d1,d2,d_nb,delta1,delta2,delta_nb,n1,n2):
    plt.figure()
    list_d = np.linspace(d1,d2,d_nb)
    list_delta = np.linspace(delta1,delta2,delta_nb)
    for n in range(n1,n2+1):
        mat_d_delta = np.zeros((d_nb,delta_nb))
        for d in list_d:
            for delta in list_delta:
                if stabilite(n,d,delta):
                    mat_d_delta[d][delta] = 1
                    # print("stable!")
                # else:
                    # print("instable!")
        print(couleurs[n])
        plt.imshow(mat_d_delta, alpha = 0.5)
    plt.show()


diagramme_turing(d_crit/100,d_crit*10,100,0.01,100,100,0,2)